"""
    Created on February 2025 by Alan D.K. for Hypernuclei project

    This code run GSM and then optimize the one-hyperon or
    nucleon-hyperon interaction in order to reproduce the experimental data. 
    You have to define the input and output files and put
    the experimental energies.
    
    An example of the input file "input.GMSHyper_Opt" can be found
    at the end of the code.
"""

import subprocess as sp
import time
import os
import re
import math as m
# from scipy.optimize import least_squares as least
from scipy.optimize import newton
from scipy.optimize import minimize
import numpy as np

# LOG NAME
logname = 'log.WD_GSMHyper_Opt.' + time.strftime( "%y.%m.%d-%H.%M", time.localtime() )
# LOG FOLDER
logfolder = os.getcwd() + '/' + logname 
if not os.path.exists(logfolder):
    os.makedirs(logfolder)
else:
    for filename in os.listdir(logfolder):
        # remove the files inside
        os.remove(f"{logfolder}/{filename}")
# LOG FILE
logfile = os.getcwd() + '/' + logname + '/' + 'log.GSMHyper_Opt.' + time.strftime( "%y.%m.%d-%H.%M", time.localtime() )

# Declaration of funcitons
def erease_output_file():
    with open(logfile,'w') as output:
        pass
def print_twice(*args,**kwargs): # Allow to print in file and in terminal at the same line
    print(*args,**kwargs)
    with open(logfile,"a") as f:  # appends to file and closes it when finished
        print(file=f,*args,**kwargs)
# .-
# Idea from https://stackoverflow.com/questions/3961265/get-line-number-of-certain-phrase-in-text-file
def searchline(file,phrase):
    with open(file, 'r', encoding='utf-8') as f:
        text = f.read()
        if phrase in text:
            phrase_index = text.index(phrase)
            l_num = text[:phrase_index].count('\n')  # Nth line has n-1 '\n's before
        else:
            l_num = None
    return l_num
# Idea from https://stackoverflow.com/questions/3873361/finding-multiple-occurrences-of-a-string-within-a-string-in-python
def searchline_all(file,phrase):
    with open(file, 'r', encoding='utf-8') as f:
        text = f.read()
        count = 0
        i = len(phrase)
        x = 0
        l_num = []
        while x < len(text) - (i-1):
            if text[x:x+i] == phrase:
                l_num.append(text[:x].count('\n'))  # Nth line has n-1 '\n's before
                x += i
                count += 1
            else:
                x += 1
    return l_num
#
def searchlinefinal(file,phrase):
    with open(file, 'r', encoding='utf-8') as f:
        text = f.read()
        if phrase in text:
            phrase_index = text.rindex(phrase)
            l_num = text[:phrase_index].count('\n')  # Nth line has n-1 '\n's before
        else:
            l_num = None
    return l_num
# .-
def extract_e_value_as_float(s):
    match = re.search(r'E:([\d.-]+)\s*MeV', s)
    if match:
        return float(match.group(1))
    return None
# .-
def extract_spin_parity_and_parenthesis(s):
    # Find all occurrences of patterns like '1/2+(0)', '0+(0)', etc.
    matches = re.findall(r'([\d/+-]+)\s*\((\d+)\)', s)
    return [(spin_parity, int(parenthesis)) for spin_parity, parenthesis in matches]
# .-
def f(x):
    # 
    # Open GSM input file
    with open(readfilename_GSM,'r') as gsmin:
        inputfile_lines = gsmin.read().split('\n')
    #
    start_value = 0
    if opt_onebaryon == 1: # Optimizing one baryon part
        print_twice('Optimizing 1 baryon part')
        # start_value = len(onebaryon_l)
        for i in range(0,baryon_n):
            onebaryon_line = onebaryon_lines[i]
            aux_onebaryon_l = onebaryon_l[i]
            aux_onebaryon_samev0 = onebaryon_samev0[i]
            for j in range(0,len(aux_onebaryon_l)):
                l = aux_onebaryon_l[j]
                aux1 = inputfile_lines_start[onebaryon_line+l].split()
                aux3 = float(aux1[onebaryon_shift])
                if same_corrective_factor.upper() == 'YES':
                    aux2 = str(x[0]*aux3)
                elif aux_onebaryon_samev0.upper() == 'YES':
                    aux2 = str(x[start_value]*aux3)
                else: aux2 = str(x[start_value+j]*aux3)
                #  
                if onebaryon_type[0] == 'WS':
                    inputfile_lines[onebaryon_line+l] = "    " + str(l) + "   " + aux1[1] + "   " + aux1[2] + "   " + aux2 + "   " + aux1[4]
                elif onebaryon_type[0] == 'SO':
                    inputfile_lines[onebaryon_line+l] = "    " + str(l) + "   " + aux1[1] + "   " + aux1[2] + "   " + aux1[3] + "   " + aux2
            if aux_onebaryon_samev0.upper() == 'NO':
                start_value += len(aux_onebaryon_l)
            else:
                start_value += 1
    if opt_yn == 1: # Optimizing YN interactions
        print_twice('Optimizing YN interactions')
        for i in range(0,yn_n):
            yn_line = yn_lines[i]
            aux1 = inputfile_lines_start[yn_line].split()
            aux2 = str(x[start_value+i]*float(aux1[0]))
            inputfile_lines[yn_line] = "  " + aux2 + "  " + aux1[1]
    # Print input array
    print_twice('\n All the interaction corrective factors strenghts:')
    print_twice(x)
    # Save and close GSM input file
    inputfile_aux = '\n'.join(inputfile_lines)
    with open(readfilename_GSM,'w') as gsmin:
        gsmin.write(inputfile_aux)
    #
    aux = time.strftime( "%y.%m.%d-%H.%M.%S", time.localtime() )
    outfilename_GSMi = logname + '/' + outfilename_GSM + '-' + aux
    start_gsm = time.time()
    print_twice('\n ' + running_prefix + running_gsm + ' < ' + readfilename_GSM + gsm_write + outfilename_GSMi)
    sp.run([running_prefix + running_gsm + ' < ' + readfilename_GSM + gsm_write + outfilename_GSMi], shell=True)
    end_gsm = time.time()
    time_gsm = end_gsm-start_gsm
    print_twice("Time to calculate: ",time_gsm, "s")
    #
    # Read results
    with open(outfilename_GSMi,'r') as gsmout:
        outputfile_lines = gsmout.read().split('\n')
    line_numbers = searchline_all(outfilename_GSMi,search)
    line = line_numbers[0] + 3
    state = []
    jpi = []
    index = []
    energy = []
    i = 0
    j = 0
    while i == 0:
        aux = outputfile_lines[line+j].split()
        if aux == []:
            i = 1
            continue
        # auxjpi = outputfile_lines[line+j].split()[3].split('(')[0]
        # auxindex = outputfile_lines[line+j].split()[3].split('(')[1].split(')')[0]
        auxjpi, auxindex = extract_spin_parity_and_parenthesis(outputfile_lines[line+j])[0]
        aux_state = '%s(%s)'% (auxjpi,auxindex)
        print(auxjpi,auxindex,aux_state,state_read)
        if aux_state in state_read:
            state.append(aux_state)
            jpi.append(auxjpi)
            index.append(auxindex)
            auxenergy = extract_e_value_as_float(outputfile_lines[line+j])
            energy.append(float(auxenergy))
        j += 1
    #
    auxiliar = list(zip(state,energy))
    print_twice('Calculated energies of the optimized states:')
    print_twice(auxiliar)
    # Compare with all the experimental energies
    if opt_separation_energy == 0:
        res_e = np.zeros(numberofstates)
        for i in range(0,numberofstates):
            expstate = state_read[i]
            expene = expene_read[i]
            for j in range(0,numberofstates):
                if state[j] == expstate:
                    res_e[i] = ( expene - energy[j] )**2 / abs(expene) * stweig_read[i]
                    print_twice('State {0:s}, E Residue = {1:10.6f}'.format(expstate,res_e[i]))
        res = np.sum(res_e)
        print_twice('X^2 = {0:10.6f}'.format(res))
    if opt_separation_energy == 1:
        res_e = np.zeros(numberofpairs)
        for i in range(0,numberofpairs):
            expstate1 = pair_state1[i]
            expstate2 = pair_state2[i]
            expsepene = pair_separation_energy[i]
            gsmstate1 = state.index(expstate1)
            gsmstate2 = state.index(expstate2)
            gsmsepene = abs(energy[gsmstate1] - energy[gsmstate2])
            res_e[i] = ( expsepene - gsmsepene )**2 / abs(expsepene)
            print_twice('States {0:s} and {1:s}, Separation Energy Residue = {2:10.6f}'.format(expstate1,expstate2,res_e[i]))
        res = np.sum(res_e)
        print_twice('X^2 = {0:10.6f}'.format(res))
    #
    # Check which optimizator we are using
    if method == 'NEWTON':
        print_twice('Finish iteration of Newton optimizer, x and f(x) must be of the same size!')
        print_twice('\n'+20*'-'+'\n')
        return res_e
    elif method == 'MINIMIZATION':
        print_twice('Finish iteration of MINIMIZATION optimizer')
        print_twice('\n'+20*'-'+'\n')
        return res
    else:
        print_twice('METHOD must be NEWTON or MINIMIZATION;(TNC or Nelder-Mead or BFGS)')
        exit()
# .-

# INPUT FILE
readfilename = "input.GSMHyper_Opt"
with open(readfilename, 'r') as readfile:
    data = readfile.read().split('\n')
# LOGILE
erease_output_file()
    
print_twice("Be sure that you are using the correct directories!")
# GSM directory
gsm_directory = os.getcwd()
# sure = input("Is %s the GSM working directory?\n YES or NO: "% gsm_directory)
# if sure.lower() == "no":
#     print("Change it in python file!")
#     exit()
# storage directory
theline = searchline(readfilename,"STORAGE-DIRECTORY:")
storage_directory = data[theline+1]
# sure = input("Is %s the storage directory?\n YES or NO: "% storage_directory)
# if sure.lower() == "no":
#     print("Change it in python file!")
#     exit()
# Checking if it is a MPI or OPENMP/secuential calculation
theline = searchline(readfilename,"PARALLELISM:")
parallelism_type = data[theline+1]
parallelism_nodes = data[theline+2]
if parallelism_type == 'MPI':
    running_prefix = 'mpirun -np ' + parallelism_nodes + ' '
elif parallelism_type == 'OPENMP':
    running_prefix = ' ./'
else:
    print_twice('Parallelism must be MPI or OPENMP')
    exit()
# Checking if we need machinefile
theline = searchline(readfilename,"MACHINEFILE:")
if theline != None:  
    machinefile_name = data[theline+1]
    running_prefix = running_prefix + '-hostfile ' + machinefile_name + ' '
#
# Looking up the optimization code to use
theline = searchline(readfilename,"OPTIMIZATIONMETHOD:")
method =  data[theline+1].split(';')[0]
if method == 'MINIMIZATION':
    mini_method = data[theline+1].split(';')[1]
#
# Look for optimized hyperons
opt_onebaryon = 0
theline = searchline(readfilename,"OPTIMIZEDBARYONS:")
if theline != None:
    print_twice('One baryon interaction will be optimized') 
    opt_onebaryon = 1
    baryon_n = int(data[theline+1])
    same_corrective_factor = data[theline+2]
    # If we are going to use the same corrective factor for all baryons, also the same for all partial waves
    if same_corrective_factor.upper() == 'YES':
        seeds_size = 1
        print_twice('Same corrective factor for all optimized one-baryon interactions')
    baryon_names = []
    for i in range(0,baryon_n):
        baryon_names.append(data[theline+3+i])
#
# Checking if one-baryon WS interaction will be optimized
for k in range(0,baryon_n):
    theline = searchline(readfilename, "OPT_" + baryon_names[k].upper() + ":")
    if theline != None:
        aux_onebaryon_type = data[theline+1] # Can be WS, SO for the moment
        if aux_onebaryon_type == 'R0' or aux_onebaryon_type == 'AA':
            print_twice('Only WS or SO fit for now')
            exit()
        aux_onebaryon_npw = int(data[theline+2])
        aux_onebaryon_l = aux_onebaryon_npw*[0]
        aux_onebaryon_samev0 = data[theline+3]
        print_twice('One-%s %s interaction of %s partial waves will be optimized'% (baryon_names[k],aux_onebaryon_type,aux_onebaryon_l))
        if aux_onebaryon_samev0.upper() == 'YES':
            aux_onebaryon_seed_n = 1
        else:
            aux_onebaryon_seed_n = aux_onebaryon_npw
        if aux_onebaryon_samev0.upper() == 'NO' and same_corrective_factor.upper() == 'YES':
            print_twice('If the same corrective factor is used for all baryons, the same should be used for all WS partial waves')
            exit()
        #
        aux_onebaryon_seed = aux_onebaryon_seed_n*[0]
        for i in range(0,aux_onebaryon_npw):
            factor = i*3
            aux_onebaryon_l[i] = int(data[theline+4+factor])
            if aux_onebaryon_samev0.upper() == 'YES' and i == 0:
                aux_onebaryon_seed[0] = float(data[theline+5+factor])
                aux_onebaryon_bounds = data[theline+6+factor]
            elif aux_onebaryon_samev0.upper() == 'NO':
                aux_onebaryon_seed[i] = float(data[theline+5+factor])
                if i == 0:
                    aux_onebaryon_bounds = data[theline+6+factor]
                else:
                    aux_onebaryon_bounds += ',' + data[theline+6+factor]
        print_twice('    Seeds: %s; Bounds: %s'% (aux_onebaryon_seed,aux_onebaryon_bounds))
        #
        if k == 0:
            onebaryon_type = [aux_onebaryon_type]
            onebaryon_l = [aux_onebaryon_l]
            onebaryon_samev0 = [aux_onebaryon_samev0]
            onebaryon_seed = aux_onebaryon_seed
            onebaryon_bounds = aux_onebaryon_bounds
        else:
            onebaryon_type.append(aux_onebaryon_type)
            onebaryon_l.append(aux_onebaryon_l)
            onebaryon_samev0.append(aux_onebaryon_samev0)
            if same_corrective_factor.upper() == 'NO':
                onebaryon_seed += aux_onebaryon_seed
                onebaryon_bounds += ',' + aux_onebaryon_bounds
# Check that the same part of the one body will be optimized
if opt_onebaryon == 1:
    aux = all(x == onebaryon_type[0] for x in onebaryon_type)
    if aux == False:
        print_twice('For OB, the same part of the interaction should be optimized on every baryon')
        exit()
    if onebaryon_type[0] == 'WS':
        onebaryon_shift = 3
    if onebaryon_type[0] == 'SO':
        onebaryon_shift = 4
#
# Checkin if YN interactions will be optimized
opt_yn = 0
yn_bounds = ''
theline = searchline(readfilename,"YNOPTIMIZATION:")
if theline != None:
    opt_yn = 1
    yn_n = int(data[theline+1])
    yn_names = []
    yn_seed = yn_n*[0]
    for i in range(0,yn_n):
        factor = i*3
        yn_names.append(data[theline+2+factor])
        yn_seed[i] = float(data[theline+3+factor])
        if i == 0:
            yn_bounds += data[theline+4+factor]
        else:
            yn_bounds += ',' + data[theline+4+factor]
#
# Reading experimental data
theline = searchline(readfilename,"EXPERIMENTALVALUES:")
print_twice('Experimental data:')
numberofstates = int(data[theline+1])
state_read = []
jpi_read = []
index_read = []
expene_read = np.zeros(numberofstates)
stweig_read = np.zeros(numberofstates)
for i in range(0,numberofstates):
    factor = i*4
    jpi_read.append(data[theline+2+factor].split()[0])
    index_read.append(data[theline+3+factor].split()[0])
    state_read.append( '%s(%s)'% (jpi_read[i],index_read[i]) )
    expene_read[i] = float(data[theline+4+factor])
    stweig_read[i] = float(data[theline+5+factor])
    print_twice('%s(%s) : E = %10.6f, Weight = %3.1f'% (jpi_read[i],index_read[i],expene_read[i],stweig_read[i]) )
# If experimental values are even, check if the separation energy should be optimized instead of the energy
opt_separation_energy = 0
if numberofstates%2 == 0:
    theline = searchline(readfilename,"SEPARATIONENERGY:")
    if theline != None:
        print_twice('Optimizing separation energy instead of energy')
        opt_separation_energy = 1
        numberofpairs = int(data[theline+1])
        pair_state1 = []
        pair_state2 = []
        pair_separation_energy = np.zeros(numberofpairs)
        for i in range(0,numberofpairs):
            factor = i*2
            aux_jpi = data[theline+2+factor].split(',')[0]
            aux_index = data[theline+3+factor].split(',')[0]
            pair_state1.append('%s(%s)'% (aux_jpi,aux_index) )
            print_twice('Pair 1: %s(%s)'% (aux_jpi,aux_index))
            aux_jpi = data[theline+2+factor].split(',')[1]
            aux_index = data[theline+3+factor].split(',')[1]
            pair_state2.append('%s(%s)'% (aux_jpi,aux_index) )
            print_twice('Pair 2: %s(%s)'% (aux_jpi,aux_index))
            #
            ene_1 = expene_read[state_read.index(pair_state1[i])]
            ene_2 = expene_read[state_read.index(pair_state2[i])]
            pair_separation_energy[i] = abs(ene_1 - ene_2)
            print_twice('Pair experimental separation energy: %10.6f', pair_separation_energy[i])
#
# Executable GSM file
theline = searchline(readfilename,"GSM-exe:")
if theline != None:  
    running_gsm = data[theline+1]
# Read-Out file name GSM
theline = searchline(readfilename,"GSM-files:")
readfilename_GSM = data[theline+1]
outfilename_GSM = data[theline+2]
gsm_write_aux = int(data[theline+3])
gsm_write = ' ' + gsm_write_aux*'>' + ' '
#
# Open GSM input file
with open(readfilename_GSM,'r') as gsmin:
    inputfile_lines_start = gsmin.read().split('\n')
#
# Defining optimization part
search = 'Spectrum'
# Lines with the 1 baryon WS part
if opt_onebaryon == 1:
    theline = searchline(readfilename_GSM,"core.potential")
    onebaryon_lines = []
    for i in range(0,baryon_n):
        shift = [x.strip(' ') for x in inputfile_lines_start[theline:theline+30]].index(baryon_names[i]) + 2
        onebaryon_lines.append(theline+shift)
# Lines with the YN interactions
if opt_yn == 1:
    theline = searchline(readfilename_GSM,"Hamiltonian.interaction")
    aux = [x.strip(' ') for x in inputfile_lines_start[theline:theline+20]]
    yn_lines = []
    for i in range(0,yn_n):
        for j in range(20):
            if yn_names[i] in aux[j]:
                shift = j
                break
        yn_lines.append(theline+shift)
#
#
# Start calculation
start_main = time.time()
print_twice("\nRunning GSM in %s"% gsm_directory)
os.chdir(gsm_directory)
seeds = []
if opt_onebaryon == 1:
    # Check sizes
    if same_corrective_factor.upper() == 'YES':
        if len(onebaryon_seed) > 1:
            print_twice('The seed size should be one if you are going to use the same corrective factor for all baryons!')
            exit()
    seeds += onebaryon_seed
if opt_yn == 1:
    seeds += yn_seed
#
if method == 'NEWTON':
    print_twice('Using Newton optimizer, x and f(x) must be of the same size!')
    opt = newton(f, seeds, tol=1e-10, maxiter=20, full_output=True)
elif method == 'MINIMIZATION':
    print_twice('Using %s optimizer'% mini_method)
    if opt_onebaryon == 1:
        bounds_opt = (eval(onebaryon_bounds + ',' + yn_bounds))
    else:
        bounds_opt = ( eval(yn_bounds + ',') )
    if mini_method == 'TNC':
        opt = minimize(f, seeds, method=mini_method, bounds=bounds_opt, jac='2-point', options={ 'xtol' : 1e-3, 'finite_diff_rel_step': 0.001 }) # idea from https://stackoverflow.com/questions/20478949/how-to-force-larger-steps-on-scipy-optimize-functions
    elif mini_method == 'Nelder-Mead' or mini_method == 'L-BFGS-B':
        opt = minimize(f, seeds, method=mini_method, bounds=bounds_opt)
else:
    print_twice('METHOD must be NEWTON or MINIMIZATION!')
    exit()
#
print_twice(opt)
#
end_main = time.time()
time_main = end_main-start_main
print_twice("\n\nAll calculations lasted: ", time_main, "s")


"""
    Example of input.GSMHyper_Opt file:
    _________________________________
    STORAGE-DIRECTORY:
    /pbs/throng/ganil/adassie/Carbon-11_Project/storage_4MPI-24.11.20_GSMOpt-24.08.26-11.00_Basis-24.12.10-13.00
        
    PARALLELISM: 1 - MPI or OPENMP; 2 - NUMBER OF NODES
    MPI
    4

    MACHINEFILE: 1 - NAME OF THE FILE
    machinefile

    OPTIMIZATIONMETHOD:  'NEWTON' or 'MINIMIZATION;' + ('TNC' or 'Nelder-Mead' or 'L-BFGS-B')
    MINIMIZATION;Nelder-Mead

    OPTIMIZEDBARYONS: 1 - NUMBER OF BARYONS TO OPTIMIZE; 2 - SAME CORR FACTOR FOR ALL BARYONS; 2,N - NAME OF EACH BARYON
    5
    YES
    neutron
    proton
    Lambda
    Sigma0
    Sigma-

    OPT_NEUTRON:  1 - OPTIMZE WS, SO, R0, AA;2 - NUMBER OF PARTIAL WAVES; 3 - SAME V0 FOR ALL PW; 4 - \ell PW; 5 - REAL SEED CF; 6 - DEFINE BOUNDS IN A FORM (MIN,MAX) (Works with MINIMIZATION)
    WS
    1
    YES
    0
    1
    (0.9,1.2)

    OPT_LAMBDA:  1 - OPTIMZE WS, SO, R0, AA;2 - NUMBER OF PARTIAL WAVES; 3 - SAME V0 FOR ALL PW; 4 - \ell PW; 5 - REAL SEED CF; 6 - DEFINE BOUNDS IN A FORM (MIN,MAX) (Works with MINIMIZATION)
    WS
    1
    YES
    0
    1
    (0.9,1.2)

    YNOPTIMIZATION: 1- Number of YN interaction; for each YN interaction -> 2 - NAME; 3 - SEED; 4 - DEFINE BOUNDS IN A FORM (MIN,MAX)
    2
    V8a.SU3.f
    -0.2
    (-1,1)
    V8s.SU3.f
    -0.25
    (-1,1)

    EXPERIMENTALVALUES: 1 - NUMBER OF EXPERIMENTAL STATES; for each state -> 2 - JPi; 3 - index; 4 - ENERGY (MeV); 5 - WEIGHT BETWEEN (0,1)
    2
    0+
    0
    -134.500
    1
    2+
    0
    -135.705
    1
    
    SEPARATIONENERGY: 1 - Number of pairs of states to optimize the separation energy instead of the energy; for each pair -> 2 - JPi_1,Jpi_2; 3 - index_1,index_2
    1
    0+,2+
    0,0

    GSM-exe:
    GSM-24.11.20-MPI.x
    GSM-files: 1-input file; 2-output file; 3-1 for overwrite or 2 for append
    IN2P3-HTC-25.01.11_40CaLambda-COSM_GSMOpt-25.01.20-20.00.in
    IN2P3-HTC-25.01.11_40CaLambda-COSM_GSMOpt-25.01.20-20.00.out
    2
    _________________________________
"""