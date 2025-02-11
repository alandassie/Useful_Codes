"""
    Created on February 2025 by Alan D.K. for Hypernuclei project

    This code run GSM and then optimize the one-hyperon or
    nucleon-hyperon interaction in order to reproduce the experimental data. 
    You have to define the input and output files and put
    the experimental energies.
    
    An example of the input file "GMSHyper_Opt.in" can be found
    at the end of the code.
"""

import subprocess as sp
import time
import os
import math as m
# from scipy.optimize import least_squares as least
from scipy.optimize import newton
from scipy.optimize import minimize
import numpy as np

# LOG FILE
logfile = os.getcwd() + '/' + 'GSMHyper_Opt_log.' + time.strftime( "%y.%m.%d-%H.%M", time.localtime() )
# LOG NAME
logname = 'WD_GSMHyper_Opt_log.' + time.strftime( "%y.%m.%d-%H.%M", time.localtime() )
# LOG FOLDER
logfolder = os.getcwd() + '/' + logname 
if not os.path.exists(logfolder):
    os.makedirs(logfolder)
else:
    for filename in os.listdir(logfolder):
        # remove the files inside
        os.remove(f"{logfolder}/{filename}")

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
def f(x):
    # Open GSM input file
    with open(readfilename_GSM,'r') as gsmin:
        inputfile_lines = gsmin.read().split('\n')
    # 
    start_value = 0
    if opt_onehyperon_ws == 1: # Optimizing WS part
        print_twice('Optimizing 1Y WS part')
        start_value = len(onehyperon_ws_l)
        for i in range(0,hyperon_names):
            onehyperon_ws_line = onehyperon_ws_lines[i]
            for j in range(0,len(onehyperon_ws_l)):
                l = onehyperon_ws_l[j]
                aux1 = inputfile_lines[onehyperon_ws_line+l].split()
                aux2 = str(x[j])
                # 
                inputfile_lines[onehyperon_ws_line+l] = "    " + str(l) + "  " + aux1[1] + "  " + aux1[2] + aux2 + "  " + aux1[4]
    if opt_yn == 1: # Optimizing YN interactions
        print_twice('Optimizing YN interactions')
        for i in range(0,yn_n):
            yn_line = yn_lines[i]
            aux1 = inputfile_lines[yn_line].split()
            aux2 = str(x[start_value+i])
            inputfile_lines[yn_line] = "  " + aux2 + "  " + aux1[1]
    # Print input array
    print_twice('\n All the interaction strenghts:')
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
    line = line_numbers[1] + 3
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
        auxjpi = outputfile_lines[line+j].split()[3].split('(')[0]
        auxindex = outputfile_lines[line+j].split()[3].split('(')[1].split(')')[0]
        aux_state = '%s(%s)'% (auxjpi,auxindex)
        if aux_state in state_read:
            state.append(aux_state)
            jpi.append(auxjpi)
            index.append(auxindex)
            auxenergy = outputfile_lines[line+j].split()[4].split(':')[1]
            energy.appen(float(auxenergy))
        j += 1
    #
    auxiliar = list(zip(state,energy))
    print_twice('Calculated energies of the optimized states:')
    print_twice(auxiliar)
    # Compare with all the experimental energies
    res_e = np.zeros(numberofstates)
    for i in range(0,numberofstates):
        expstate = state_read[i]
        expene = expene_read[i]
        for j in range(0,numberofstates):
            if state[j] == expstate:
                res_e[i] = ( expene - energy[j] )**2 / abs(expene) * stweig_read[i]
                print_twice('State {0:s}, E Residue = {2:10.6f}'.format(expstate,res_e[i]))
    res = np.sum(res_e)
    print_twice('X^2 = {0:10.6f}'.format(res))
    #
    # Check which optimizator we are using
    if method == 'NEWTON':
        print_twice('Finish iteration of Newton optimizer, x and f(x) must be of the same size!')
        print_twice('\n'+20*'-'+'\n')
        return res_e
    elif method == 'MINIMIZATION':
        print_twice('Finish iteration of TNC or Nelder-Mead optimizer')
        print_twice('\n'+20*'-'+'\n')
        return res
    else:
        print_twice('METHOD must be NEWTON or MINIMIZATION;(TNC or Nelder-Mead)')
        exit()
# .-

# INPUT FILE
readfilename = "GSMHyper_Opt.in"
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
used_icf = 0
icf_type = 'NONE'
icf_bounds = ''
theline = searchline(readfilename,"OPTIMIZEDHYPERONS:")
hyeron_n = int(data[theline+1])
hyperon_names = []
for i in range(0,hyeron_n):
    hyperon_names.append(data[theline+2+i])
#
# Checking if one-hyperon WS interaction will be optimized
opt_onehyperon_ws = 0
onehyperon_ws_bounds = ''
theline = searchline(readfilename,"VWSOPTIMIZATION:")
if theline != None:
    opt_onehyperon_ws = 1
    onehyperon_ws_npw = int(data[theline+1])
    onehyperon_ws_l = onehyperon_ws_npw*[0]
    onehyperon_ws_seed = onehyperon_ws_npw*[0]
    #
    for i in range(0,onehyperon_ws_npw):
        factor = i*3
        onehyperon_ws_l[i] = int(data[theline+2+factor])
        onehyperon_ws_seed[i] = float(data[theline+3+factor])
        if i == 0:
            onehyperon_ws_bounds += data[theline+4+factor]
        else:
            onehyperon_ws_bounds += ',' + data[theline+4+factor]
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
numberofstates = int(data[theline+1])
jpi_read = []
expene_read = np.zeros(numberofstates)
stweig_read = np.zeros(numberofstates)
for i in range(0,numberofstates):
    factor = i*2
    jpi_read.append(data[theline+2+factor])
    expene_read[i] = float(data[theline+3+factor])
    stweig_read[i] = float(data[theline+4+factor])
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
# Defining optimization part
search = 'Spectrum'
# Lines with the 1Y WS part
if opt_onehyperon_ws == 1:
    theline = searchline(readfilename_GSM,"core.potential")
    onehyperon_ws_lines = []
    for i in range(0,hyeron_n):
        shift = [x.strip(' ') for x in inputfile_lines[theline:theline+30]].index(hyperon_names[i]) + 2
        onehyperon_ws_lines.append(theline+shift)
# Lines with the YN interactions
if opt_yn == 1:
    theline = searchline(readfilename_GSM,"Hamiltonian.interaction")
    yn_lines = []
    for i in range(0,yn_n):
        shift = [x.strip(' ') for x in inputfile_lines[theline:theline+20]].index(yn_names[i])
        yn_lines.append(theline+shift)
#
#
# Start calculation
start_main = time.time()
print_twice("\nRunning GSM in %s"% gsm_directory)
os.chdir(gsm_directory)
seeds = []
if opt_onehyperon_ws == 1:
    seeds += onehyperon_ws_seed
if opt_yn == 1:
    seeds += yn_seed
#
if method == 'NEWTON':
    print_twice('Using Newton optimizer, x and f(x) must be of the same size!')
    opt = newton(f, seeds, tol=1e-10, maxiter=20, full_output=True)
elif method == 'MINIMIZATION':
    print_twice('Using %s optimizer'% mini_method)
    if mini_method == 'TNC':
        opt = minimize(f, seeds, method=mini_method, jac='2-point', options={ 'xtol' : 1e-3, 'finite_diff_rel_step': 0.001 }) # idea from https://stackoverflow.com/questions/20478949/how-to-force-larger-steps-on-scipy-optimize-functions
    elif mini_method == 'Nelder-Mead':
        if opt_onehyperon_ws == 1:
            bounds_opt = (eval(onehyperon_ws_bounds + ',' + yn_bounds))
        else:
            bounds_opt = ( eval(yn_bounds + ',') )
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
    Example of GSMHyper_Opt.in file:
    _________________________________
    STORAGE-DIRECTORY:
    /pbs/throng/ganil/adassie/Carbon-11_Project/storage_4MPI-24.11.20_GSMOpt-24.08.26-11.00_Basis-24.12.10-13.00
        
    PARALLELISM: 1 - MPI or OPENMP; 2 - NUMBER OF NODES
    MPI
    4

    MACHINEFILE: 1 - NAME OF THE FILE
    machinefile

    OPTIMIZATIONMETHOD:  NEWTON or ('MINIMIZATION;' + 'TNC' or 'Nelder-Mead' for the moment)
    MINIMIZATION;Nelder-Mead

    OPTIMIZEDHYPERONS: 1 - NUMBER OF HYPERONS TO OPTIMIZE; 2,N - NAME OF EACH HYPERON
    Lambda
    Sigma0
    Sigma-

    VWSOPTIMIZATION: 1 - REAL SEED; 2 - IF Nelder-Mead, DEFINE BOUNDS IN A FORM (MIN,MAX)
    40
    (30,60)

    YNOPTIMIZATION: for each YN interaction -> 1 - NAME; 2 - SEED; 3 - IF Nelder-Mead, DEFINE BOUNDS IN A FORM (MIN,MAX)
    V8a.SU3.f
    -0.2
    (-1,1)
    V8s.SU3.f
    -0.25
    (-1,1)

    EXPERIMENTALVALUES: 1 - NUMBER OF EXPERIMENTAL STATES; for each state -> 2 - JPi; 3 - ENERGY (MeV); 4 - WEIGHT BETWEEN (0,1)
    2
    0+
    -36.446
    1
    2+
    -35.945
    1

    GSM-exe:
    GSM-24.11.20-MPI.x
    GSM-file: 1-input file; 2-output file; 3-1 for overwrite or 2 for append
    IN2P3-HTC-25.01.11_40CaLambda-COSM_GSMOpt-25.01.20-20.00.in
    IN2P3-HTC-25.01.11_40CaLambda-COSM_GSMOpt-25.01.20-20.00.out
    2
    _________________________________
"""