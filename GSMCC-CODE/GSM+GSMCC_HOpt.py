"""
    Created on May 2025 by Alan D.K. for 11C_Project

    This code run GSM (if needed) and GSMCC and then
    optimize the One- and/or Two-body Hamiltonian interaction
    so as to reproduce the ENERGY of one state, or the SEP ENERGY
    between two of them.
    
    For reading file, the information about working directories and files
    is the same as in the code GSM+GSMCC_run.py.
    You have to define the experimental values
    and the state/s to optimize in "input.GMS+GSMCC_HOpt". 
    An example can be found at the end of the code.
    
    For the moment, this will only modify GSMCC input file.
"""

import subprocess as sp
import time
import os
import math as m
# from scipy.optimize import least_squares as least
from scipy.optimize import newton
from scipy.optimize import minimize
import numpy as np
import re
import ast

# LOG FILE
logfile = os.getcwd() + '/' + 'log.GSM+GSMCC_HOpt.' + time.strftime( "%y.%m.%d-%H.%M", time.localtime() )
# LOG NAME
logname = 'log.WD_GSM+GSMCC_HOpt.' + time.strftime( "%y.%m.%d-%H.%M", time.localtime() )
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
    # Print input array
    print_twice('\n OPTIMIZATION ITERATION')
    print_twice('\n X TIMES ALL THE OPTIMIZED PARAMETERS:')
    print_twice(x)
    # Open GSMCC input file
    with open(readfilename_CC,'r') as gsmin:
        inputfile_lines = gsmin.read().split('\n')
    # 
    start_value = 0
    # The order of X is [ONEPROTON_VAL, ONENEUTRON_VAL, TB_VAL]
    k = 0
    if onebody_opt == 1:
        if one_proton_opt == 1:
            # EDIT ONE PROTON
            theline = searchline(readfilename_CC,"core.potential")
            shift = [line_aux.strip(' ') for line_aux in inputfile_lines[theline:theline+10]].index('proton') + 2
            i = 0
            k = 0
            while i == 0:
                aux = inputfile_lines_start[theline + shift + k].split()
                if int(aux[0]) in one_proton_partialwaves:
                    if one_proton_independent == 'YES':
                        lwave_index = one_proton_partialwaves.index(int(aux[0]))
                        aux_vo = float(aux[3])*x[start_value+lwave_index]
                    else:
                        aux_vo = float(aux[3])*x[start_value]
                    inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+str(aux_vo)+'  '+aux[4]
                k += 1
                if inputfile_lines[theline + shift + k].split() == []:
                    i = 1
            start_value += 1
        if one_neutron_opt == 1:
            # EDIT ONE NEUTRON
            theline = searchline(readfilename_CC,"core.potential")
            shift = [x.strip(' ') for x in inputfile_lines[theline:theline+k+10]].index('neutron') + 2
            i = 0
            k = 0
            while i == 0:
                aux = inputfile_lines_start[theline + shift + k].split()
                if int(aux[0]) in one_neutron_partialwaves:
                    if one_neutron_independent == 'YES':
                        lwave_index = one_neutron_partialwaves.index(int(aux[0]))
                        if one_proton_opt == 1 and one_proton_independent == 'YES':
                            aux_vo = float(aux[3])*x[one_proton_npartialwaves+lwave_index]
                        elif one_proton_opt == 1 and one_proton_independent == 'NO':
                            aux_vo = float(aux[3])*x[start_value+lwave_index]
                        else:
                            aux_vo = float(aux[3])*x[lwave_index]
                    inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+str(aux_vo)+'  '+aux[4]
                k += 1
                if inputfile_lines[theline + shift + k].split() == []:
                    i = 1
            start_value += 1
    if twobody_opt == 1:
        # EDIT TWO BODY INTERACTIONS
        for i in range(n_tb_interactions):
            theline = searchline(readfilename_CC, tb_interactions[i])
            shift_from_nucleons = 0
            if one_proton_opt == 1 and one_proton_independent == 'YES': shift_from_nucleons += one_proton_npartialwaves
            if one_neutron_opt == 1 and one_neutron_independent == 'YES': shift_from_nucleons += one_neutron_npartialwaves
            if shift_from_nucleons == 0: shift_from_nucleons = start_value
            aux = float(inputfile_lines_start[theline].split()[0]) * x[shift_from_nucleons+i]
            inputfile_lines[theline] = '  ' + str(aux) + ' ' + inputfile_lines_start[theline].split()[1]
    #
    # Save and close GSMCC input file
    inputfile_aux = '\n'.join(inputfile_lines)
    with open(readfilename_CC,'w') as gsmin:
        gsmin.write(inputfile_aux)
    #
    aux = time.strftime( "%y.%m.%d-%H.%M.%S", time.localtime() )
    outfilename_CCi = logname + '/' + outfilename_CC + '-' + aux
    # outfilename_CCi = 're.52'
    start_gsmcc = time.time()
    print_twice('\n ' + running_prefix + running_cc + ' < ' + readfilename_CC + cc_write + outfilename_CCi)
    sp.run([running_prefix + running_cc + ' < ' + readfilename_CC + cc_write + outfilename_CCi], shell=True)
    end_gsmcc = time.time()
    time_gsmcc = end_gsmcc-start_gsmcc
    print_twice("Time to calculate: ",time_gsmcc, "s")
    #
    # Read results
    with open(outfilename_CCi,'r') as gsmout:
        outputfile_lines = gsmout.read().split('\n')
    line_numbers = searchline_all(outfilename_CCi,search)
    jpi_i = []
    ene_wid = []
    for line in line_numbers:
        auxjpi = outputfile_lines[line-10].split()[0]
        auxindex = re.findall(r'\((.*?)\)',outputfile_lines[line-10])[0]
        jpi_i.append([auxjpi,auxindex])
        auxene = float(outputfile_lines[line].split(':')[1].split(' ')[1])
        wid = float(outputfile_lines[line+1].split(':')[1].split(' ')[1])
        ene_wid.append([auxene,wid])
    #
    print_twice('Calculated energies and widths:')
    for i, state in enumerate(zip(jpi_i,ene_wid)):
        print_twice('  {0:s} ({1:s}) : E = {2:10.6f}, W = {3:10.6f};  '.format(state[0][0], state[0][1], state[1][0], state[1][1]))
    print_twice(' ')
    #
    # Compare with the experimental energy or separation energy of the selected JPi states
    if opt_sepenergy == 1:
        print("DEPRECATED!")
        exit()
        # res = 0
        # for i in range(n_jpi_states):
        #     if jpi_states[i] in jpi:
        #         related_index = jpi.index(jpi_states[i])
        #         calc_sep_energy = ene_wid[related_index][jpi_states_index[i][0]][0] - ene_wid[related_index][jpi_states_index[i][1]][0]
        #         print_twice('Calculated separation energy for JPi={0:s}: {1:10.6f} MeV'.format(jpi_states[i],calc_sep_energy))
        #         res_i = m.sqrt( abs(calc_sep_energy**2 - exp_value[i]**2) )
        #         print_twice('Residue as sqrt(calc**2 - exp**2): {0:10.6f}'.format(res_i))
        #         res += res_i
        #     else:
        #         print_twice('The state JPi=%s has not been calculated!!' % jpi_states[i])
        # print_twice('Sum of calculated residues: {0:10.6f}'.format(res))
    if opt_energy == 1:
        res = 0
        for i, jpi in enumerate(jpi_states):
            index = jpi_states_index[i]
            state = [jpi, str(index)]
            if state in jpi_i:
                calc_energy = ene_wid[jpi_i.index(state)][0]
                delta_e = jpi_states_energy[i] - calc_energy
                res_i = m.sqrt(abs(jpi_states_energy[i]**2 - calc_energy**2))
                res += res_i
                print_twice('For JPi: %s (%s) : Delta E = %s, Res = %s'% (jpi, index, delta_e, res_i))
    #
    print_twice('Sum of calculated residues: {0:10.6f}'.format(res))
    # Check which optimizator we are using
    if method == 'NEWTON':
        print_twice('Finish iteration of Newton optimizer, x and f(x) must be of the same size!')
        print_twice('\n'+20*'-'+'\n')
        return [res]
    elif method == 'MINIMIZATION':
        print_twice('Finish iteration of MINIMIZATION optimizer')
        print_twice('\n'+20*'-'+'\n')
        return res
    else:
        print_twice('METHOD must be NEWTON or MINIMIZATION;(TNC or Nelder-Mead or BFGS)')
        exit()
# .-

# INPUT FILE
readfilename = os.getcwd() + '/input.GSM+GSMCC_run'
with open(readfilename, 'r') as readfile:
    data = readfile.read().split('\n')
# LOGILE
erease_output_file()
    
print_twice("Be sure that you are using the correct directories!")
# GSMCC directory
gsmcc_directory = os.getcwd()
# GSM directory
theline = searchline(readfilename,"GSM-DIRECTORY:")
gsm_directory = data[theline+1]
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
    running_prefix = 'mpirun -np ' + parallelism_nodes + ' ./'
elif parallelism_type == 'OPENMP':
    running_prefix = ' ./'
elif parallelism_type == 'SLURM':
    running_prefix = 'srun -N ' + parallelism_nodes + ' '
else:
    print_twice('Parallelism must be MPI or OPENMP')
# Checking if we need machinefile
theline = searchline(readfilename,"MACHINEFILE:")
if theline != None:  
    machinefile_name = data[theline+1]
    running_prefix = running_prefix[:-2] + '-hostfile ' + machinefile_name + ' '
# Using experimental thresholds in GSMCC calculations
exp_thr = 0
theline = searchline(readfilename,"EXPERIMENTAL-THRESHOLDS:")
if theline != None:
    exp_thr = int(data[theline+1])
# Executable GSMCC file
theline = searchline(readfilename,"GSMCC-exe:")
running_cc = data[theline+1]
# Read-Out file name CC
theline = searchline(readfilename,"GSMCC-files:")
readfilename_CC = data[theline+1]
outfilename_CC = data[theline+2]
cc_write_aux = int(data[theline+3])
cc_write = ' ' + cc_write_aux*'>' + ' '
# Executable GSM file
theline = searchline(readfilename,"GSM-exe:")
if theline != None:  
    running_gsm = data[theline+1]
# Read-Out file name GSM
theline = searchline(readfilename,"GSM-files:")
if theline != None:  
    gsm_files = int(data[theline+1])
    gsm_write_aux = int(data[theline+2])
    gsm_write = ' ' + gsm_write_aux*'>' + ' '
    readfilename_GSM = data[theline+3:(theline+3)+2*gsm_files:2]
    outfilename_GSM = data[theline+4:(theline+4)+2*gsm_files:2]
    if len(readfilename_GSM) > gsm_files:
        readfilename_GSM = readfilename_GSM[:-1]
    
    
# CALCULATING FILE
calcfilename = os.getcwd() + '/input.GSM+GSMCC_HOpt'
with open(calcfilename, 'r') as readfile:
    data = readfile.read().split('\n')
#
# Looking up the optimization code to use
theline = searchline(calcfilename,"OPTIMIZATIONMETHOD:")
method =  data[theline+1].split(';')[0]
if method == 'MINIMIZATION':
    mini_method = data[theline+1].split(';')[1]
if method == 'MINIMIZATION':
    print_twice("Optimization method: %s; %s"% (method,mini_method))
else:
    print_twice("Optimization method: %s"% method)
#
# Will one-body be optimized?
onebody_opt = 0
theline = searchline(calcfilename,"OPT_ONEBODY:")
if theline != None:
    onebody_opt = 1
    ##
    one_proton_opt = 0
    one_neutron_opt = 0
    theline = searchline(calcfilename,"OPT_ONEPROTON:")
    if theline != None:
        one_proton_opt = 1
        # FOR THE MOMENT, ONLY WS WILL BE OPTIMIZED
        one_proton_type = data[theline+1]
        if one_proton_type == 'WS':
            print_twice('One proton WS will be optimized:')
            one_proton_partialwaves = data[theline+2]
            print_twice('  Partial waves to optimize: %s'% one_proton_partialwaves)
            one_proton_partialwaves = ast.literal_eval(one_proton_partialwaves)
            one_proton_npartialwaves = len(one_proton_partialwaves)
            one_proton_independent = data[theline+3].upper()
            one_proton_seed = float(data[theline+4])
            one_proton_bounds = data[theline+5]
    theline = searchline(calcfilename,"OPT_ONENEUTRON:")
    if theline != None:
        one_neutron_opt = 1
        # FOR THE MOMENT, ONLY WS WILL BE OPTIMIZED
        one_neutron_type = data[theline+1]
        if one_neutron_type == 'WS':
            print_twice('One neutron WS will be optimized:')
            one_neutron_partialwaves = data[theline+2]
            print_twice('  Partial waves to optimize: %s'% one_neutron_partialwaves)
            one_neutron_partialwaves = ast.literal_eval(one_neutron_partialwaves)
            one_neutron_npartialwaves = len(one_neutron_partialwaves)
            one_neutron_independent = data[theline+3].upper()
            one_neutron_seed = float(data[theline+4])
            one_neutron_bounds = data[theline+5]
# Will two-body be optimized?
twobody_opt = 0
theline = searchline(calcfilename,"OPT_TWOBODY:")
if theline != None:
    twobody_opt = 1
    print_twice('Two-body interactions will be optimized:')
    n_tb_interactions = int(data[theline+1])
    print_twice('  Number of two-body interactions to optimize: %d'% n_tb_interactions)
    tb_interactions = []
    for i in range(n_tb_interactions):
        print_twice('  Interaction %d: %s'% (i+1, data[theline+i+2]))
        tb_interactions.append( data[theline+i+2] )
    tb_seed = float(data[theline+2+n_tb_interactions])
    tb_bounds = data[theline+2+n_tb_interactions+1]
# Will SEPENERGY or ENERGY be optimized?
# theline = searchline(calcfilename,"OPTIMIZATION_OF:")
# opt_sepenergy = 0
# opt_energy = 0
# if data[theline+1] == 'SEPENERGY':
#     opt_sepenergy = 1
# elif data[theline+1] == 'ENERGY':
#     opt_energy = 1
# New version, only energy will be optimized
opt_energy = 1
opt_sepenergy = 0
#
# Reading JPi states to optimize their energies
theline = searchline(calcfilename,"JPI_STATES:")
jpi_states = data[theline+1].split(',')
jpi_states_index = ast.literal_eval( data[theline+2] )
print_twice("JPi states to optimize: ")
print_twice(list(zip(jpi_states,jpi_states_index)))
n_jpi_states = len(jpi_states)
#
# Reading experimental data
theline = searchline(calcfilename,"EXPERIMENTALVALUES:")
jpi_states_energy = [float(x) for x in data[theline+1].split(',')]
print_twice("Experimental separation energy to optimize: %s"% jpi_states_energy)
#
#
#
#
# Start calculation
theline = searchline(readfilename,"GSM-files:")
start_main = time.time()
# RUN GSM
if theline != None:
    print_twice("\nRunning GSM in %s"% gsm_directory)
    os.chdir(gsm_directory)
    #
    for i in range(0,gsm_files):
        inp = readfilename_GSM[i]
        out = outfilename_GSM[i]
        start_gsm = time.time()
        print_twice('\n ' + running_prefix + running_gsm  + ' < ' + inp + gsm_write + out)
        sp.run([running_prefix + running_gsm + ' < ' + inp + gsm_write + out], shell=True)
        end_gsm = time.time()
        time_gsm = end_gsm-start_gsm
        print_twice("Time to calculate: ",time_gsm, "s")
    #
else:
    print_twice("\nSkip GSM part, only GSMCC calculation!")
#
# Edit thresholds
if exp_thr == 0:
    print_twice("\nCalculations of GSMCC using GSM thresholds")
elif exp_thr == 1:
    print_twice("\nCalculations of GSMCC using experimental thresholds")
    print_twice("\nEdit thresholds in %s"% storage_directory)
    print_twice("\nBe sure that the file EditThreshold.in is in the directory!")
    os.chdir(storage_directory)
    sp.run(['python3 Useful_Codes/GSMCC-CODE/EditThresholds.py'], shell=True)
#
print_twice("\nRunning GSMCC in %s"% gsmcc_directory)
os.chdir(gsmcc_directory)
#
# Defining optimization part
search = 'E(reference frame) :'
#
# The seed values are defined in the GSMCC input file
# And since we are going to optimize an X value of proportion
# to the GSMCC values, the seed for X is one
# The order of the array is 
# [ONEPROTONSEED x NUMBER_OF_PARTIAL_WAVES, ONENEUTRONSEED x NUMBER_OF_PARTIAL_WAVES, TO_SEED, TE_SEED, SO_SEED, SE_SEED, SOTO_SEED, SOTE_SEED, TTO_SEED, TTE_SEED]
# For ONEBODY_SEED
onebody_seed = []
onebody_bounds = ''
if onebody_opt == 1 :
    if one_proton_opt == 1:
        onebody_bounds += one_proton_bounds
        if one_proton_independent == 'NO':
            onebody_seed.append(one_proton_seed)
        else:
            for i in range(one_proton_npartialwaves):
                onebody_seed.append(one_proton_seed)
    if one_neutron_opt == 1:
        if onebody_bounds != '':
            onebody_bounds += ','
        onebody_bounds += one_neutron_bounds
        if one_neutron_independent == 'NO':
            onebody_seed.append(one_neutron_seed)
        else:
            for i in range(one_neutron_npartialwaves):
                onebody_seed.append(one_neutron_seed)
# For TWOBODY_SEED
if twobody_opt == 1:
    twobody_seed = n_tb_interactions * [1]
else:
    twobody_seed = []
#
seeds = onebody_seed + twobody_seed
if onebody_bounds != '':
    bounds = (eval(onebody_bounds + ',' + tb_bounds))
else:
    bounds = (eval(tb_bounds + ',' ))
#
# Read GSMCC input file
with open(readfilename_CC,'r') as gsmin:
    inputfile_lines_start = gsmin.read().split('\n')
# Then calculation
if method == 'NEWTON':
    print_twice('Using Newton optimizer, x and f(x) must be of the same size!')
    opt = newton(f, seeds, tol=1e-10, maxiter=20, full_output=True)
elif method == 'MINIMIZATION':
    print_twice('Using %s optimizer'% mini_method)
    if mini_method == 'TNC':
        opt = minimize(f, seeds, method=mini_method, jac='2-point', options={ 'xtol' : 1e-3, 'finite_diff_rel_step': 0.001 }) # idea from https://stackoverflow.com/questions/20478949/how-to-force-larger-steps-on-scipy-optimize-functions
    elif mini_method == 'Nelder-Mead' or mini_method == 'BFGS':
        opt = minimize(f, seeds, method=mini_method, bounds=bounds)
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
    Example of input.GSM+GSMCC_HOpt file:
    _________________________________
    OPTIMIZATIONMETHOD:  'NEWTON' or 'MINIMIZATION;' + ('TNC' or 'Nelder-Mead' or 'BFGS')
    MINIMIZATION;BFGS
    
    OPT_ONEBODY: Define which one-body will be optimized
    OPT_ONEPROTON: 1 - Part to optimize WS, SO, R0, AA; 2 - l-waves to optimize [0,1...]; 3 - INDEPENDENT l-waves; 4 - GLOBAL CF SEED; 5 - Bounds
    WS
    [0,1]
    YES
    1
    (0.8,1.2)
    OPT_ONENEUTRON: 1 - Part to optimize WS, SO, R0, AA; 2 - l-waves to optimize [0,1...]; 3 - INDEPENDENT l-waves; 4 - GLOBAL CF SEED; 5 - Bounds
    WS
    [0,1,2,3]
    NO
    1
    (0.8,1.2)

    OPT_TWOBODY: 1 - number of TB interaction to opt; 2,.. - name of each as in GSMCC input; N+1 - GLOBAL CF SEED; N+2 - Bounds
    4
    (V0.NN.central.odd.triplet(S=1,T=1))
    (V0.NN.central.even.triplet(S=1,T=0))
    (V0.NN.central.odd.singlet(S=0,T=0))
    (V0.NN.central.even.singlet(S=0,T=1))
    1
    (0.8,1.2)

    JPI_STATES: 1 - JPi of each state to be optimized; 2 - index of the desired state
    3/2+,5/2+,5/2+,7/2+
    [0,0,1,0]

    EXPERIMENTALVALUES: Values of the energy in MeV of the desired JPi's
    6,2.5,2.9
    _________________________________
"""