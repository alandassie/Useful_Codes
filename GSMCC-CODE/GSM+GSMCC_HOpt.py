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
    if separate_optimization == 0:
        # The order of X is [ONEPROTON_VAL, ONENEUTRON_VAL, TB_VAL]
        k = 0
        if one_proton_opt == 1:
            # EDIT ONE PROTON
            theline = searchline(readfilename_CC,"core.potential")
            shift = [x.strip(' ') for x in inputfile_lines[theline:theline+10]].index('proton') + 2
            i = 0
            k = 0
            while i == 0:
                aux = inputfile_lines_start[theline + shift + k].split()
                if one_proton_partialwaves == 'ALL':
                    aux_vo = float(aux[3])*x[start_value]
                    inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+str(aux_vo)+'  '+aux[4]
                elif int(aux[0]) in one_proton_partialwaves:
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
                if one_neutron_partialwaves == 'ALL':
                    aux_vo = float(aux[3])*x[start_value]
                    inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+str(aux_vo)+'  '+aux[4]
                elif int(aux[0]) in one_neutron_partialwaves:
                    aux_vo = float(aux[3])*x[start_value]
                    inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+str(aux_vo)+'  '+aux[4]
                k += 1
                if inputfile_lines[theline + shift + k].split() == []:
                    i = 1
            start_value += 1
    #
    # Save and close GSMCC input file
    inputfile_aux = '\n'.join(inputfile_lines)
    with open(readfilename_CC,'w') as gsmin:
        gsmin.write(inputfile_aux)
    #
    aux = time.strftime( "%y.%m.%d-%H.%M.%S", time.localtime() )
    outfilename_CCi = logname + '/' + outfilename_CC + '-' + aux
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
    energy = []
    width = []
    j = 0
    for line in line_numbers:
        auxe = outputfile_lines[line]
        ene = float(auxe.split(':')[1].split(' ')[1])
        auxe_pole = outputfile_lines[line-4]
        ene_pole = float(auxe_pole.split(':')[1].split(' ')[1])
        auxw = outputfile_lines[line+1]
        wid = float(auxw.split(':')[1].split(' ')[1])
        # Only append when ene and ene_pole are close
        if abs(ene - ene_pole) < 8: # Limit on 8 MeV
            energy.append( ene )
            width.append( wid )
        else:
            print_twice('OBS! The state with index %s differs from the pole approx more than 8 MeV!'% j)
            print_twice('E_pole = {0:10.6f}, E = {1:10.6f}'.format(ene_pole,ene))
            energy.append( ene )
            width.append( wid )
            # print_twice('It will be discarded!')
        j += 1
    #
    auxiliar = list(zip(energy,width))
    auxiliar.sort()
    print_twice('Calculated energies and widths:')
    print_twice(auxiliar)
    # Compare with the experimental energy or separation energy
    if opt_sepenergy == 1:
        calc_sep_energy = auxiliar[0][0] - auxiliar[1][0]
        print_twice('Calculated separation energy: {0:10.6f} MeV'.format(calc_sep_energy))
        res = m.sqrt(abs(calc_sep_energy**2 - exp_value**2))
    if opt_energy == 1:
        res = m.sqrt(abs(auxiliar[0][0]**2 - exp_value**2))
    #
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
#
# Will one-body be optimized?
onebody_opt = 0
theline = searchline(calcfilename,"OPT_ONEBODY:")
if theline != None:
    onebody_opt = 1
    separate_optimization = int(data[theline+1])
    ##
    if separate_optimization == 1:
        print_twice('l-WAVE SEPARATE OPTIMIZATION IS NOT CODED YET!')
        exit()
    ##
    one_proton_opt = 0
    one_neutron_opt = 0
    theline = searchline(calcfilename,"OPT_ONEPROTON:")
    if theline != None:
        one_proton_opt = 1
        # FOR THE MOMENT, ONLY WS WILL BE OPTIMIZED
        one_proton_type = data[theline+1]
        if one_proton_type == 'WS':
            one_proton_partialwaves = data[theline+2]
            if one_proton_partialwaves == 'ALL' and separate_optimization == 1:
                one_proton_npartialwaves = data[theline+3]
            elif one_proton_partialwaves != 'ALL':
                one_proton_partialwaves = ast.literal_eval(one_proton_partialwaves)
                one_proton_npartialwaves = len(one_proton_partialwaves)
    theline = searchline(calcfilename,"OPT_ONENEUTRON:")
    if theline != None:
        one_neutron_opt = 1
        # FOR THE MOMENT, ONLY WS WILL BE OPTIMIZED
        one_neutron_type = data[theline+1]
        if one_neutron_type == 'WS':
            one_neutron_partialwaves = data[theline+2]
            if one_neutron_partialwaves == 'ALL' and separate_optimization == 1:
                one_neutron_npartialwaves = data[theline+3]
            elif one_neutron_partialwaves != 'ALL':
                one_neutron_partialwaves = ast.literal_eval(one_neutron_partialwaves)
                one_neutron_npartialwaves = len(one_neutron_partialwaves)
# Will two-body be optimized?
twobody_opt = 0
theline = searchline(calcfilename,"OPT_TWOBODY:")
if theline != None:
    print("NOT PROGRAMED YET!")
    exit()
# Will SEPENERGY or ENERGY be optimized?
theline = searchline(calcfilename,"OPTIMIZATION_OF:")
opt_sepenergy = 0
opt_energy = 0
if data[theline+1] == 'SEPENERGY':
    opt_sepenergy = 1
elif data[theline+1] == 'ENERGY':
    opt_energy = 1
#
# Reading experimental data
theline = searchline(calcfilename,"EXPERIMENTALVALUE:")
exp_value= float(data[theline+1])
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
# Open GSMCC input file
# with open(readfilename_CC,'r') as gsmin:
#     inputfile_lines = gsmin.read().split('\n')
# if one_proton_opt == 1:
#     # Line with the first proton partial wave
#     aux = searchline(readfilename_CC,"core.potential")
#     shift = [x.strip(' ') for x in inputfile_lines[theline:theline+10]].index('proton') + 2
#     one_proton_line = aux + shift
# # Line with the interaction corrective factor
# if used_icf == 1:
#     icf_line = searchline(readfilename_CC,"CC.interaction.corrective.factor.composite(s)")+2
# # Lines with the cluster corrective factor
# if used_ccf == 1:
#     ccf_lines = np.zeros(ccf_numberofcluster,dtype=int)
#     for i in range(0,ccf_numberofcluster):
#         ccf_lines[i] = searchline(readfilename_CC,"CC.corrective.factor.%s.composite(s)"% ccf_clusters[i])+2
#
# The seed values are defined in the GSMCC input file
# And since we are going to optimize an X value of proportion
# to the GSMCC values, the seed for X is one
# If separate optimization, the order of the array is 
# [ONEPROTONSEED x NUMBER_OF_PARTIAL_WAVES, ONENEUTRONSEED x NUMBER_OF_PARTIAL_WAVES, TO_SEED, TE_SEED, SO_SEED, SE_SEED, SOTO_SEED, SOTE_SEED, TTO_SEED, TTE_SEED]
# If not, the order of the array is
# [ONEPROTONSEED, ONENEUTRONSEED, TB_SEED]
# For ONEBODY_SEED
if separate_optimization == 1:
    if one_proton_opt == 1 and one_neutron_opt == 1:
        onebody_seed = [1 for i in range(0,one_proton_npartialwaves)] + [1 for i in range(0,one_neutron_npartialwaves)] 
    elif one_proton_opt == 1:
        onebody_seed = [1 for i in range(0,one_proton_npartialwaves)]
    elif one_neutron_opt == 1:
        onebody_seed = [1 for i in range(0,one_neutron_npartialwaves)] 
else:
    if one_proton_opt == 1 and one_neutron_opt == 1:
        onebody_seed = [1, 1]
    else:
        onebody_seed = [1]
# FOR THE MOMENT ONLY ONE-BODY OPT
seeds = onebody_seed # + twobody_seed
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
        opt = minimize(f, seeds, method=mini_method)
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
    
    OPT_ONEBODY: 0 means all the l-waves have the same proportion, 1 on contrary
    0
    OPT_ONEPROTON: 1 - Part to optimize WS, SO, R0, AA; 2 - l-waves to optimize ALL or [0,1...]
    WS
    ALL
    OPT_ONENEUTRON: 1 - Part to optimize WS, SO, R0, AA; 2 - l-waves to optimize ALL or [0,1...]
    WS
    ALL
        
    OPTIMIZATION_OF: ENERGY for only energy of one state or SEPENERGY for the distance between two calculated states
    SEPENERGY
    
    EXPERIMENTALVALUE: Value of the energy or separation energy (no sign) in MeV
    0.5
    _________________________________
"""