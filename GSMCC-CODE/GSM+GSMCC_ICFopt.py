"""
    Created on October 2024 by Alan D.K. for 11C_Project

    This code run GSM (if needed) and GSMCC and then
    optimize the complex (or not) interaction corrective factor
    and/or cluster corrective factor in order
    to reproduce the experimental data. 
    You have to define the input and output files and put
    the experimental values (energy and width).
    
    An example of the input file "input.GMS+GSMCC_ICF" can be found
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
import re

# LOG NAME
logname = 'log.WD_GSM+GSMCC_ICF.' + time.strftime( "%y.%m.%d-%H.%M", time.localtime() )
# LOG FOLDER
logfolder = os.getcwd() + '/' + logname 
if not os.path.exists(logfolder):
    os.makedirs(logfolder)
else:
    for filename in os.listdir(logfolder):
        # remove the files inside
        os.remove(f"{logfolder}/{filename}")
# LOG FILE
logfile = logfolder + '/' + 'log.GSM+GSMCC_ICF.' + time.strftime( "%y.%m.%d-%H.%M", time.localtime() )

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
    print_twice('\n All the corrective factors:')
    print_twice(x)
    # Open GSMCC input file
    with open(readfilename_CC,'r') as gsmin:
        inputfile_lines = gsmin.read().split('\n')
    # 
    start_value = 0
    if used_icf == 1: # Using interaction corrective factors
        aux1 = inputfile_lines[icf_line].split()
        if icf_type == 'COMPLEX': # Complex interaction corrective factors
            start_value = 2
            aux2 = str(x[0])
            aux3 = str(x[1])
        elif icf_type == 'REAL': # Only real part optimization
            start_value = 1
            aux2 = str(x[0])
            aux3 = str(icf_imag_seed)
        elif icf_type == 'IMAG': # Only imag part optimization
            start_value = 1
            aux2 = str(icf_real_seed)
            aux3 = str(x[0])
        else:
            print_twice('PROBLEM WITH CCF_TYPE, MUST BE REAL, IMAG OR COMPLEX')
            exit()
        # 
        inputfile_lines[icf_line] = "  " + aux1[0] + " " + aux2 + " " + aux3
    if used_ccf == 1: # Using Cluster corrective factors
        if ccf_type == 'COMPLEX':
            # Now edit each ccf for all the clusters
            for i in range(0,ccf_numberofcluster):
                aux2 = inputfile_lines[ccf_lines[i]].split()
                inputfile_lines[ccf_lines[i]] = "      " + aux2[0] + " " + str(x[start_value+i*2]) + " " + str(x[start_value+1+i*2])
        elif ccf_type == 'REAL':
            # Now edit each ccf for all the clusters
            for i in range(0,ccf_numberofcluster):
                aux2 = inputfile_lines[ccf_lines[i]].split()
                inputfile_lines[ccf_lines[i]] = "      " + aux2[0] + " " + str(x[start_value+i]) + " " + str(ccf_imag_seed[i])
        elif ccf_type == 'IMAG':
            # Now edit each ccf for all the clusters
            for i in range(0,ccf_numberofcluster):
                aux2 = inputfile_lines[ccf_lines[i]].split()
                inputfile_lines[ccf_lines[i]] = "      " + aux2[0] + " " + str(ccf_real_seed[i]) + " " + str(x[start_value+i])
        else:
            print_twice('PROBLEM WITH CCF_TYPE, MUST BE REAL, IMAG OR COMPLEX')
            exit()
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
    # New version to avoid problems with doublets or loss of states
    res_e = np.zeros(numberofstates)
    res_w = np.zeros(numberofstates)
    indexmin_res_e = np.array(3*[-100])
    for i in range(0,numberofstates):
        expene = expene_read[i]
        expwid = expwid_read[i]
        stweig = stweig_read[i]
        numberindex = index_read[i]
        if index_search == 'YES':
            # Test with the index and index\pm1 states 
            if numberindex >= 1:
                index0_res_e = expene - float(auxiliar[numberindex-1][0])
                index1_res_e = expene - float(auxiliar[numberindex][0])
                index2_res_e = expene - float(auxiliar[numberindex+1][0])
                #
                aux = numberindex + np.argmin( [ abs(index0_res_e), abs(index1_res_e), abs(index2_res_e) ] ) - 1
                if np.size(np.argwhere(indexmin_res_e == aux)) == 0:
                    indexmin_res_e[i] = aux
                else:
                    # Second minimum to avoid two states assigned to a same eigenvalue
                    aux2 = np.argsort( [ abs(index0_res_e), abs(index1_res_e), abs(index2_res_e) ] )[1] - 1
                    indexmin_res_e[i] = numberindex + aux2
            else:
                if len(line_numbers) > 1:
                    index1_res_e = expene - float(auxiliar[numberindex][0])
                    index2_res_e = expene - float(auxiliar[numberindex+1][0])
                    #
                    aux = [ abs(index1_res_e), abs(index2_res_e) ]
                    indexmin_res_e[i] = numberindex + np.argmin( aux )
                else:
                    indexmin_res_e[i] = numberindex
            
            res_e_aux = expene - float(auxiliar[indexmin_res_e[i]][0])
            res_e[i] = ( res_e_aux )**2 / abs(expene) * stweig
            res_w_aux = expwid - float(auxiliar[indexmin_res_e[i]][1])
            if expwid < 1e-10: # Zero widht!
                res_w[i] = ( res_w_aux )**2 * stweig
            else:
                res_w[i] = ( res_w_aux )**2 / abs(expwid) * stweig
            print_twice('State selected Index : {0:d}, Real index : {1:d}\n  E Residue = {2:10.6f} MeV, W Residue = {3:10.6f} keV'.format(numberindex,indexmin_res_e[i],res_e_aux,res_w_aux))
        elif index_search == 'NO':
            res_e_aux = expene - float(auxiliar[numberindex][0])
            res_e[i] = ( res_e_aux )**2 / abs(expene) * stweig
            res_w_aux = expwid - float(auxiliar[numberindex][1])
            if expwid < 1e-10: # Zero widht!
                res_w[i] = ( res_w_aux )**2 * stweig
            else:
                res_w[i] = ( res_w_aux )**2 / abs(expwid) * stweig
            print_twice('State selected Index : {0:d}\n  E Residue = {1:10.6f} MeV, W Residue = {2:10.6f} keV'.format(numberindex,res_e_aux,res_w_aux))
        else:
            print_twice('YES or NO for index searching!')
            exit()
    if optimize_width == 0:
        res = np.sum(res_e)
    elif optimize_energy == 0:
        res = np.sum(res_w)
    else:
        res = np.sum(res_e) + np.sum(res_w)/1000
    print_twice('X^2 = {0:10.6f}'.format(res))
    #
    # Check which optimizator we are using
    if method == 'NEWTON':
        print_twice('Finish iteration of Newton optimizer, x and f(x) must be of the same size!')
        if optimize_width == 0:
            return_residue = res_e
        elif optimize_energy == 0:
            return_residue = res_w
        else:
            aux = list(res_e) + list(res_w)
            aux[::2] = list(res_e)
            aux[1::2] = list(res_w)
            return_residue = aux
        print_twice('\n'+20*'-'+'\n')
        return return_residue
    elif method == 'MINIMIZATION':
        print_twice('Finish iteration of MINIMIZATION optimizer')
        print_twice('\n'+20*'-'+'\n')
        return res
    else:
        print_twice('METHOD must be NEWTON or MINIMIZATION;(TNC or Nelder-Mead or BFGS)')
        exit()
# .-

# INPUT FILE
readfilename = "input.GSM+GSMCC_ICF"
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
#
# Looking up the optimization code to use
theline = searchline(readfilename,"OPTIMIZATIONMETHOD:")
method =  data[theline+1].split(';')[0]
if method == 'MINIMIZATION':
    mini_method = data[theline+1].split(';')[1]
#
# The code will optimize energy, width or both
theline = searchline(readfilename,"OPTIMIZATION_OF:")
optimize_energy = 0
optimize_width = 0
optimize_read = data[theline + 1 ]
if optimize_read == 'ENERGY':
    optimize_energy = 1
    print_twice('Only the energies will be adjusted')
elif optimize_read == 'WIDTH':
    optimize_width = 1
    print_twice('Only the widths will be adjusted')
elif optimize_read == 'BOTH':
    optimize_energy = 1
    optimize_width = 1
    print_twice('The energies and widths will be adjusted')
else:
    print_twice('ENERGY, WIDHT or BOTH must be the optimization code!')
    exit()
#
# Checking if real or complex interaction corrective factors will be used
used_icf = 0
icf_type = 'NONE'
icf_bounds = ''
theline = searchline(readfilename,"CORRECTIVEFACTORS:")
if theline != None:
    used_icf = 1
    icf_type = data[theline+1]
    icf_real_seed = float(data[theline+2])
    icf_imag_seed = float(data[theline+3])
    icf_bounds = data[theline+4]
#
# Checking if real or complex clusters corrective factors will be used
used_ccf = 0
ccf_type = 'NONE'
ccf_bounds = ''
theline = searchline(readfilename,"CLUSTERCORRFACTORS:")
if theline != None:
    used_ccf = 1
    ccf_type = data[theline+1]
    ccf_numberofcluster = int(data[theline+2])
    ccf_clusters = ccf_numberofcluster*['0']
    ccf_real_seed = ccf_numberofcluster*[0]
    ccf_imag_seed = ccf_numberofcluster*[0]
    #
    for i in range(0,ccf_numberofcluster):
        factor = i*4
        ccf_clusters[i] = data[theline+3+factor]
        ccf_real_seed[i] = float(data[theline+4+factor])
        ccf_imag_seed[i] = float(data[theline+5+factor])
        if i == 0:
            ccf_bounds = ccf_bounds + data[theline+6+factor]
        else:
            ccf_bounds = ccf_bounds + ',' + data[theline+6+factor]
#
# Reading experimental data
theline = searchline(readfilename,"EXPERIMENTALVALUES:")
numberofstates = int(data[theline+1])
index_search = data[theline+2]
expene_read = np.zeros(numberofstates)
expwid_read = np.zeros(numberofstates)
stweig_read = np.zeros(numberofstates)
index_read = np.zeros(numberofstates, dtype=int)
for i in range(0,numberofstates):
    expene_read[i] = float(data[theline+3+i*4])
    expwid_read[i] = float(data[theline+4+i*4])
    index_read[i] = int(data[theline+5+i*4])
    stweig_read[i] = float(data[theline+6+i*4])
#

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
#
# Start calculation
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
# Line with the interaction corrective factor
if used_icf == 1:
    icf_line = searchline(readfilename_CC,"CC.interaction.corrective.factor.composite(s)")+2
# Lines with the cluster corrective factor
if used_ccf == 1:
    ccf_lines = np.zeros(ccf_numberofcluster,dtype=int)
    for i in range(0,ccf_numberofcluster):
        ccf_lines[i] = searchline(readfilename_CC,"CC.corrective.factor.%s.composite(s)"% ccf_clusters[i])+2
#
# if icf_type == 'COMPLEX' and ccf_type == 'COMPLEX':
# Define seed value first
if used_ccf == 0 and used_icf == 1:
    if icf_type == 'COMPLEX':
        seeds = [icf_real_seed,icf_imag_seed]
        print_twice('Optimizing real and imaginary parts of ICF')
    elif icf_type == 'REAL':
        seeds = [icf_real_seed]
        print_twice('Optimizing real part of ICF')
    elif icf_type == 'IMAG':
        seeds = [icf_imag_seed]
        print_twice('Optimizing imaginary part of ICF')
    else:
        print_twice('PROBLEM WITH ICF_TYPE, MUST BE REAL, IMAG OR COMPLEX')
        exit()
elif used_ccf == 1:
    if used_icf == 1:
        if icf_type == 'COMPLEX':
            # The first two seeds are always the real and/or imaginary part of the interaction corrective factors
            seeds_aux1 = [icf_real_seed,icf_imag_seed]
            print_twice('Optimizing real and imaginary parts of ICF')
        elif icf_type == 'REAL':
            seeds_aux1 = [icf_real_seed]
            print_twice('Optimizing real part of ICF')
        elif icf_type == 'IMAG':
            seeds_aux1 = [icf_imag_seed]
            print_twice('Optimizing imaginary part of ICF')
        else:
            print_twice('PROBLEM WITH ICF_TYPE, MUST BE REAL, IMAG OR COMPLEX')
            exit()
    else:
        seeds_aux1 = []
    if ccf_type == 'COMPLEX':
        # Then, we define pairs of ccf for all the clusters
        # Idea from https://www.geeksforgeeks.org/python-interleave-multiple-lists-of-same-length/
        seeds_aux2 = ccf_real_seed + ccf_imag_seed 
        seeds_aux2[::2] = ccf_real_seed
        seeds_aux2[1::2] = ccf_imag_seed
        print_twice('Optimizing real and imaginary parts of CCF')
    elif ccf_type == 'REAL':
        seeds_aux2 = ccf_real_seed
        print_twice('Optimizing real part of CCF')
    elif ccf_type == 'IMAG':
        seeds_aux2 = ccf_imag_seed
        print_twice('Optimizing imaginary part of CCF')
    else:
        print_twice('PROBLEM WITH CCF_TYPE, MUST BE REAL, IMAG OR COMPLEX')
        exit()
    #
    seeds = seeds_aux1 + seeds_aux2
# Then calculation
if method == 'NEWTON':
    print_twice('Using Newton optimizer, x and f(x) must be of the same size!')
    opt = newton(f, seeds, tol=1e-10, maxiter=20, full_output=True)
elif method == 'MINIMIZATION':
    print_twice('Using %s optimizer'% mini_method)
    if mini_method == 'TNC':
        opt = minimize(f, seeds, method=mini_method, jac='2-point', options={ 'xtol' : 1e-3, 'finite_diff_rel_step': 0.001 }) # idea from https://stackoverflow.com/questions/20478949/how-to-force-larger-steps-on-scipy-optimize-functions
    elif mini_method == 'Nelder-Mead' or mini_method == 'BFGS':
        if used_icf == 1:
            bounds_opt = (eval(icf_bounds + ',' + ccf_bounds))
        else:
            bounds_opt = ( eval(ccf_bounds + ',') )
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
    Example of input.GSM+GSMCC_ICF file:
    _________________________________
    GSM-DIRECTORY:
    /home/dassie/2024/Carbon-11_Porject/GSM-24.02/GSM_dir_2D/GSM_dir
    STORAGE-DIRECTORY:
    /home/dassie/2024/Carbon-11_Porject/GSM-24.02/GSM_dir_2D/storage_11C_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30
    
    PARALLELISM: 1 - MPI or OPENMP; 2 - NUMBER OF NODES
    MPI
    2
    MACHINEFILE: DEFINE THE FILE FOR THE EXECUTION, IF NEEDED
    machinefile
    
    OPTIMIZATIONMETHOD:  'NEWTON' or 'MINIMIZATION;' + ('TNC' or 'Nelder-Mead' or 'BFGS')
    MINIMIZATION;TNC
    
    OPTIMIZATION_OF: 'ENERGY', 'WIDTH' or 'BOTH'
    ENERGY
    
    CORRECTIVEFACTORS: 1 - COMPLEX, REAL or IMAG; 2 - REAL SEED; 3 - IMAG SEED; 4 - IF Nelder-Mead, DEFINE COMPLEX OR REAL BOUNDS IN A FORM (R.MIN,R.MAX),(I.MIN,I.MAX)
    REAL
    1.025
    0.0
    (0.5,1.5)

    CLUSTERCORRFACTORS: 1 - COMPLEX, REAL or IMAG; 2 - NUMBER OF CLUSTERS; for each cluster -> 3 - NAME OF THE CLUSTER; 4 - REAL SEED; 5 - IMAG SEED; 6 - IF Nelder-Mead, DEFINE COMPLEX OR REAL BOUNDS IN A FORM (R.MIN,R.MAX),(I.MIN,I.MAX)
    COMPLEX
    1
    alpha
    1.0
    0.0
    (0,8),(-1,1)

    EXPERIMENTALVALUES: 1 - NUMBER OF EXPERIMENTAL STATES; 2 - "YES" FOR SEARCH INDEX OR "NO" FOR NOT; for each state -> 3 - ENERGY (MeV); 4 - WIDTH (keV); 5 - ESTIMATED ENERGY ORDERED INDEX; 6 - WEIGHT BETWEEN (0,1)
    2
    NO
    -36.446
    15.0
    1
    1
    -36.908
    12.
    2
    0.5
    
    
    GSM-exe:
    GSM-24.11.20-MPI.x
    GSM-files: 1-n input files; 2-1 for overwrite or 2 for append; 3,n-input and output names
    4
    1
    CLUSPHY_targforCC_10B_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.in
    CLUSPHY_targforCC_10B_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.out
    CLUSPHY_targforCC_7Be_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.in
    CLUSPHY_targforCC_7Be_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.out
    CLUSPHY_projforCC_protonnocore_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.in
    CLUSPHY_projforCC_protonnocore_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.out
    CLUSPHY_projforCC_alphanocore_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.in
    CLUSPHY_projforCC_alphanocore_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.out
    
    EXPERIMENTAL-THRESHOLDS: Using experimental thresholds, 0 means NO and 1 means YES
    0
    
    GSMCC-exe:
    CC-24.11.20-MPI.x
    GSMCC-files: 1-input file; 2-output file; 3-1 for overwrite or 2 for append
    CLUSPHY_11C_CC_GSMOpt-24.08.26-11.00_Basis-24.10.24-17.00_ICF-24.10.24-17.30.in
    CLUSPHY_11C_CC_GSMOpt-24.08.26-11.00_Basis-24.10.24-17.00_ICF-24.10.24-17.30_7I2+.out
    2
    _________________________________
"""