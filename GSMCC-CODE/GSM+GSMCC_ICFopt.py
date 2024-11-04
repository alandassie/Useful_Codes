"""
    Created on October 2024 by Alan D.K. for 11C_Project

    This code run GSM (if needed) and GSMCC and then
    optimize the complex (or not) interaction corrective factor
    to reproduce the experimental data. 
    You have to define the input and output files and put
    the experimental values (energy and width).
    
    An example of the input file "GMS+GSMCC_ICF.in" can be found
    at the end of the code.
"""

import subprocess as sp
import time
import os
import math as m
# from scipy.optimize import least_squares as least
from scipy.optimize import newton

# LOG FILE
logfile = os.getcwd() + '/GSM+GSMCC_ICF.log'

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
    print_twice('\n Correction factors:')
    print_twice(x)
    # Open GSMCC input file
    with open(readfilename_CC,'r') as gsmin:
        inputfile_lines = gsmin.read().split('\n')
    aux = inputfile_lines[icf_line].split()
    if len(x) == 2: # Complex corrective factor
        inputfile_lines[icf_line] = "  " + aux[0] + " " + str(x[0]) + " " + str(x[1])
    else: # Real corrective factor
        inputfile_lines[icf_line] = "  " + aux[0] + " " + str(x[0]) + " 0.0 "
    # Save and close GSMCC input file
    inputfile_aux = '\n'.join(inputfile_lines)
    with open(readfilename_CC,'w') as gsmin:
        gsmin.write(inputfile_aux)
    #
    outfilename_CCi = outfilename_CC + str(int(time.time()))
    start_gsmcc = time.time()
    print_twice('\n ' + running_prefix + './CC_exe < '+readfilename_CC+' > '+outfilename_CCi)
    sp.run([running_prefix + './CC_exe < '+readfilename_CC+' > '+outfilename_CCi], shell=True)
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
    for line in line_numbers:
        auxe = outputfile_lines[line]
        ene = float(auxe.split(':')[1].split(' ')[1])
        auxw = outputfile_lines[line+1]
        wid = float(auxw.split(':')[1].split(' ')[1])
        energy.append( ene )
        width.append( wid )
    #
    auxiliar  = list(zip(energy,width))
    auxiliar.sort()
    print_twice('Calculated energies and widths:')
    print_twice(auxiliar)
    res_e = expene - float(auxiliar[numberindex][0])
    res_w = expwid - float(auxiliar[numberindex][1])
    print_twice('E Residue = {0:7.3f}, W Residue = {1:10.6f}'.format(res_e,res_w))
    res = m.sqrt(res_e**2 + res_w**2)
    print_twice('Sum^2 Residue = {0:7.3f}'.format(res))
    print_twice('\n'+20*'-'+'\n')
    #
    if len(x) == 2: # Complex corrective factor
        return_residue = [res_e,res_w]
    else:
        return_residue = res_e
    #
    return return_residue
# .-

# INPUT FILE
readfilename = "GSM+GSMCC_ICF.in"
with open(readfilename, 'r') as readfile:
    data = readfile.read().split('\n')
# LOGILE
erease_output_file()
    
print_twice("Be sure that you are using the correct directories!")
# GSMCC directory
gsmcc_directory = os.getcwd()
# GSM directory
theline = searchline(readfilename,"GSM-DIRECTORY")
gsm_directory = data[theline+1]
# sure = input("Is %s the GSM working directory?\n YES or NO: "% gsm_directory)
# if sure.lower() == "no":
#     print("Change it in python file!")
#     exit()
# storage directory
theline = searchline(readfilename,"STORAGE-DIRECTORY")
storage_directory = data[theline+1]
# sure = input("Is %s the storage directory?\n YES or NO: "% storage_directory)
# if sure.lower() == "no":
#     print("Change it in python file!")
#     exit()
# Checking if it is a MPI or OPENMP/secuential calculation
theline = searchline(readfilename,"PARALLELISM")
parallelism_type = data[theline+1]
parallelism_nodes = data[theline+2]
if parallelism_type == 'MPI':
    running_prefix = 'mpirun -np ' + parallelism_nodes + ' '
else:
    running_prefix = ' '
# Checking if we need machinefile
theline = searchline(readfilename,"MACHINEFILE")
if theline != None:  
    machinefile_name = data[theline+1]
    running_prefix = running_prefix + '-hostfile ' + machinefile_name + ' '
# Checking if real or complex corrective factors will be used
theline = searchline(readfilename,"CORRECTIVEFACTORS")
icf_type = data[theline+1]
icf_real_seed = float(data[theline+2])
if icf_type == 'COMPLEX':
    icf_imag_seed = float(data[theline+3])
# Reading experimental data
theline = searchline(readfilename,"EXPERIMENTALVALUES")
expene_read = float(data[theline+1])
expwid_read = float(data[theline+2])
index_read = int(data[theline+3])

# Read-Out file name CC
theline = searchline(readfilename,"GSMCC-files")
readfilename_CC = data[theline+1]
outfilename_CC = data[theline+2]
# Read-Out file name GSM
theline = searchline(readfilename,"GSM-files")
if theline != None:  
    readfilename_GSM = data[theline+1::2]
    outfilename_GSM = data[theline+2::2]
    gsm_files = len(outfilename_GSM)
    if len(readfilename_GSM) > gsm_files:
        readfilename_GSM = readfilename_GSM[:-1]

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
        print_twice('\n ' + running_prefix + './GSM_exe < '+inp+' > '+out)
        sp.run([running_prefix + './GSM_exe < '+inp+' > '+out], shell=True)
        end_gsm = time.time()
        time_gsm = end_gsm-start_gsm
        print_twice("Time to calculate: ",time_gsm, "s")
    #
else:
    print_twice("\nSkip GSM part, only GSMCC calculation!")
#
# Edit thresholds
print_twice("\nEdit thresholds in %s"% storage_directory)
os.chdir(storage_directory)
sp.run(['python3 Useful_Codes/GSMCC-CODE/EditThresholds.py'], shell=True)
#
print_twice("\nRunning GSMCC in %s"% gsmcc_directory)
os.chdir(gsmcc_directory)
# Experimental state to optimize with interaction.corrective.factor
expene = expene_read # Energy in MeV
expwid = expwid_read # Width in keV
numberindex = index_read
search = 'E(reference frame) :'
# Line with the corrective factor
icf_line = searchline(readfilename_CC,"CC.interaction.corrective.factor.composite(s)")+2
#
# opt = least(f,[1.048,0],diff_step=[0.02,0.001],gtol=1e-3,max_nfev=30, bounds=([0.8,1.2],[-0.1,0.1]))
if icf_type == 'COMPLEX':
    opt = newton(f,[icf_real_seed,icf_imag_seed],tol=5e-5,maxiter=30, full_output=True)
else:
    opt = newton(f,[icf_real_seed],tol=5e-5,maxiter=30, full_output=True)
#
print_twice(opt)
#
end_main = time.time()
time_main = end_main-start_main
print_twice("\n\nAll calculations lasted: ", time_main, "s")

"""
    Example of GSM+GSMCC_ICF.in file:
    _________________________________
    GSM-DIRECTORY
    /home/dassie/2024/Carbon-11_Porject/GSM-24.02/GSM_dir_2D/GSM_dir
    STORAGE-DIRECTORY
    /home/dassie/2024/Carbon-11_Porject/GSM-24.02/GSM_dir_2D/storage_11C_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30
    
    PARALLELISM
    MPI
    2
    MACHINEFILE
    machinefile
    
    CORRECTIVEFACTORS
    COMPLEX
    1.025
    0.0
    
    EXPERIMENTALVALUES
    -36.491
    5.0
    0

    GSMCC-files
    CLUSPHY_11C_CC_GSMOpt-24.08.26-11.00_Basis-24.10.24-17.00_ICF-24.10.24-17.30.in
    CLUSPHY_11C_CC_GSMOpt-24.08.26-11.00_Basis-24.10.24-17.00_ICF-24.10.24-17.30_7I2+.out

    GSM-files
    CLUSPHY_targforCC_10B_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.in
    CLUSPHY_targforCC_10B_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.out
    CLUSPHY_targforCC_7Be_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.in
    CLUSPHY_targforCC_7Be_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.out
    CLUSPHY_projforCC_protonnocore_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.in
    CLUSPHY_projforCC_protonnocore_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.out
    CLUSPHY_projforCC_alphanocore_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.in
    CLUSPHY_projforCC_alphanocore_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30.out
    _________________________________
"""