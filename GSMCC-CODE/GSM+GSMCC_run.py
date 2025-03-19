"""
    Created on August 2024 by Alan D.K. for 11C_Project

    This code run GSM and GSMC. The running information must be in
    a file called input.GSM+GSMCC_run located in the main folder of the
    GSMCC code (i.e. CC_dir/GSM+GSMCC_run.py).
    An example of the input file can be founded at the end of this code.
"""

import subprocess as sp
import time
import os

# LOG FILE
logfile = os.getcwd() + '/log.GSM+GSMCC_run.' + time.strftime( "%y.%m.%d-%H.%M", time.localtime() )

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
else:
    running_prefix = ' ./'
# Checking if we need machinefile
theline = searchline(readfilename,"MACHINEFILE:")
if theline != None:  
    machinefile_name = data[theline+1]
    running_prefix = running_prefix[:-2] + '-hostfile ' + machinefile_name + ' ./'
    
# Using experimental thresholds in GSMCC calculations
exp_thr = 0
theline = searchline(readfilename,"EXPERIMENTAL-THRESHOLDS:")
if theline != None:
    exp_thr = int(data[theline+1])

# Executable GSMCC file
calc_cc = 0
theline = searchline(readfilename,"GSMCC-exe:")
if theline != None:  
    calc_cc = 1
    running_cc = data[theline+1]
# Read-Out file name CC
theline = searchline(readfilename,"GSMCC-files:")
if theline != None:  
    readfilename_CC = data[theline+1]
    outfilename_CC = data[theline+2]
    cc_write_aux = int(data[theline+3])
    cc_write = ' ' + cc_write_aux*'>' + ' '
#
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
    print_twice("\nSkip GSM part!")
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
# RUN CC
if calc_cc == 1:
    print_twice("\nRunning GSMCC in %s"% gsmcc_directory)
    os.chdir(gsmcc_directory)
    start_gsmcc = time.time()
    print_twice('\n ' + running_prefix + running_cc + ' < ' + readfilename_CC + cc_write+outfilename_CC)
    sp.run([running_prefix + running_cc + ' < ' + readfilename_CC + cc_write+outfilename_CC], shell=True)
    end_gsmcc = time.time()
    time_gsmcc = end_gsmcc-start_gsmcc
    print_twice("Time to calculate: ",time_gsmcc, "s")
    #
    end_main = time.time()
    time_main = end_main-start_main
    print_twice("\n\nAll calculations lasted: ", time_main, "s")
else:
    print_twice("\nSkip GSMCC part!")

"""
    Example of input.GSM+GSMCC_run file:
    _________________________________
    GSM-DIRECTORY:
    /home/dassie/2024/Carbon-11_Porject/GSM-24.02/GSM_dir_2D/GSM_dir
    STORAGE-DIRECTORY:
    /home/dassie/2024/Carbon-11_Porject/GSM-24.02/GSM_dir_2D/storage_11C_GSMOpt-24.08.26-11.00_Basis-24.10.17-17.30
    
    PARALLELISM:
    MPI
    2
    MACHINEFILE:
    machinefile
    
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
    1

    GSMCC-exe:
    CC-24.11.20-MPI.x
    GSMCC-files: 1-input file; 2-output file; 3-1 for overwrite or 2 for append
    CLUSPHY_11C_CC_GSMOpt-24.08.26-11.00_Basis-24.10.24-17.00_ICF-24.10.24-17.30.in
    CLUSPHY_11C_CC_GSMOpt-24.08.26-11.00_Basis-24.10.24-17.00_ICF-24.10.24-17.30_7I2+.out
    2
    _________________________________
"""