"""
    Created on August 2024 by Alan D.K. for 11C_Project

    This code run multiples GSM and/or GSMCC calculations in order to calculate 
    the GSMCC spectrum with different Basis.parameters.
    
    For reading file, the information about working directories and files
    is the same as in the code GSM+GSMCC_run.py.
    A particular file called BasisTest_alt.in is used to define the parameters to
    test, the range of test and the step for each one. An example can be found at the end of the code.
"""

import numpy as np
import subprocess as sp
import os
import time
import json


# LOG FILE
logfile = os.getcwd() + '/BasisTest_log.' + time.strftime( "%y.%m.%d-%H.%M", time.localtime() )

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
# Idea from https://stackoverflow.com/questions/176918/how-to-find-the-index-for-a-given-item-in-a-list
def indexof (obj, elem, offset=0):
    if elem in obj[offset:]:
        return offset + obj[offset:].index(elem)
    return -1
# .-

# INPUT FILE
readfilename = os.getcwd() + '/GSM+GSMCC_run.in'
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
    running_prefix = 'mpirun -np ' + parallelism_nodes + ' '
else:
    running_prefix = ' '
# Checking if we need machinefile
theline = searchline(readfilename,"MACHINEFILE:")
if theline != None:  
    machinefile_name = data[theline+1]
    running_prefix = running_prefix + '-hostfile ' + machinefile_name + ' '

# Executable GSMCC file
theline = searchline(readfilename,"GSMCC-exe:")
running_cc = data[theline+1]
# Read-Out file name CC
theline = searchline(readfilename,"GSMCC-files:")
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

# CALCULATING FILE
calcfilename = os.getcwd() + '/BasisTest_alt.in'
with open(calcfilename, 'r') as readfile:
    data = readfile.read().split('\n')
numfile = int(data[1])

# Start calculation
start_main = time.time()
# RUN GSM
theline = searchline(readfilename,"GSM-files:")
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
print_twice("\nEdit thresholds in %s"% storage_directory)
os.chdir(storage_directory)
sp.run(['python3 Useful_Codes/GSMCC-CODE/EditThresholds.py'], shell=True)
#
print_twice("\nRunning GSMCC in %s"% gsmcc_directory)
os.chdir(gsmcc_directory)
# Calculating
print_twice('Start calculations')
# Start calculations
for j in range(0,numfile):
    # Open GSMCC input file
    readfilename_CC_j = readfilename_CC[:-3] + '_%s.in'% (j+10)
    outfilename_CC_j = outfilename_CC[:-4] + '_%s.out'% (j+10)
    # Runnning the code
    start_gsmcc = time.time()
    print_twice('\n ' + running_prefix + running_cc + ' < ' + readfilename_CC_j + cc_write+outfilename_CC_j)
    sp.run([running_prefix + running_cc + ' < ' + readfilename_CC_j + cc_write+outfilename_CC_j], shell=True)
    end_gsmcc = time.time()
    time_gsmcc = end_gsmcc-start_gsmcc
    print_twice("Time to calculate: ",time_gsmcc, "s")
    # if j == 16 :
    #     readfilename_CC_j = readfilename_CC[:-3] + '_%sa.in'% (j+10)
    #     outfilename_CC_j = outfilename_CC[:-4] + '_%sa.out'% (j+10)
    #     # Runnning the code
    #     start_gsmcc = time.time()
    #     print_twice('\n ' + running_prefix + running_cc + ' < ' + readfilename_CC_j + cc_write+outfilename_CC_j)
    #     sp.run([running_prefix + running_cc + ' < ' + readfilename_CC_j + cc_write+outfilename_CC_j], shell=True)
    #     end_gsmcc = time.time()
    #     time_gsmcc = end_gsmcc-start_gsmcc
    #     print_twice("Time to calculate: ",time_gsmcc, "s")
        

#
end_main = time.time()
time_main = end_main-start_main
print_twice("\n\nAll calculations lasted: ", time_main, "s")

"""
    Example of BasisTest.in file:
    _________________________________
    # Remove the lines that you are not using, this example is with the complete test
    # Number of input files to use
    40
    _________________________________
"""