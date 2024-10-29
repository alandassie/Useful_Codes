"""
    Created on August 2024 by Alan D.K. for 11C_Project

    This code run GSM and GSMC. You have to define the input and output files.
"""

import subprocess as sp
import numpy as np
import time
import os

# LOG FILE
logfile = os.getcwd() + '/ThrTest_run_3I2-.log'

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
readfilename = "GSM+GSMCC_run_3I2-_4.in"
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
# theline = None
# RUN GSM
if theline != None:
    print_twice("\nRunning GSM in %s"% gsm_directory)
    os.chdir(gsm_directory)
    #
    for i in range(0,gsm_files):
        inp = readfilename_GSM[i]
        out = outfilename_GSM[i]
        start_gsm = time.time()
        print_twice('\nmpirun -np 4 ./GSM_exe < '+inp+' > '+out)
        sp.run(['mpirun -np 4 ./GSM_exe < '+inp+' > '+out], shell=True)
        end_gsm = time.time()
        time_gsm = end_gsm-start_gsm
        print_twice("Time to calculate: ",time_gsm, "s")
    #
else:
    print_twice("\nSkip GSM part, only GSMCC calculation!")
#
alpha_exp = -28.296
be7_exp = -9.305
# Read alpha gs
with open(storage_directory+'/eigenvector_E_averaged_n_scat_Z2_N2_0+_0','r') as f:
    line = f.read().split('\n')
alpha_init = float(line[0].split(',')[0].split('(')[1])
#
# RUN CC Theoretical values
# print_twice("\nRunning GSMCC in %s"% gsmcc_directory)
# os.chdir(gsmcc_directory)
# start_gsmcc = time.time()
# print_twice('\nmpirun -np 4 ./CC_exe < '+readfilename_CC+' > '+outfilename_CC)
# sp.run(['mpirun -np 4 ./CC_exe < '+readfilename_CC+' >> '+outfilename_CC], shell=True)
# end_gsmcc = time.time()
# time_gsmcc = end_gsmcc-start_gsmcc
# print_twice("Time to calculate: ",time_gsmcc, "s")
#
# Prepare alpha array
alpha_array = np.arange(round(-28,1),round(alpha_exp,1)-0.1,-0.3)
# Loop over thresholds
for ene_alpha in alpha_array:
    # Edit thresholds
    print_twice("\nEdit thresholds in %s"% storage_directory)
    os.chdir(storage_directory)
    # Change file with thresholds
    aux = np.array([ene_alpha,be7_exp])
    np.savetxt('EditThresholdsVal.in',aux,fmt="%.3f")
    sp.run(['python3 EditThresholds_alt2.py'], shell=True)
    # RUN CC
    print_twice("\nRunning GSMCC in %s"% gsmcc_directory)
    os.chdir(gsmcc_directory)
    start_gsmcc = time.time()
    print_twice('\nmpirun -np 4 ./CC_exe < '+readfilename_CC+' > '+outfilename_CC)
    sp.run(['mpirun -np 4 ./CC_exe < '+readfilename_CC+' >> '+outfilename_CC], shell=True)
    end_gsmcc = time.time()
    time_gsmcc = end_gsmcc-start_gsmcc
    print_twice("Time to calculate: ",time_gsmcc, "s")
#
#
end_main = time.time()
time_main = end_main-start_main
print_twice("\n\nAll calculations lasted: ", time_main, "s")
