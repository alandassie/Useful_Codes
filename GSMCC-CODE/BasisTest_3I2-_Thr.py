"""
    Created on August 2024 by Alan D.K. for 11C_Project

    This code run multiples GSM and GSMCC calculations in order to calculate 
    the GSMCC spectrum with different Basis.parameters
"""

import numpy as np
import subprocess as sp
import os
import time


# LOG FILE
logfile = os.getcwd() + '/BasisTest.log'

# GSMCC directory
gsmcc_directory = os.getcwd()
# GSM directory
gsm_directory = "/pbs/home/a/adassie/Carbon-11_Project/GSM-MPI/GSM_dir_2D/GSM_dir"
# storage directory
storage_directory = "/pbs/throng/ganil/adassie/storage_GSM-24.07.24-10.00_BasisTest2"

# Output file name
# now = datetime.now()
# optout = 'conttest%s.out'% (now.strftime("%y.%m.%d-%H:%M"))
# # Erease output file
# with open(optout,'w') as output:
#     pass
# INPUT FILES
# Read file name CC
readfilenameCC_array = ['IN2P3_11C_CC_GSMOpt-24.07.24-10.00_Basis-24.08.13-09.00-%s_3I2-.in'% i for i in range(1,10)]
# Read file name GSM 11C
readfilename11C_array = ['IN2P3_11C_GSMOpt-24.07.24-10.00_Basis-24.08.13-09.00_%s.in'% i for i in range(1,10)]
# Read file name GSM 7Be
readfilename7Be_array = ['IN2P3_7Be_targforCC_GSMOpt-24.07.24-10.00_Basis-24.08.13-09.00_%s.in'% i for i in range(1,10)]
# Read file name GSM 10B
readfilename10B_array = ['IN2P3_10B_targforCC_GSMOpt-24.07.24-10.00_Basis-24.08.13-09.00_%s.in'% i for i in range(1,10)]
# Read file name GSM alpha
readfilenamealpha = 'IN2P3_alphanocore_projforCC_GSMOpt-24.07.24-10.00_Basis-24.08.13-09.00.in'
# Read file name GSM proton
readfilenameproton = 'IN2P3_protonnocore_projforCC_GSMOpt-24.07.24-10.00_Basis-24.08.13-09.00.in'
# OUTPUT FILES
outfilenameCC_array = ['IN2P3_11C_CC_GSMOpt-24.07.24-10.00_Basis-24.08.13-09.00-%s_3I2-_ExpThreshold.out'% i for i in range(1,10)]
# out file name GSM 11C
outfilename11C_array = ['IN2P3_11C_GSMOpt-24.07.24-10.00_Basis-24.08.13-09.00_%s.out'% i for i in range(1,10)]
# out file name GSM 7Be
outfilename7Be_array = ['IN2P3_7Be_targforCC_GSMOpt-24.07.24-10.00_Basis-24.08.13-09.00_%s.out'% i for i in range(1,10)]
# out file name GSM 10B
outfilename10B_array = ['IN2P3_10B_targforCC_GSMOpt-24.07.24-10.00_Basis-24.08.13-09.00_%s.out'% i for i in range(1,10)]
# out file name GSM alpha
# outfilenamealpha_array = ['IN2P3_alphanocore_projforCC_GSMOpt-24.07.24-10.00_Basis-24.08.08-11.40.out',
#                           'IN2P3_alphanocore_projforCC_GSMOpt-24.07.24-10.00_Basis-24.08.08-12.00.out',
#                           'IN2P3_alphanocore_projforCC_GSMOpt-24.07.24-10.00_Basis-24.08.08-12.20.out']
outfilenamealpha = 'IN2P3_alphanocore_projforCC_GSMOpt-24.07.24-10.00_Basis-24.08.12-16.00.out'
# out file name GSM proton
# outfilenameproton_array = ['IN2P3_protonnocore_projforCC_GSMOpt-24.07.24-10.00_Basis-24.08.08-11.40.out',
#                            'IN2P3_protonnocore_projforCC_GSMOpt-24.07.24-10.00_Basis-24.08.08-12.00.out',
#                            'IN2P3_protonnocore_projforCC_GSMOpt-24.07.24-10.00_Basis-24.08.08-12.20.out']
outfilenameproton = 'IN2P3_protonnocore_projforCC_GSMOpt-24.07.24-10.00_Basis-24.08.12-16.00.out'

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

# LOGILE
erease_output_file()

# Define limits and step
numberoffiles = len(readfilenameCC_array)

# Start calculation loop
start_main = time.time()
for i in range(0,numberoffiles):
    start_one = time.time()
    # RUN GSM
    os.chdir(gsm_directory)
    print_twice("\nRunning GSM in %s"% gsm_directory)
    #
    readfilename7Be = readfilename7Be_array[i]
    outfilename7Be = outfilename7Be_array[i]
    start_gsm = time.time()
    print_twice('\nmpirun -np 4 ./GSM_exe < '+readfilename7Be+' > '+outfilename7Be)
    sp.run(['mpirun -np 4 ./GSM_exe < '+readfilename7Be+' > '+outfilename7Be], shell=True)
    end_gsm = time.time()
    time_gsm = end_gsm-start_gsm
    print_twice("Time to calculate: ",time_gsm, "s")
    #
    readfilename10B = readfilename10B_array[i]
    outfilename10B = outfilename10B_array[i]
    start_gsm = time.time()
    print_twice('\nmpirun -np 4 ./GSM_exe < '+readfilename10B+' > '+outfilename10B)
    sp.run(['mpirun -np 4 ./GSM_exe < '+readfilename10B+' > '+outfilename10B], shell=True)
    end_gsm = time.time()
    time_gsm = end_gsm-start_gsm
    print_twice("Time to calculate: ",time_gsm, "s")
    # 
    readfilename11C = readfilename11C_array[i]
    outfilename11C = outfilename11C_array[i]
    start_gsm = time.time()
    print_twice('\nmpirun -np 4 ./GSM_exe < '+readfilename11C+' > '+outfilename11C)
    sp.run(['mpirun -np 4 ./GSM_exe < '+readfilename11C+' > '+outfilename11C], shell=True)
    end_gsm = time.time()
    time_gsm = end_gsm-start_gsm
    print_twice("Time to calculate: ",time_gsm, "s")
    #
    # readfilenamealpha = readfilenamealpha_array[i]
    # outfilenamealpha = outfilenamealpha_array[i]
    start_gsm = time.time()
    print_twice('\nmpirun -np 4 ./GSM_exe < '+readfilenamealpha+' > '+outfilenamealpha)
    sp.run(['mpirun -np 4 ./GSM_exe < '+readfilenamealpha+' > '+outfilenamealpha], shell=True)
    end_gsm = time.time()
    time_gsm = end_gsm-start_gsm
    print_twice("Time to calculate: ",time_gsm, "s")
    #
    # readfilenameproton = readfilenameproton_array[i]
    # outfilenameproton = outfilenameproton_array[i]
    start_gsm = time.time()
    print_twice('\nmpirun -np 4 ./GSM_exe < '+readfilenameproton+' > '+outfilenameproton)
    sp.run(['mpirun -np 4 ./GSM_exe < '+readfilenameproton+' > '+outfilenameproton], shell=True)
    end_gsm = time.time()
    time_gsm = end_gsm-start_gsm
    print_twice("Time to calculate: ",time_gsm, "s")
    #
    # Edit thresholds
    os.chdir(storage_directory)
    print_twice("\nEdit thresholds in %s"% storage_directory)
    sp.run(['python3 EditThresholds_alt.py'], shell=True)
    #
    # RUN CC
    os.chdir(gsmcc_directory)
    print_twice("\nRunning GSMCC in %s"% gsmcc_directory)
    readfilenameCC = readfilenameCC_array[i]
    outfilenameCC = outfilenameCC_array[i]
    start_gsmcc = time.time()
    print_twice('\nmpirun -np 4 ./CC_exe < '+readfilenameCC+' > '+outfilenameCC)
    sp.run(['mpirun -np 4 ./CC_exe < '+readfilenameCC+' > '+outfilenameCC], shell=True)
    end_gsmcc = time.time()
    time_gsmcc = end_gsmcc-start_gsmcc
    print_twice("Time to calculate: ",time_gsmcc, "s")
    #
    end_one = time.time()
    time_one = end_one-start_one
    print_twice("\nNumber %s calculations lasted: "% i, time_one, "s")
    
end_main = time.time()
time_main = end_main-start_main
print_twice("\n\nAll calculations lasted: ", time_main, "s")
