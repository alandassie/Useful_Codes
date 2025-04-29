"""
    Created on June 2024 by Alan D.K.

    This code run multiples GSM and/or GSMCC calculations in order to calculate 
    the GSMCC spectrum with different contour.
    
    For reading file, the information about working directories and files
    is the same as in the code GSM+GSMCC_run.py.
    A particular file called input.ContourTest is used to define the parameters of the test.
    An example can be found at the end of the code.
"""

import numpy as np
import subprocess as sp
import os
import time

# LOG FILE
logfile = os.getcwd() + '/log.ContourTest.' + time.strftime( "%y.%m.%d-%H.%M", time.localtime() )
# LOG NAME
logname = 'log.WD_ContourTest.' + time.strftime( "%y.%m.%d-%H.%M", time.localtime() )
# LOG FOLDER
logfolder = os.getcwd() + '/' + logname 
if not os.path.exists(logfolder):
    os.makedirs(logfolder)
else:
    for filename in os.listdir(logfolder):
        # remove the files inside
        os.remove(f"{logfolder}/{filename}")

# Declaration of funcitons
def print_twice(*args,**kwargs): # Allow to print in file and in terminal at the same line
    print(*args,**kwargs)
    with open(optout,"a") as f:  # appends to file and closes it when finished
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
# storage directory
theline = searchline(readfilename,"STORAGE-DIRECTORY:")
storage_directory = data[theline+1]
# Checking if it is a MPI or OPENMP/secuential calculation
theline = searchline(readfilename,"PARALLELISM:")
parallelism_type = data[theline+1]
parallelism_nodes = data[theline+2]
if parallelism_type == 'MPI':
    running_prefix = 'mpirun -np ' + parallelism_nodes + ' '
else:
    running_prefix = './'
# Checking if we need machinefile
theline = searchline(readfilename,"MACHINEFILE:")
if theline != None:  
    machinefile_name = data[theline+1]
    running_prefix = running_prefix + '-hostfile ' + machinefile_name + ' '
    
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
calcfilename = os.getcwd() + '/input.ContourTest'
with open(calcfilename, 'r') as readfile:
    data = readfile.read().split('\n')

# IDENTIFY NUMBER OF CLUSTER CONTOUR TO CHANGE
theline = searchline(calcfilename,"N_CLUSTERS:")
n_clusters = int(data[theline+1])
# NAME OF EACH CLUSTER
theline = searchline(calcfilename,"NAMES_CLUSTERS:")
names_clusters = []
for i in range(n_clusters):
    names_clusters.append( data[theline +1+i] )
# PARTIAL WAVE FOR EACH CLUSTER
n_pw_clusters = []
pw_clusters = []
for partition in names_clusters:
    theline = searchline(calcfilename,"%s_PW:"% partition)
    n_pw = int( data[theline + 1] )
    n_pw_clusters.append( n_pw )
    aux = []
    for i in range(n_pw):
        aux.append( data[theline +2+i] )
    pw_clusters.append( aux )
# Kpeak, Kmiddle OR Kmax FOR MODIFICATION:
theline = searchline(calcfilename,"K_TO_MODIFY:")
kmod = int(data[theline + 1]) # 1-Kpeak, 2-Kmiddle OR 3-Kmax
ktype = data[theline + 2] # REAL, IMAG
#
kpeak_real_min = []
kpeak_real_max = []
kpeak_imag_min = []
kpeak_imag_max = []
#
for partition in names_clusters:
    theline = searchline(calcfilename,"%s_K:"% partition)
    # if kmod == 1:
    #     if ktype == 'REAL':
            

# MODIFY PARTIAL WAVES CONTOUR OF EACH CLUSTER AT THE SAME TIME
theline = searchline(calcfilename,"SAME_TIME_CALC:")
same_time = int( data[theline + 1] )


# Define Kpeak, Kmiddle and Kmax limits
kpeak_real_min = 0.2
kpeak_imag_min = 0.06
kpeak_real_max = 0.2
# kpeak_imag_max will be allways between kpeak_imag_min and kpeak_real_max
kpeak_step = 0.02
#
kmiddle_real_min = 1
kmiddle_real_max = 1
# kmiddle_imag and kmiddle_imag will be allways equal to kpeak_imag
kmiddle_step = 0.1
#
kmax_min = 1.5
kmax_max = 3
kmax_step = 0.5
#
# Define number of scattering states
Nmin = 20
Nmax = 25
Nstep = 5
#
# Start calculation loop
# kpeak_real will be the outher loop
kpeak_real = np.arange(kpeak_real_min, kpeak_real_max+kpeak_step/2, kpeak_step)
kmiddle_real = np.arange(kmiddle_real_min, kmiddle_real_max+kmiddle_step/2, kmiddle_step)
N_array = np.arange(Nmin, Nmax + 1, Nstep)
# N_array = np.array([10,15,25])
for x1 in kpeak_real:
    for x2 in kmiddle_real:
        # Build kpeak_imag array
        kpeak_imag = np.arange(kpeak_imag_min, 0.14+kpeak_step/2, kpeak_step)
        # kpeak_imag = np.array([0.05])
        for readfilename in readfilename_array:
            for y1 in kpeak_imag:
                y2 = y1
                for N in N_array:
                    #
                    with open(readfilename, 'r') as readfile:
                        data = readfile.read().split('\n')
                        #
                    # Search the contour line
                    theline = searchline(readfilename,"N.K.CM.max")
                    datain = ["proton         s1/2                   1                              1                   ({0:.2f},-{1:.3f})        ({2:.2f},-{3:.3f})         2.0                 {4:2d}          {4:2d}            {4:2d} ".format(x1,y1,x2,y2,N),
                              "proton         p1/2                   1                              0                   ({0:.2f},-{1:.3f})        ({2:.2f},-{3:.3f})         2.0                 {4:2d}          {4:2d}            {4:2d} ".format(x1,y1,x2,y2,N),
                              "proton         p3/2                   1                              0                   ({0:.2f},-{1:.3f})        ({2:.2f},-{3:.3f})         2.0                 {4:2d}          {4:2d}            {4:2d} ".format(x1,y1,x2,y2,N),
                              "proton         d3/2                   0                              0                   ({0:.2f},-{1:.3f})        ({2:.2f},-{3:.3f})         2.0                 {4:2d}          {4:2d}            {4:2d} ".format(x1,y1,x2,y2,N),
                              "proton         d5/2                   0                              0                   ({0:.2f},-{1:.3f})        ({2:.2f},-{3:.3f})         2.0                 {4:2d}          {4:2d}            {4:2d} ".format(x1,y1,x2,y2,N),
                              "proton         f5/2                   0                              0                   ({0:.2f},-{1:.3f})        ({2:.2f},-{3:.3f})         2.0                 {4:2d}          {4:2d}            {4:2d} ".format(x1,y1,x2,y2,N),
                              "proton         f7/2                   0                              0                   ({0:.2f},-{1:.3f})        ({2:.2f},-{3:.3f})         2.0                 {4:2d}          {4:2d}            {4:2d} ".format(x1,y1,x2,y2,N),
                              "proton         g7/2                   0                              0                   ({0:.2f},-{1:.3f})        ({2:.2f},-{3:.3f})         2.0                 {4:2d}          {4:2d}            {4:2d} ".format(x1,y1,x2,y2,N),
                              "proton         g9/2                   0                              0                   ({0:.2f},-{1:.3f})        ({2:.2f},-{3:.3f})         2.0                 {4:2d}          {4:2d}            {4:2d} ".format(x1,y1,x2,y2,N),
                              " ",
                              "alpha          S                      4                              0                   ({0:.2f},-{1:.3f})        ({2:.2f},-{3:.3f})         2.0                 {4:2d}          {4:2d}            {4:2d} ".format(x1,y1,x2,y2,N),
                              "alpha          P                      4                              0                   ({0:.2f},-{1:.3f})        ({2:.2f},-{3:.3f})         2.0                 {4:2d}          {4:2d}            {4:2d} ".format(x1,y1,x2,y2,N),
                              "alpha          D                      3                              0                   ({0:.2f},-{1:.3f})        ({2:.2f},-{3:.3f})         2.0                 {4:2d}          {4:2d}            {4:2d} ".format(x1,y1,x2,y2,N),
                              "alpha          F                      2                              0                   ({0:.2f},-{1:.3f})        ({2:.2f},-{3:.3f})         2.0                 {4:2d}          {4:2d}            {4:2d} ".format(x1,y1,x2,y2,N),
                              "alpha          G                      1                              0                   ({0:.2f},-{1:.3f})        ({2:.2f},-{3:.3f})         2.0                 {4:2d}          {4:2d}            {4:2d} ".format(x1,y1,x2,y2,N)]
                    data[theline+1:theline+len(datain)+1] = datain
                    inputfile_aux = '\n'.join(data)
                    with open(readfilename,'w') as fileout:
                        fileout.write(inputfile_aux)
                    #
                    sp.run(['mpirun -np 4 ./CC_exe < '+readfilename+' >> IN2P3_11C_CC_GSMOpt-24.07.24-10.00_Basis-24.08.07-19.00_3I2-.out'], shell=True)
                    #


# Start Calculations
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
# Saving startgin point
with open(readfilename,'r') as gsmin:
    startingpoint_lines = gsmin.read().split('\n')
theline = searchline(readfilename,"N.K.CM.max")
i = 1
datain = []
while i > 0:
    aux_line = startingpoint_lines[theline+i]
    if aux_line.split() == []:
        break
    aux_proj = aux_line.split()[0]
    if aux_proj in names_clusters:
        datain.append(aux_line)
    i += 1
input_n_contour = i

# Calculations for K_Peak
# Kmiddle, Kmax and N are fixed from input 