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
print_twice('START CONTOUR TEST\n')
theline = searchline(calcfilename,"N_CLUSTERS:")
n_clusters = int(data[theline+1])
print_twice('Number of clusters/projectiles in the test: %s'% n_clusters)
# NAME OF EACH CLUSTER
theline = searchline(calcfilename,"NAMES_CLUSTERS:")
print_twice('Names:')
names_clusters = []
for i in range(n_clusters):
    names_clusters.append( data[theline +1+i] )
    print_twice(names_clusters[i])
# PARTIAL WAVE FOR EACH CLUSTER
n_pw_clusters = []
pw_clusters = []
for partition in names_clusters:
    print_twice('Partial waves to modify in the %s cluster/projectile'% partition)
    theline = searchline(calcfilename,"%s_PW:"% partition)
    n_pw = int( data[theline + 1] )
    n_pw_clusters.append( n_pw )
    aux = []
    for i in range(n_pw):
        aux.append( data[theline +2+i] )
    print_twice(" ".join(aux))
    pw_clusters.append( aux )
#
#
# Read contours from file first:
# Search the contour line
# theline = searchline(readfilename_CC,".peak(fm^(-1))")
# # Prepare arrays
# kpeak = []
# kmiddle = []
# kmax = []
# # Read information from the specific mass particion and partial width
# for i in range(n_clusters):
#     partition = names_clusters[i]
#     aux_kpeak = []
#     aux_kmiddle = []
#     aux_kmax = []
#     #
#     j = 0
#     while j >= 0:
#         aux = data[theline + 1 + j]
#         if aux.split() == []:
#             contourlines = j-1
#             j = -1
#             continue
#         aux_p = aux.split()[0]
#         if aux_p == partition:
#             aux_pw = aux.split()[1]
#             if aux_pw in pw_clusters:
#                 aux_kpeak.append(eval( aux.split()[4] ))
#                 aux_kmiddle.append(eval( aux.split()[5] ))
#                 aux_kmax.append(eval( aux.split()[6] ))
#         elif aux_p in partial_wave: # If is not a partition, could be directly a partial wave for the nucleon cases
#             aux_kpeak.append(eval( aux.split()[4] ))
#             aux_kmiddle.append(eval( aux.split()[5] ))
#             aux_kmax.append(eval( aux.split()[6] ))
#         j += 1
#     kpeak.append(aux_kpeak)
#     kmiddle.append(aux_kmiddle)
#     kmax.append(aux_kmax)
#
# Kpeak, Kmiddle OR Kmax FOR MODIFICATION:
theline = searchline(calcfilename,"K_TO_MODIFY:")
kmod = int(data[theline + 1]) # 1-Kpeak, 2-Kmiddle OR 3-Kmax
ktype = data[theline + 2] # REAL, IMAG
# Read
if kmod == 1:
    print_twice('\n The %s kpeak will be tested between %s and %s'% (ktype, float(data[theline + 3]), float(data[theline + 4])))
    kpeak_array = np.arange( float(data[theline + 3]), float(data[theline + 4])+float(data[theline + 5])/2, float(data[theline + 5]) )
if kmod == 2:
    print_twice('\n The %s kmiddle will be tested between %s and %s'% (ktype, float(data[theline + 3]), float(data[theline + 4])))
    kmiddle_array = np.arange( float(data[theline + 3]), float(data[theline + 4])+float(data[theline + 5])/2, float(data[theline + 5]) )
if kmod == 3:
    print_twice('\n The %s kmax will be tested between %s and %s'% (ktype, float(data[theline + 3]), float(data[theline + 4])))
    kmax_array = np.arange( float(data[theline + 3]), float(data[theline + 4])+float(data[theline + 5])/2, float(data[theline + 5]) )
            
# MODIFY PARTIAL WAVES CONTOUR OF EACH CLUSTER AT THE SAME TIME: ONLY OPTION FOR THE MOMENT
# theline = searchline(calcfilename,"SAME_TIME_CALC:")
# same_time = int( data[theline + 1] )
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
with open(readfilename_CC, 'r') as readfile:
    data = readfile.read().split('\n')
    #
theline = searchline(readfilename_CC,".peak(fm^(-1))")
# Start calculation loop
if kmod == 1:
    k=0
    for x1 in kpeak_array:
        j = 0
        while j >= 0:
            aux = data[theline + 1 + j].split()
            if aux == []:
                j = -1
                continue
            aux_p = aux[0]
            if aux_p in partition:
                aux_pw = aux[1]
                if any(aux_pw in t for t in pw_clusters):
                    if ktype == 'REAL':
                        auxi = list(eval(aux[4]))
                        auxi[0] = x1
                        aux[4] = "".join(repr(tuple(auxi)).split())
                    if ktype == 'IMAG':
                        auxi = list(eval(aux[4]))
                        aux2 = list(eval(aux[5]))
                        auxi[1] = x1
                        aux2[1] = x1 # Kpeak and Kmiddle must have the same imaginary value
                        aux[4] = "".join(repr(tuple(auxi)).split())
                        aux[5] = "".join(repr(tuple(aux2)).split())
            if any(aux_p in t for t in pw_clusters):
                if ktype == 'REAL':
                    auxi = list(eval(aux[4]))
                    auxi[0] = x1
                    aux[4] = "".join(repr(tuple(auxi)).split())
                if ktype == 'IMAG':
                    auxi = list(eval(aux[4]))
                    aux2 = list(eval(aux[5]))
                    auxi[1] = x1
                    aux2[1] = x1 # Kpeak and Kmiddle must have the same imaginary value
                    aux[4] = "".join(repr(tuple(auxi)).split())
                    aux[5] = "".join(repr(tuple(aux2)).split())
            data[theline + 1 + j] = "  ".join(aux)
            j += 1
            #
        inputfile_aux = '\n'.join(data)
        with open(readfilename_CC,'w') as fileout:
            fileout.write(inputfile_aux)
        # Runnning the code
        start_gsmcc = time.time()
        outfilename_CC_j = logname + '/' + outfilename_CC[:-4] + '_%s.out'% (k)
        print_twice('\n ' + running_prefix + running_cc + ' < ' + readfilename_CC + cc_write+outfilename_CC_j)
        sp.run([running_prefix + running_cc + ' < ' + readfilename_CC + cc_write+outfilename_CC_j], shell=True)
        end_gsmcc = time.time()
        time_gsmcc = end_gsmcc-start_gsmcc
        print_twice("Time to calculate: ",time_gsmcc, "s")
        k += 1
            
if kmod == 2:
    k = 0
    for x1 in kmiddle_array:
        datain = []
        j = 0
        while j >= 0:
            aux = data[theline + 1 + j].split()
            if aux == []:
                j = -1
                continue
            aux_p = aux[0]
            if aux_p in partition:
                aux_pw = aux[1]
                if any(aux_pw in t for t in pw_clusters):
                    if ktype == 'REAL':
                        auxi = list(eval(aux[5]))
                        auxi[0] = x1
                        aux[5] = "".join(repr(tuple(auxi)).split())
                    if ktype == 'IMAG':
                        auxi = list(eval(aux[5]))
                        aux2 = list(eval(aux[4]))
                        auxi[1] = x1
                        aux2[1] = x1 # Kpeak and Kmiddle must have the same imaginary value
                        aux[5] = "".join(repr(tuple(auxi)).split())
                        aux[4] = "".join(repr(tuple(aux2)).split())
            if any(aux_p in t for t in pw_clusters):
                if ktype == 'REAL':
                    auxi = list(eval(aux[5]))
                    auxi[0] = x1
                    aux[5] = "".join(repr(tuple(auxi)).split())
                if ktype == 'IMAG':
                    auxi = list(eval(aux[5]))
                    aux2 = list(eval(aux[4]))
                    auxi[1] = x1
                    aux2[1] = x1 # Kpeak and Kmiddle must have the same imaginary value
                    aux[5] = "".join(repr(tuple(auxi)).split())
                    aux[4] = "".join(repr(tuple(aux2)).split())
            data[theline + 1 + j] = "  ".join(aux)
            j += 1
            #
        inputfile_aux = '\n'.join(data)
        with open(readfilename_CC,'w') as fileout:
            fileout.write(inputfile_aux)
        # Runnning the code
        start_gsmcc = time.time()
        outfilename_CC_j = logname + '/' + outfilename_CC[:-4] + '_%s.out'% (k)
        print_twice('\n ' + running_prefix + running_cc + ' < ' + readfilename_CC + cc_write+outfilename_CC_j)
        sp.run([running_prefix + running_cc + ' < ' + readfilename_CC + cc_write+outfilename_CC_j], shell=True)
        end_gsmcc = time.time()
        time_gsmcc = end_gsmcc-start_gsmcc
        print_twice("Time to calculate: ",time_gsmcc, "s")
        k += 1

if kmod == 3:
    k = 0
    for x1 in kmax_array:
        j = 0
        while j >= 0:
            aux = data[theline + 1 + j].split()
            if aux == []:
                j = -1
                continue
            aux_p = aux[0]
            if aux_p in partition:
                aux_pw = aux[1]
                if any(aux_pw in t for t in pw_clusters):
                    if ktype == 'REAL':
                        aux[6] = str(x1)
            if any(aux_p in t for t in pw_clusters):
                if ktype == 'REAL':
                    aux[6] = str(x1)
            data[theline + 1 + j] = "  ".join(aux)
            j += 1
            #
        inputfile_aux = '\n'.join(data)
        with open(readfilename_CC,'w') as fileout:
            fileout.write(inputfile_aux)
        # Runnning the code
        start_gsmcc = time.time()
        outfilename_CC_j = logname + '/' + outfilename_CC[:-4] + '_%s.out'% (k)
        print_twice('\n ' + running_prefix + running_cc + ' < ' + readfilename_CC + cc_write + outfilename_CC_j)
        sp.run([running_prefix + running_cc + ' < ' + readfilename_CC + cc_write+outfilename_CC_j], shell=True)
        end_gsmcc = time.time()
        time_gsmcc = end_gsmcc-start_gsmcc
        print_twice("Time to calculate: ",time_gsmcc, "s")
        k += 1

#
end_main = time.time()
time_main = end_main-start_main
print_twice("\n\nAll calculations lasted: ", time_main, "s")

# # Saving startgin point
# with open(readfilename,'r') as gsmin:
#     startingpoint_lines = gsmin.read().split('\n')
# theline = searchline(readfilename,"N.K.CM.max")
# i = 1
# datain = []
# while i > 0:
#     aux_line = startingpoint_lines[theline+i]
#     if aux_line.split() == []:
#         break
#     aux_proj = aux_line.split()[0]
#     if aux_proj in names_clusters:
#         datain.append(aux_line)
#     i += 1
# input_n_contour = i

# Calculations for K_Peak
# Kmiddle, Kmax and N are fixed from input 

"""
    Example of input.ContourTest file:
    _________________________________
    N_CLUSTERS: # Number of clusters/projectiles
    1

    NAMES_CLUSTERS: # Name of each cluster/projectile
    proton

    proton_PW: # Partial waves for modify contour of each cluster/projectile; 1-number of pw, 2-then name of each (same order as GSMCC file)
    5
    s1/2
    p3/2
    p1/2
    d5/2
    d3/2

    K_TO_MODIFY: # The k to modify, for the moment is the same for each cluster and partial wave. 1-KPEAK; 2-KMIDDLE; 3-KMAX; 2-then, REAL or IMAG part; 3-min; 4-max; 5-step
    3
    REAL
    1.6
    5
    0.2
    _________________________________
"""