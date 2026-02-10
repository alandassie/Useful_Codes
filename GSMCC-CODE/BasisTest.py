"""
    Created on August 2024 by Alan D.K. for 11C_Project

    This code run multiples GSM and/or GSMCC calculations in order to calculate 
    the GSMCC spectrum with different Basis.parameters.
    
    For reading file, the information about working directories and files
    is the same as in the code GSM+GSMCC_run.py.
    A particular file called BasisTest.in is used to define the parameters to
    test, the range of test and the step for each one. An example can be found at the end of the code.
"""

import numpy as np
import subprocess as sp
import os
import time
import ast


# LOG NAME
logname = 'log.WD_BasisTest.' + time.strftime( "%y.%m.%d-%H.%M", time.localtime() )
# LOG FOLDER
logfolder = os.getcwd() + '/' + logname 
if not os.path.exists(logfolder):
    os.makedirs(logfolder)
else:
    for filename in os.listdir(logfolder):
        # remove the files inside
        os.remove(f"{logfolder}/{filename}")
# LOG FILE
logfile = logfolder + '/log.BasisTest.' + time.strftime( "%y.%m.%d-%H.%M", time.localtime() )

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
#
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
# .-
# Idea from https://thelearninghub.in/get-list-of-indexes-of-a-repeated-item-in-python-list/
def indexof (my_list, key):
    index_val = [] 
    i=0 #counter
    while i != len(my_list):
        try:
            ind = my_list[i:len(my_list)].index(key) #checking the index 
            index_val.append(ind+i)
            i = ind+1+i
        except ValueError:
            i = len(my_list)
    return index_val
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
    running_prefix = 'mpirun -np ' + parallelism_nodes + ' '
elif parallelism_type == 'OPENMP':
    running_prefix = ' ./'
elif parallelism_type == 'SLURM':
    running_prefix = 'srun -N ' + parallelism_nodes + ' '
else:
    print_twice('Parallelism must be MPI, OPENMP or SRUN')
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
calcfilename = os.getcwd() + '/input.BasisTest'
with open(calcfilename, 'r') as readfile:
    data = readfile.read().split('\n')

# Generalization to all possible projectiles at the same time
# Preparing input arrays
theline = searchline(calcfilename,"PROJECTILE_TYPE:")
if theline != None:
    projectile_type = data[theline+1] # Only one projectile at a time, name as in GSMCC input file = "proton", "neutron", "deuteron", "triton", "3He", "alpha"
    projectile_type_list = ['proton', 'neutron', 'deuteron', 'triton', '3He', 'alpha']
    if projectile_type not in projectile_type_list:
        print_twice('Projectile type must be one of the following: ', projectile_type_list)
        exit()
    print_twice('\n Doing Woods-Saxon+SO test for %s'% projectile_type)
    parameters_type = int(data[theline+2]) # Parameter to change: 0 - diffuseness (a0), 1 - radius (R0), 2 - WS depth (VO), 3 - spin-orbit depth (VSO)
    parameters_list = ['diffuseness (a0)', 'radius (R0)', 'WS depth (VO)', 'spin-orbit depth (VSO)']
    if parameters_type < 0 or parameters_type > 3:
        print_twice('Parameter to change must be 0, 1, 2 or 3')
        exit()
    print_twice('  Changing parameter: %s'% parameters_list[parameters_type])
    theline = searchline(calcfilename,"PARTIAL_WAVES:")
    projectile_partialwaves_type = data[theline+1]
    if projectile_partialwaves_type == 'ALL':
        print_twice('  All WS partial waves at the same time!')
        theline = searchline(calcfilename,"STRENGTHWS")
        projectile_lwave_n = 0
        projectile_step = float(data[theline+1])
        projectile_n = int(data[theline+2])
        print_twice('    Step: %s, Times: %s', (projectile_step, projectile_n))
    else:
        print_twice('  Selected specific partial waves:')
        projectile_lwave = ast.literal_eval(projectile_partialwaves_type)
        print_twice('    ', projectile_lwave)
        projectile_lwave_n = len(projectile_lwave)
        theline = searchline(calcfilename,"STRENGTHWS")
        projectile_step = float(data[theline+1])
        projectile_n = int(data[theline+2])
        print_twice('    Step: %s, Times: %s', (projectile_step, projectile_n))


# Saving WF
theline = searchline(calcfilename,"SAVE_WF:")
if theline != None:
    savingwf_info = data[theline+1]
    if savingwf_info == 'NO':
        print_twice('\nINFO: The channels WF of each calculation will not be saved')
    else:
        theline = searchline(calcfilename,"WF_JPI_INDEX:")
        n_wf = int(data[theline+1])
        jpi_wf = []
        index_wf = []
        for i in range(0,n_wf):
            jpi_wf.append( data[theline+2*(i+1)] )
            index_wf.append( data[theline+2*(i+1)+1] )

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
with open(readfilename_CC,'r') as gsmin:
    startingpoint_lines = gsmin.read().split('\n')
#


# Calculating projectile WS
theline = searchline(calcfilename,"PROJECTILE_TYPE:")
if theline != None:
    projectile_type = data[theline+1]
    print_twice('Start GSMCC target+%s calculations for the parameter %s'% (projectile_type, parameters_list[parameters_type]))
    #
    # Start calculations
    for j in range(projectile_n+1):
        # Open GSMCC input file
        with open(readfilename_CC,'r') as gsmin:
            inputfile_lines = gsmin.read().split('\n')
        # Find the basis.parameters line
        if projectile_type == 'proton' or projectile_type == 'neutron':
            theline = searchline_all(readfilename_CC,"Basis.WS.parameters")[0] # The nucleons are at the beginning of the input file
            shift = [x.strip(' ') for x in inputfile_lines[theline:theline+20]].index(projectile_type) + 2
        elif projectile_type in ['deuteron', 'triton', '3He', 'alpha']:
            theline = searchline(readfilename_CC,"cluster(s)")
            thelineend = searchline(readfilename_CC,"cluster.type")
            aux = indexof([x.strip(' ') for x in inputfile_lines[theline:thelineend]],'Basis.WS.parameters')[0]
            shift = [x.strip(' ') for x in inputfile_lines[theline+aux:thelineend]].index(projectile_type) + 2 + aux
        i = 0
        k = 0
        while i == 0:
            aux = inputfile_lines[theline + shift + k].split()
            if projectile_partialwaves_type == 'ALL':
                if j == 0:
                    aux_vo = float(aux[parameters_type + 1])
                    print_twice('INFO: Starting %s from %s for ALL l-waves'% (parameters_list[parameters_type], aux_vo))
                else:
                    aux_vo = float(aux[parameters_type + 1]) + projectile_step
                for ii in range(4):
                    if parameters_type == ii:
                        aux[ii] = str(aux_vo)
                        inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+aux[3]+'  '+aux[4]
            elif int(aux[0]) in projectile_lwave:
                if j == 0:
                    aux_vo = float(aux[parameters_type + 1])
                    print_twice('INFO: Starting %s from %s for l=%s'% (parameters_list[parameters_type], aux_vo, aux[0]))
                else:
                    aux_vo = float(aux[parameters_type + 1]) + projectile_step
                for ii in range(4):
                    if parameters_type == ii:
                        aux[i] = str(aux_vo)
                        inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+aux[3]+'  '+aux[4]
            k += 1
            if inputfile_lines[theline + shift + k].split() == []:
                i = 1
        # Save and close GSMCC input file
        inputfile_aux = '\n'.join(inputfile_lines)
        with open(readfilename_CC,'w') as gsmin:
            gsmin.write(inputfile_aux)
        #
        # Runnning the code
        start_gsmcc = time.time()
        outfilename_CC_j = logname + '/' + outfilename_CC[:-4] + '_%s.out'% (j+1)
        print_twice('\n ' + running_prefix + running_cc + ' < ' + readfilename_CC + cc_write+outfilename_CC_j)
        sp.run([running_prefix + running_cc + ' < ' + readfilename_CC + cc_write+outfilename_CC_j], shell=True)
        end_gsmcc = time.time()
        time_gsmcc = end_gsmcc-start_gsmcc
        print_twice("Time to calculate: ",time_gsmcc, "s")


#
end_main = time.time()
time_main = end_main-start_main
print_twice("\n\nAll calculations lasted: ", time_main, "s")

"""
    Example of input.BasisTest file:
    _________________________________
    # Remove the lines that you are not using, this example is with the complete test
    PROJECTILE_TYPE: name of the projectile as in GSMCC input file= "proton", "neutron", "deuteron", "triton", "3He", "alpha"
    alpha
    PARAMETER: options: 0 - diffuseness (a0), 1 - radius (R0), 2 - WS depth (VO), 3 - spin-orbit depth (VSO)
    2
    PARTIAL_WAVES: can be "ALL" for changing all the l-wave at the same time or "[0,2,3]" being 0,2,3 the partial waves to test
    ALL
    STRENGTHWS: if ALL, line 1 - step, line 2 - n points; if not, 1-2 for each partial wave defined before
    0.5
    10
    
    # Save WF?
    SAVE_WF: YES or NO
    YES
    WF_JPI_INDEX: If YES, indicate 1-number of states, and for each, 2-JPi, 3-Index
    8
    3/2-
    0
    3/2-
    1
    3/2-
    2
    3/2-
    3
    7/2+
    0
    7/2+
    1
    7/2+
    2
    7/2+
    3
    _________________________________
"""