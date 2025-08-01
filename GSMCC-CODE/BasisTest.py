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


# LOG FILE
logfile = os.getcwd() + '/log.BasisTest.' + time.strftime( "%y.%m.%d-%H.%M", time.localtime() )
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

# Preparing proton WS
theline = searchline(calcfilename,"PROTONWS")
if theline != None:
    print_twice('Doing Woods-Saxon proton test')
    proton_type = data[theline+1]
    if proton_type == 'ALL':
        print_twice('All WS partial waves equal')
        theline = searchline(calcfilename,"PROTONSTRENGTHWS")
        protonws_lwave_n = 0
        # protonws_starting_point = [float(data[theline+1])] # Read from the input file
        protonws_step = float(data[theline+1])
        protonws_n = int(data[theline+2])
    else:
        print_twice('Selected specific partial waves')
        protonws_lwave = ast.literal_eval(proton_type)
        protonws_lwave_n = len(protonws_lwave)
        theline = searchline(calcfilename,"PROTONSTRENGTHWS")
        # protonws_starting_point = [float(data[theline+1])] # Read from the input file
        protonws_step = float(data[theline+1])
        protonws_n = int(data[theline+2])
# Preparing proton SO 
theline = searchline(calcfilename,"PROTONSO")
if theline != None:
    print_twice('Doing spin-orbit proton test')
    proton_type = data[theline+1]
    if proton_type == 'ALL':
        print_twice('All SO partial waves equal')
        theline = searchline(calcfilename,"PROTONSTRENGTHSO")
        protonso_lwave_n = 0
        # protonso_starting_point = [float(data[theline+1])] # Read from the input file
        protonso_step = float(data[theline+1])
        protonso_n = int(data[theline+2])
    else:
        print_twice('Selected specific partial waves')
        protonso_lwave = ast.literal_eval(proton_type)
        protonso_lwave_n = len(protonso_lwave)
        theline = searchline(calcfilename,"PROTONSTRENGTHSO")
        protonso_lwave_n = 0
        # protonso_starting_point = [float(data[theline+1])] # Read from the input file
        protonso_step = float(data[theline+1])
        protonso_n = int(data[theline+2])

# Preparing neutron WS
theline = searchline(calcfilename,"NEUTRONWS")
if theline != None:
    print_twice('Doing Woods-Saxon neutron test')
    neutron_type = data[theline+1]
    if neutron_type == 'ALL':
        print_twice('All WS partial waves equal')
        theline = searchline(calcfilename,"NEUTRONSTRENGTHWS")
        neutronws_lwave_n = 0
        # neutronws_starting_point = [float(data[theline+1])] # Read from the input file
        neutronws_step = float(data[theline+1])
        neutronws_n = int(data[theline+2])
    else:
        print_twice('Selected specific partial waves')
        neutronws_lwave = ast.literal_eval(neutron_type)
        neutronws_lwave_n = len(neutronws_lwave)
        theline = searchline(calcfilename,"NEUTRONSTRENGTHWS")
        # neutronws_starting_point = [float(data[theline+1])] # Read from the input file
        neutronws_step = float(data[theline+1])
        neutronws_n = int(data[theline+2])
# Preparing neutron SO 
theline = searchline(calcfilename,"NEUTRONSO")
if theline != None:
    print_twice('Doing spin-orbit neutron test')
    neutron_type = data[theline+1]
    if neutron_type == 'ALL':
        print_twice('All SO partial waves equal')
        theline = searchline(calcfilename,"NEUTRONSTRENGTHSO")
        neutronso_lwave_n = 0
        # neutronso_starting_point = [float(data[theline+1])] # Read from the input file
        neutronso_step = float(data[theline+1])
        neutronso_n = int(data[theline+2])
    else:
        print_twice('Selected specific partial waves')
        neutronso_lwave = ast.literal_eval(neutron_type)
        neutronso_lwave_n = len(neutronso_lwave)
        theline = searchline(calcfilename,"NEUTRONSTRENGTHSO")
        # neutronso_starting_point = [float(data[theline+1])] # Read from the input file
        neutronso_step = float(data[theline+1])
        neutronso_n = int(data[theline+2])
        
# Preparing both WS
theline = searchline(calcfilename,"BOTHWS")
if theline != None:
    print_twice('Doing Woods-Saxon test for proton and neutron at the same time')
    both_proton_type = data[theline+1]
    both_neutron_type = data[theline+2]
    if both_proton_type == 'ALL' and both_neutron_type == 'ALL':
        print_twice('All WS partial waves equal in neutron and proton')
        theline = searchline(calcfilename,"BOTH_STRENGTHWS")
        both_lwave_n = 0
        # protonws_starting_point = [float(data[theline+1])] # Read from the input file
        both_proton_step = float(data[theline+1])
        both_neutron_step = float(data[theline+2])
        both_n = int(data[theline+3])
    else:
        print_twice('Selected specific partial waves')
        both_protonws_lwave = ast.literal_eval(both_proton_type)
        both_protonws_lwave_n = len(both_protonws_lwave)
        both_neutronws_lwave = ast.literal_eval(both_neutron_type)
        both_neutronws_lwave_n = len(both_neutronws_lwave)
        theline = searchline(calcfilename,"BOTH_STRENGTHWS")
        both_proton_step = float(data[theline+1])
        both_neutron_step = float(data[theline+2])
        both_n = int(data[theline+3])
        
# Preparing alpha WS
theline = searchline(calcfilename,"ALPHAWS")
if theline != None:
    print_twice('Doing Woods-Saxon test for alpha')
    alpha_type = data[theline+1]
    if alpha_type == 'ALL':
        print_twice('All WS partial waves at the same time')
        theline = searchline(calcfilename,"ALPHA_STRENGTHWS")
        alpha_lwave_n = 0
        # protonws_starting_point = [float(data[theline+1])] # Read from the input file
        alpha_step = float(data[theline+1])
        alpha_n = int(data[theline+2])
    else:
        print_twice('Selected specific partial waves')
        alpha_lwave = ast.literal_eval(alpha_type)
        alpha_lwave_n = len(alpha_lwave)
        theline = searchline(calcfilename,"ALPHA_STRENGTHWS")
        alpha_step = float(data[theline+1])
        alpha_n = int(data[theline+2])

        
# Saving WF
theline = searchline(calcfilename,"SAVE_WF:")
if theline != None:
    savingwf_info = data[theline+1]
    if savingwf_info == 'NO':
        print_twice('The channels WF of each calculation will not be saved')
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
# Calculating neutron WS
theline = searchline(calcfilename,"NEUTRONWS")
if theline != None:
    print_twice('Start neutron WS calculations')
    neutron_type = data[theline+1]
    # Start calculations
    for j in range(neutronws_n+1):
        # Open GSMCC input file
        with open(readfilename_CC,'r') as gsmin:
            inputfile_lines = gsmin.read().split('\n')
        # Find the basis.parameters line
        theline = searchline_all(readfilename_CC,"Basis.WS.parameters")[0]
        shift = [x.strip(' ') for x in inputfile_lines[theline:theline+20]].index('neutron') + 2
        i = 0
        k = 0
        while i == 0:
            aux = inputfile_lines[theline + shift + k].split()
            if neutron_type == 'ALL':
                if j == 0:
                    aux_vo = float(aux[3])
                else:
                    aux_vo = float(aux[3]) + neutronws_step
                inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+str(aux_vo)+'  '+aux[4]
                if inputfile_lines[theline + shift + k].split() == []:
                    i = 1
            elif int(aux[0]) in neutronws_lwave:
                if j == 0:
                    aux_vo = float(aux[3])
                else:
                    aux_vo = float(aux[3]) + neutronws_step
                inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+str(aux_vo)+'  '+aux[4]
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
        outfilename_CC_j = logname + '/' + outfilename_CC[:-4] + '_%s.out'% j
        print_twice('\n ' + running_prefix + running_cc + ' < ' + readfilename_CC + cc_write+outfilename_CC_j)
        sp.run([running_prefix + running_cc + ' < ' + readfilename_CC + cc_write+outfilename_CC_j], shell=True)
        end_gsmcc = time.time()
        time_gsmcc = end_gsmcc-start_gsmcc
        print_twice("Time to calculate: ",time_gsmcc, "s")
#
# Calculating proton WS
theline = searchline(calcfilename,"PROTONWS")
if theline != None:
    print_twice('Start proton WS calculations')
    proton_type = data[theline+1]
    # Start calculations
    for j in range(protonws_n+1):
        # Open GSMCC input file
        with open(readfilename_CC,'r') as gsmin:
            inputfile_lines = gsmin.read().split('\n')
        # Find the basis.parameters line
        theline = searchline_all(readfilename_CC,"Basis.WS.parameters")[0]
        shift = [x.strip(' ') for x in inputfile_lines[theline:theline+20]].index('proton') + 2
        i = 0
        k = 0
        while i == 0:
            aux = inputfile_lines[theline + shift + k].split()
            if proton_type == 'ALL':
                if j == 0:
                    aux_vo = float(aux[3])
                else:
                    aux_vo = float(aux[3]) + protonws_step
                inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+str(aux_vo)+'  '+aux[4]
            elif int(aux[0]) in protonws_lwave:
                if j == 0:
                    aux_vo = float(aux[3])
                else:
                    aux_vo = float(aux[3]) + protonws_step
                inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+str(aux_vo)+'  '+aux[4]
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
# Calculating both WS
theline = searchline(calcfilename,"BOTHWS")
if theline != None:
    print_twice('Start WS test')
    both_proton_type = data[theline+1]
    both_neutron_type = data[theline+2]
    # Start calculations
    for j in range(both_n+1):
        # Open GSMCC input file
        with open(readfilename_CC,'r') as gsmin:
            inputfile_lines = gsmin.read().split('\n')
        # Find the basis.parameters line
        theline = searchline_all(readfilename_CC,"Basis.WS.parameters")[0]
        # First edit proton
        shift = [x.strip(' ') for x in inputfile_lines[theline:theline+20]].index('proton') + 2
        i = 0
        k = 0
        while i == 0:
            aux = inputfile_lines[theline + shift + k].split()
            if both_proton_type == 'ALL':
                if j == 0:
                    aux_vo = float(aux[3])
                else:
                    aux_vo = float(aux[3]) + both_proton_step
                inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+str(aux_vo)+'  '+aux[4]
            elif int(aux[0]) in both_protonws_lwave:
                if j == 0:
                    aux_vo = float(aux[3])
                else:
                    aux_vo = float(aux[3]) + both_proton_step
                inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+str(aux_vo)+'  '+aux[4]
            k += 1
            if inputfile_lines[theline + shift + k].split() == []:
                i = 1
        # Then edit neutron
        shift = [x.strip(' ') for x in inputfile_lines[theline:theline+20]].index('neutron') + 2
        i = 0
        k = 0
        while i == 0:
            aux = inputfile_lines[theline + shift + k].split()
            if both_neutron_type == 'ALL':
                if j == 0:
                    aux_vo = float(aux[3])
                else:
                    aux_vo = float(aux[3]) + both_neutron_step
                inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+str(aux_vo)+'  '+aux[4]
            elif int(aux[0]) in both_neutronws_lwave:
                if j == 0:
                    aux_vo = float(aux[3])
                else:
                    aux_vo = float(aux[3]) + both_neutron_step
                inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+str(aux_vo)+'  '+aux[4]
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
# Calculating alpha WS
theline = searchline(calcfilename,"ALPHAWS")
if theline != None:
    print_twice('Start alpha WS calculations')
    alpha_type = data[theline+1]
    # Start calculations
    for j in range(alpha_n+1):
        # Open GSMCC input file
        with open(readfilename_CC,'r') as gsmin:
            inputfile_lines = gsmin.read().split('\n')
        # Find the basis.parameters line
        theline = searchline_all(readfilename_CC,"Basis.WS.parameters")[-1]
        thelineend = searchline(readfilename_CC,"cluster.type")
        shift = [x.strip(' ') for x in inputfile_lines[theline:thelineend]].index('alpha') + 2
        i = 0
        k = 0
        while i == 0:
            aux = inputfile_lines[theline + shift + k].split()
            if alpha_type == 'ALL':
                if j == 0:
                    aux_vo = float(aux[3])
                else:
                    aux_vo = float(aux[3]) + alpha_step
                inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+str(aux_vo)+'  '+aux[4]
            elif int(aux[0]) in alpha_lwave:
                if j == 0:
                    aux_vo = float(aux[3])
                else:
                    aux_vo = float(aux[3]) + alpha_step
                inputfile_lines[theline + shift + k] = '    '+aux[0]+'   '+aux[1]+'   '+aux[2]+'    '+str(aux_vo)+'  '+aux[4]
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
# Restore GSMCC file
startingpoint_aux = '\n'.join(startingpoint_lines)
with open(readfilename_CC,'w') as gsmin:
    gsmin.write(startingpoint_aux)


#
end_main = time.time()
time_main = end_main-start_main
print_twice("\n\nAll calculations lasted: ", time_main, "s")

"""
    Example of input.BasisTest file:
    _________________________________
    # Remove the lines that you are not using, this example is with the complete test
    PROTONWS: can be "ALL" for changing all the l-wave at the same time or "[0,2,3]" being 0,2,3 the partial waves to test
    ALL
    PROTONSTRENGTHWS: if ALL, 1 - step, 2 - n points; if not, 1-2 for each partial wave defined before
    0.2
    10
    PROTONSO: can be "ALL" for changing all the l-wave at the same time or "[0,2,3]" being 0,2,3 the partial waves to test
    [0,2]
    PROTONSTRENGTHSO: if ALL, 1 - step, 2 - n points; if not, 1-2 for each partial wave defined before
    0.1
    3
    0.05
    2
    
    NEUTRONWS: can be "ALL" for changing all the l-wave at the same time or "[0,2,3]" being 0,2,3 the partial waves to test
    ALL
    NEUTRONSTRENGTHWS: if ALL, 1 - step, 2 - n points; if not, 1-2 for each partial wave defined before
    54
    0.2
    10
    NEUTRONSO: can be "ALL" for changing all the l-wave at the same time or "[0,2,3]" being 0,2,3 the partial waves to test
    ALL
    NEUTRONSTRENGTHSO: if ALL, 1 - step, 2 - n points; if not, 1-2 for each partial wave defined before
    0.1
    3
    
    BOTHWS: can be "ALL" for changing all the l-wave at the same time or "[0,2,3]" being 0,2,3 the partial waves to test
    [0]
    [0]
    BOTH_STRENGTHWS: if ALL, 1 - proton step, 2 - neutron step, 3 - n points; if not, 1-2-3 for each partial wave defined before
    1
    1
    20
    
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