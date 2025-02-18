"""
    Created on February 2024 by Alan D.K. for 11C_Project,
    in particular to calculate the radiative cross sections
    for an specific range of E_cm energies.

    Then, the code run GSMCC multiple times in order to calculate
    the radiative Cross Section with different energies. 
    The data is read from a file called "GSMCC_rCS.in". 
    An example can be found at the end of the code.
    
    The code also use the input file "GSM+GSMCC_run.in" of the code
    GSM+GSMCC_run.py as guide of calculations. An example of that file
    can be found at the end of the GSM+GSMCC_run.py code.
"""

import subprocess as sp
import time
import os
import math as m
import numpy as np

# LOG FILE
logfile = os.getcwd() + '/log.GSMCC_rCS'

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

# CALCULATION INPUT FILE
readfilename = "GSM+GSMCC_run.in"
with open(readfilename, 'r') as readfile:
    data = readfile.read().split('\n')
# DATA INPUT FILE
readexpername = "input.GSMCC_rCS"
theline = searchline(readexpername,"ENERGY START POINT (MeV):")
init_ene = float(data[theline+1])
theline = searchline(readexpername,"ENERGY END (MeV):")
end_ene = float(data[theline+1])
theline = searchline(readexpername,"NUMBER OF ENERGY POINTS:")
n_ene = int(data[theline+1])
energy = np.linspace(init_ene,end_ene,n_ene)
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
    running_prefix = ' ./'
# Executable GSMCC file
theline = searchline(readfilename,"GSMCC-exe:")
running_cc = data[theline+1]
# Checking if we need machinefile
theline = searchline(readfilename,"MACHINEFILE:")
if theline != None:  
    machinefile_name = data[theline+1]
    running_prefix = running_prefix + running_cc + '-hostfile ' + machinefile_name + ' '

# Read-Out file name CC
theline = searchline(readfilename,"GSMCC-files:")
readfilename_CC = data[theline+1]
outfilename_CC = data[theline+2]
cc_write_aux = int(data[theline+3])
cc_write = ' ' + cc_write_aux*'>' + ' '
#
logfolder = os.getcwd() + '/' + outfilename_CC[:-4]
if not os.path.exists(logfolder):
    os.makedirs(logfolder)

# Start calculation
start_main = time.time()
# Edit thresholds
print_twice("\nEdit thresholds in %s"% storage_directory)
os.chdir(storage_directory)
sp.run(['python3 Useful_Codes/GSMCC-CODE/EditThresholds.py'], shell=True)
#
# Run GSMCC for each pair ene-ang, and move the files to different folders
print_twice("\nRunning GSMCC in %s"% gsmcc_directory)
os.chdir(gsmcc_directory)
#
# Line with energy, angle
energy_line = searchline(readfilename_CC,"MeV(E.kinetic.total.system.CM)")
for i in range(0,n_ene):
    energy_i = energy[i]
    # Open GSMCC input file
    with open(readfilename_CC,'r') as gsmin:
        inputfile_lines = gsmin.read().split('\n')
    inputfile_lines[energy_line] = str(energy_i) + " MeV(E.kinetic.total.system.CM)"
    # Save and close GSMCC input file
    inputfile_aux = '\n'.join(inputfile_lines)
    with open(readfilename_CC,'w') as gsmin:
        gsmin.write(inputfile_aux)
    #
    outfilename_CCi = str(i) + '_E=' + str(round(energy_i,3))
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

"""
    Example of input.GSMCC_rCS file:
    _________________________________
    ENERGY START POINT (MeV):
    1
    ENERGY END POINT (MeV):
    4
    NUMBER OF ENERGY POINTS:
    30
    _________________________________
"""