"""
    Created on November 2024 by Alan D.K. for 11C_Project,
    in particular to calculate the alpha elastic dif. cross
    section with different pair (ene,ang) each time, in order
    to compare with the experiment H. Yamaguchi et al., Phys. Rev. C 87, 034303 (2013).

    This code run GSMCC multiple times in order to calculate
    the Differential Cross Section with different combinations
    of energies and angles. The data is read from a file called
    "GSMCC_DCS_ene-ang.in". An example can be found at the end of the code.
    
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
logfile = os.getcwd() + '/GSMCC_DCS_ene-ang.log'

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
readexpername = "GSMCC_DCS_ene-ang.in"
exp_data = np.genfromtxt(readexpername,skip_header=3,skip_footer=1)
energy = exp_data[:,0]
angle = exp_data[:,2]
with open(readexpername, 'r') as readfile:
    aux = readfile.read().split('\n')[0]
ncol = int(aux.split()[1])
nene = int(aux.split()[2])
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
angle_line = energy_line+3
for i in range(0,nene):
    energy_i = energy[i]
    angle_i = angle[i]
    # Open GSMCC input file
    with open(readfilename_CC,'r') as gsmin:
        inputfile_lines = gsmin.read().split('\n')
    inputfile_lines[energy_line] = str(energy_i) + " MeV(E.kinetic.total.system.CM)"
    inputfile_lines[angle_line] = "    " + str(angle_i)
    # Save and close GSMCC input file
    inputfile_aux = '\n'.join(inputfile_lines)
    with open(readfilename_CC,'w') as gsmin:
        gsmin.write(inputfile_aux)
    #
    start_gsmcc = time.time()
    print_twice('\n ' + running_prefix + running_cc + ' < ' + readfilename_CC + cc_write+outfilename_CC)
    sp.run([running_prefix + running_cc + ' < ' + readfilename_CC + cc_write+outfilename_CC], shell=True)
    end_gsmcc = time.time()
    time_gsmcc = end_gsmcc-start_gsmcc
    print_twice("Time to calculate: ",time_gsmcc, "s")
    # Save files
    folder = logfolder + '/' + str(round(energy_i,3)) + '+' + str(round(angle_i,3))
    if not os.path.exists(folder):
        os.makedirs(folder)
    filesincwd = os.listdir()
    for f in filesincwd:
        if (f.startswith('analyzing_power') or f.startswith('scattering_differential_cross_section') ):
            os.rename(f,folder + '/' + f)
#
end_main = time.time()
time_main = end_main-start_main
print_twice("\n\nAll calculations lasted: ", time_main, "s")

"""
    Example of GSMCC_DCS_ene-ang.in file:
    _________________________________
    DATA                 5        135
    EN-CM      EN-ERR    ANG-CM  DATA-CM    DATA-ERR
    MEV        MEV       ADEG    MB/SR      MB/SR
    1.23       0.10      142.3   7.89E+01   5.33E+00
    1.25       0.10      143.4   8.01E+01   5.48E+00
    1.28       0.10      145.0   7.08E+01   5.69E+00
    1.32       0.10      146.6   7.69E+01   5.74E+00
    1.35       0.09      148.2   7.86E+01   6.08E+00
    1.38       0.09      149.8   8.69E+01   6.23E+00
    1.41       0.09      151.4   7.74E+01   6.75E+00
    1.45       0.09      153.0   7.52E+01   6.62E+00
    1.48       0.09      154.6   6.63E+01   6.66E+00
    1.51       0.09      155.7   6.70E+01   6.79E+00
    1.54       0.09      156.1   6.26E+01   6.61E+00
    1.58       0.09      156.5   6.48E+01   6.70E+00
    1.61       0.09      156.9   6.85E+01   7.37E+00
    1.64       0.09      157.4   7.64E+01   7.01E+00
    1.67       0.09      157.8   6.71E+01   7.46E+00
    1.71       0.08      158.2   6.15E+01   7.39E+00
    1.74       0.08      158.6   6.39E+01   7.51E+00
    1.77       0.08      159.0   6.15E+01   7.34E+00
    1.80       0.08      159.4   6.08E+01   8.35E+00
    1.84       0.08      159.8   5.80E+01   7.20E+00
    1.87       0.08      160.2   6.39E+01   7.46E+00
    1.90       0.08      160.6   5.33E+01   8.67E+00
    1.93       0.08      161.0   5.59E+01   8.01E+00
    1.97       0.08      161.4   5.86E+01   8.81E+00
    2.00       0.08      161.8   5.64E+01   9.66E+00
    2.03       0.08      162.0   6.67E+01   8.15E+00
    2.06       0.08      162.1   6.08E+01   9.09E+00
    2.10       0.08      162.3   6.02E+01   9.02E+00
    2.13       0.08      162.5   7.74E+01   9.22E+00
    2.16       0.07      162.7   7.28E+01   9.86E+00
    2.19       0.07      162.8   8.45E+01   1.12E+01
    2.23       0.07      163.0   8.95E+01   1.13E+01
    2.26       0.07      163.2   1.02E+02   1.22E+01
    2.29       0.07      163.4   1.03E+02   1.29E+01
    2.32       0.07      163.5   1.30E+02   1.20E+01
    2.36       0.07      163.7   1.15E+02   1.23E+01
    2.39       0.07      163.9   1.21E+02   1.22E+01
    2.42       0.07      164.1   1.09E+02   1.28E+01
    2.45       0.07      164.2   7.70E+01   9.49E+00
    2.49       0.07      164.4   7.29E+01   9.72E+00
    2.52       0.07      164.6   5.31E+01   8.11E+00
    2.55       0.07      164.8   4.35E+01   8.04E+00
    2.58       0.07      165.0   2.27E+01   8.04E+00
    2.62       0.07      165.1   9.05E+00   5.24E+00
    2.65       0.07      165.3   7.39E+00   4.50E+00
    2.68       0.07      165.5   1.58E+01   5.14E+00
    2.71       0.07      165.7   1.32E+01   8.12E+00
    2.75       0.07      165.8   2.22E+01   8.04E+00
    2.78       0.06      166.0   2.90E+01   9.10E+00
    2.81       0.06      166.2   3.42E+01   8.64E+00
    2.84       0.06      166.4   2.43E+01   8.61E+00
    2.88       0.06      166.5   2.85E+01   7.89E+00
    2.91       0.06      166.7   2.68E+01   8.79E+00
    2.94       0.06      166.9   3.10E+01   1.04E+01
    2.97       0.06      167.1   2.72E+01   1.08E+01
    3.01       0.06      167.2   2.65E+01   1.24E+01
    3.04       0.06      167.3   2.98E+01   1.56E+01
    3.07       0.06      167.5   4.66E+01   1.62E+01
    3.10       0.06      167.6   4.53E+01   1.79E+01
    3.14       0.06      167.7   6.83E+01   1.76E+01
    3.17       0.06      167.8   9.91E+01   1.75E+01
    3.20       0.06      167.9   7.73E+01   1.98E+01
    3.23       0.06      168.0   1.09E+02   1.81E+01
    3.27       0.06      168.2   9.86E+01   2.05E+01
    3.30       0.06      168.3   8.74E+01   2.20E+01
    3.33       0.06      168.4   1.03E+02   2.10E+01
    3.36       0.06      168.5   8.12E+01   2.18E+01
    3.40       0.06      168.6   1.05E+02   2.33E+01
    3.43       0.06      168.7   1.13E+02   2.20E+01
    3.46       0.06      168.9   8.32E+01   2.29E+01
    3.49       0.06      169.0   6.87E+01   2.21E+01
    3.53       0.06      169.1   5.15E+01   2.10E+01
    3.56       0.06      169.2   3.75E+01   2.27E+01
    3.59       0.06      169.3   3.83E+01   1.94E+01
    3.62       0.06      169.4   3.85E+01   1.77E+01
    3.66       0.06      169.6   3.53E+01   1.59E+01
    3.69       0.06      169.7   2.79E+01   2.09E+01
    3.72       0.06      169.8   6.24E+01   2.86E+01
    3.75       0.06      169.9   3.01E+01   4.08E+01
    3.79       0.06      170.0   5.74E+01   3.84E+01
    3.82       0.06      170.1   4.34E+01   2.75E+01
    3.85       0.06      170.3   2.57E+01   1.77E+01
    3.88       0.05      170.4   2.81E+01   7.12E+00
    3.92       0.05      170.5   1.68E+01   9.95E+00
    3.95       0.05      170.6   1.54E+01   7.08E+00
    3.98       0.05      170.7   2.88E+01   1.04E+01
    4.01       0.05      170.8   3.48E+01   8.81E+00
    4.05       0.05      170.9   2.46E+01   9.68E+00
    4.08       0.05      171.0   2.84E+01   1.02E+01
    4.11       0.05      171.1   1.84E+01   8.53E+00
    4.14       0.05      171.1   3.68E+01   9.96E+00
    4.18       0.05      171.2   2.82E+01   9.13E+00
    4.21       0.05      171.3   2.26E+01   1.07E+01
    4.24       0.05      171.4   2.34E+01   7.96E+00
    4.27       0.05      171.5   2.85E+01   8.68E+00
    4.31       0.05      171.5   2.80E+01   8.68E+00
    4.34       0.05      171.6   2.44E+01   8.31E+00
    4.37       0.05      171.7   2.26E+01   9.04E+00
    4.40       0.05      171.8   2.31E+01   8.29E+00
    4.44       0.05      171.8   2.90E+01   1.16E+01
    4.47       0.05      171.9   2.81E+01   8.11E+00
    4.50       0.05      172.0   3.23E+01   1.06E+01
    4.53       0.05      172.1   6.92E+01   1.42E+01
    4.57       0.05      172.2   6.54E+01   1.49E+01
    4.60       0.05      172.2   8.81E+01   1.78E+01
    4.63       0.05      172.3   1.30E+02   1.91E+01
    4.66       0.05      172.4   1.43E+02   2.17E+01
    4.70       0.05      172.5   1.63E+02   2.20E+01
    4.73       0.05      172.5   2.07E+02   2.42E+01
    4.76       0.05      172.6   1.83E+02   2.21E+01
    4.79       0.05      172.7   2.36E+02   2.60E+01
    4.83       0.05      172.8   2.54E+02   2.77E+01
    4.86       0.05      172.9   2.34E+02   2.77E+01
    4.89       0.05      172.9   2.44E+02   2.86E+01
    4.92       0.05      173.0   2.44E+02   2.77E+01
    4.96       0.05      173.1   2.41E+02   3.12E+01
    4.99       0.05      173.2   2.44E+02   3.00E+01
    5.02       0.05      173.2   2.82E+02   2.93E+01
    5.05       0.05      173.3   2.98E+02   3.25E+01
    5.09       0.05      173.3   3.13E+02   3.19E+01
    5.12       0.05      173.4   2.92E+02   3.27E+01
    5.15       0.05      173.4   3.32E+02   3.38E+01
    5.18       0.05      173.4   3.30E+02   3.22E+01
    5.22       0.05      173.5   3.17E+02   3.31E+01
    5.25       0.05      173.5   3.02E+02   3.18E+01
    5.28       0.05      173.6   2.91E+02   3.15E+01
    5.31       0.05      173.6   2.54E+02   2.97E+01
    5.35       0.05      173.7   2.72E+02   3.09E+01
    5.38       0.05      173.7   2.07E+02   2.65E+01
    5.41       0.05      173.7   2.20E+02   2.75E+01
    5.44       0.05      173.8   1.94E+02   2.59E+01
    5.48       0.05      173.8   1.76E+02   2.48E+01
    5.51       0.05      173.9   1.62E+02   2.40E+01
    5.54       0.05      173.9   1.20E+02   2.29E+01
    5.57       0.05      174.0   1.34E+02   2.21E+01
    ENDDATA            137
    _________________________________
"""