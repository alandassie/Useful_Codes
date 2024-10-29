"""
    Created on August 2024 by Alan D.K. for 11C_Project

    This code run GSM and GSMC. You have to define the input and output files.
"""

import subprocess as sp
import time
import os
from scipy.optimize import least_squares as least
from scipy.optimize import newton

# LOG FILE
logfile = os.getcwd() + '/GSM+GSMCC_run.log'

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
def f(x):
    # Print input array
    print_twice('\n Correction factors:')
    print_twice(x)
    # Open GSMCC input file
    with open(readfilename_CC,'r') as gsmin:
        inputfile_lines = gsmin.read().split('\n')
    aux = inputfile_lines[icf_line].split()
    # inputfile_lines[icf_line] = "  " + aux[0] + " " + str(x[0])
    inputfile_lines[icf_line] = "  " + aux[0] + " " + str(x)
    # Save and close GSMCC input file
    inputfile_aux = '\n'.join(inputfile_lines)
    with open(readfilename_CC,'w') as gsmin:
        gsmin.write(inputfile_aux)
    #
    outfilename_CCi = outfilename_CC + str(int(time.time()))
    start_gsmcc = time.time()
    print_twice('\n mpirun -np 2 ./CC_exe < '+readfilename_CC+' > '+outfilename_CCi)
    sp.run(['mpirun -np 2 ./CC_exe < '+readfilename_CC+' > '+outfilename_CCi], shell=True)
    end_gsmcc = time.time()
    time_gsmcc = end_gsmcc-start_gsmcc
    print_twice("Time to calculate: ",time_gsmcc, "s")
    #
    # Read results
    with open(outfilename_CCi,'r') as gsmout:
        outputfile_lines = gsmout.read().split('\n')
    line_numbers = searchline_all(outfilename_CCi,search)
    energy = []
    width = []
    for line in line_numbers:
        auxe = outputfile_lines[line]
        ene = float(auxe.split(':')[1].split(' ')[1])
        auxw = outputfile_lines[line+1]
        wid = float(auxw.split(':')[1].split(' ')[1])
        energy.append( ene )
        width.append( wid )
    #
    energy.sort()
    print_twice('Calculated energies:')
    print_twice(energy)
    res = expene - float(energy[numberindex])
    print_twice('Residue = {0:7.3f}'.format(res))
    print_twice('\n'+20*'-'+'\n')
    #
    return res
# .-

# INPUT FILE
readfilename = "GSM+GSMCC_run_3I2-.in"
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
# RUN GSM
if theline != None:
    print_twice("\nRunning GSM in %s"% gsm_directory)
    os.chdir(gsm_directory)
    #
    for i in range(0,gsm_files):
        inp = readfilename_GSM[i]
        out = outfilename_GSM[i]
        start_gsm = time.time()
        print_twice('\n mpirun -np 2 ./GSM_exe < '+inp+' > '+out)
        sp.run(['mpirun -np 2 ./GSM_exe < '+inp+' > '+out], shell=True)
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
sp.run(['python3 EditThresholds_alt.py'], shell=True)
#
print_twice("\nRunning GSMCC in %s"% gsmcc_directory)
os.chdir(gsmcc_directory)
# Experimental state to optimize with interaction.corrective.factor
expene = -36.491
numberindex = 0
search = 'E(reference frame) :'
# Line with the corrective factor
icf_line = searchline(readfilename_CC,"CC.interaction.corrective.factor.composite(s)")+2
#
# opt = least(f,[1],diff_step=[0.02],gtol=1e-3,max_nfev=30, bounds=(0.8,1.2))
opt = newton(f,1,tol=1e-4,maxiter=30, full_output=True)
print_twice(opt)
#
end_main = time.time()
time_main = end_main-start_main
print_twice("\n\nAll calculations lasted: ", time_main, "s")
