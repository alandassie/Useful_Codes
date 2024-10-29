"""
    Created on June 2024 by Alan D.K.

    This code Test different shapes of the CC contours
"""

import numpy as np
import subprocess as sp


# Output file name
# now = datetime.now()
# optout = 'conttest%s.out'% (now.strftime("%y.%m.%d-%H:%M"))
# # Erease output file
# with open(optout,'w') as output:
#     pass
# Read file name
readfilename_array = ['IN2P3_11C_CC_GSMOpt-24.07.24-10.00_Basis-24.08.07-19.00_3I2-.in']

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


