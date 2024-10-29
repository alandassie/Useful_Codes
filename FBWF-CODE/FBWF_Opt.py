"""
    Created on 2022 by Alan D.K.

    This code optimize the four body energy over a core,
    by adjusting the potential strengths of the np interaction
    to get the experimental value. For this matter, it use the
    least_squares minimization package from the scipy.optimize library.

    It can be set to optimize only the T=0 strengths, or the full Vnp. 

    The necessary data are the initial values, the experimental energy,
    the angular momentum J and parity, and the initial chi value, 
    all read from the file "Read/alwf*.in".

    You have to change the lines 29 and 30 to set the corresponding
    filename and spaces for the running version.
"""

import numpy as np
import os
import subprocess as sp
from scipy.optimize import least_squares as least

# Output file name
optout = 'Out/opt.out'
# Erease output file
with open(optout,'w') as output:
    pass

# Input file, depends on the running version
namefile = 'alwf_Dif.J.in' # The options are 'alwf.in' and 'alwf_Dif.J.in'
ls = "          " # The options are "" or "          "

# Declaration of funcitons
def print_twice(*args,**kwargs): # Allow to print in file and in terminal at the same line
    print(*args,**kwargs)
    with open(optout,"a") as f:  # appends to file and closes it when finished
        print(file=f,*args,**kwargs)
# .-
def f(x):

    print_twice('\nNew \u03C7:')
    print_twice('\n    \u03C7 = {0:9.4f}'.format(x[0]))

    if optmode == 1:
        vse = vse_ini*x
        vto = vto_ini*x
        print_twice('\nNew T=1 strengths:')
        print_twice('    V_se = '+' '.join('{0:9.4f}'.format(val) for val in vse) )
        print_twice('    V_to = '+' '.join('{0:9.4f}'.format(val) for val in vto) )
    vso = vso_ini*x
    vte = vte_ini*x

    print_twice('\nNew T=0 strengths:')
    print_twice('    V_so = {0:9.4f}'.format(vso[0]))
    print_twice('    V_te = {0:9.4f}'.format(vte[0]))

    # Changing the Vso and Vte Potentials
    with open(readfilename,'r') as alwfin:
        alwflines = alwfin.read().split('\n')
    #
    if optmode == 1:
        alwflines[11] = ' '.join('{0:9.4f}'.format(val) for val in vse) + ' {0:5.3f}  !V_se beta_se'.format(betase)
        alwflines[12] = ' '.join('{0:9.4f}'.format(val) for val in vto) + ' {0:5.3f}  !V_to beta_to'.format(betato)
    alwflines[13] = ls+'{0:9.4f} {1:5.3f}  !V_so beta_so'.format(vso[0],betaso)
    alwflines[14] = ls+'{0:9.4f} {1:5.3f}  !V_te beta_te'.format(vte[0],betate)
    #
    alwf = '\n'.join(alwflines)
    with open(readfilename,'w') as alwfout:
        alwfout.write(alwf)

    #Executing the code
    sp.run(['./alwf'])

    #________________________________________________
    # Read the ground state result
    gs = np.min(np.genfromtxt('Out/eigen.out'))
    # Compare with experimental energy
    res = gs-exp
    print_twice('\nCalculated ground state energy = {0:7.3f}'.format(gs))
    print_twice('Residue = {0:7.3f}'.format(res))
    print_twice('\n'+20*'-'+'\n')

    return res
# .-

# Reading part
# .-
# From the file
readfilename = 'Read/'+namefile 
# read the strengths initial values, the experimental energy, 
# the angular momentum J and the initial chi value
with open(readfilename, 'r') as readfile:
    data = readfile.read().split('\n')
    #
# SINGLET EVEN
array   = np.fromstring(data[11].split('!')[0], dtype=float, sep=' ')
lim     = len(array)-1
vse_ini = array[:len(array)-1] # TB optimized singlet even Strength(s)
betase  = array[len(array)-1]  # TB optimized singlet even Beta
# TRIPLET ODD
array   = np.fromstring(data[12].split('!')[0], dtype=float, sep=' ')
lim     = len(array)-1
vto_ini = array[:len(array)-1] # TB optimized triplet odd Strength
betato  = array[len(array)-1]  # TB optimized triplet odd Beta
# SINGLET ODD
array   = np.fromstring(data[13].split('!')[0], dtype=float, sep=' ')
lim     = len(array)-1
vso_ini = array[:len(array)-1] # TB optimized singlet odd Strength(s)
betaso  = array[len(array)-1]  # TB optimized singlet odd Beta
# TRIPLET EVEN
array   = np.fromstring(data[14].split('!')[0], dtype=float, sep=' ')
lim     = len(array)-1
vte_ini = array[:len(array)-1] # TB optimized triplet even Strength
betate  = array[len(array)-1]  # TB optimized triplet even Beta
# 
J  = np.fromstring(data[24].split('!')[0], dtype=int, sep=',')[0] # FB angular momentum 
Pi = np.fromstring(data[24].split('!')[0], dtype=int, sep=',')[1] # FB parity
if Pi == 1 :
    parity = '+'
else:
    parity = '-'
#
exp = float(data[28].split('=')[1]) # Experimental energy
#
seed    = float(data[27].split('=')[1]) # \chi initial value
optmode = int(data[26].split('=')[1].split('!')[0]) # Optimization mode
# .-

print_twice('---Start optimization FBWF code---\n')
print_twice(30*'-'+'\n')

print_twice('---Read data from file Read/opt.in---')
#
print_twice('    Opt four body state J^Pi = {0:2d}^{1:1s}'.format(J, parity))
print_twice('    Experimental energy = {0:7.3f} MeV'.format(exp))
print_twice('    Potentials optimized for two body np levels:')
print_twice('        V_se = '+' '.join('{0:9.4f}'.format(val) for val in vse_ini) )
print_twice('        V_to = '+' '.join('{0:9.4f}'.format(val) for val in vto_ini) )
print_twice('        V_so = '+ls+'{0:9.4f}'.format(vso_ini[0]))
print_twice('        V_te = '+ls+'{0:9.4f}'.format(vte_ini[0]))
print_twice('    \u03C7 seed = {0:7.3f}'.format(seed))
if optmode == 0:
    print_twice('    Only Vnp T=0 will be optimize')
else:
    print_twice('    Full Vnp will be optimize')

print_twice('\n'+30*'-'+'\n')


# Start optimization
print_twice('---Start optimization part---')
opt = least(f,seed,gtol=5e-4,diff_step=[5e-5],max_nfev=10)

print_twice('\n Optimized \u03C7 = {0:9.4f}'.format(opt.x[0]))

if optmode == 1:
    chit1 = opt.x[0]
else: 
    chit1 = 1
print_twice('\nFinal potentials strengths:')
print_twice('    V_se = '+' '.join('{0:9.4f}'.format(val) for val in vse_ini*chit1) )
print_twice('    V_to = '+' '.join('{0:9.4f}'.format(val) for val in vto_ini*chit1) )
print_twice('    V_so = '+ls+'{0:9.4f}'.format(vso_ini[0]*opt.x[0]))
print_twice('    V_te = '+ls+'{0:9.4f}'.format(vte_ini[0]*opt.x[0]))

print_twice('\n'+30*'-'+'\n')

print_twice('\nValue of the cost function: {0:9.6e}'.format(opt.cost))
print_twice('\nVector of residuals: {0:9.6e}'.format(opt.fun[0]))
print_twice('\nGradient of the cost function: {0:9.6e}'.format(opt.grad[0]))
print_twice('\nThe reason for algorithm termination: {0:2d}'.format(opt.status))

# Restart the input file with the original values
with open(readfilename,'r') as alwfin:
    alwflines = alwfin.read().split('\n')
#
if optmode == 1:
    alwflines[11] = ' '.join('{0:9.4f}'.format(val) for val in vse_ini) + ' {0:5.3f}  !V_se beta_se'.format(betase)
    alwflines[12] = ' '.join('{0:9.4f}'.format(val) for val in vto_ini) + ' {0:5.3f}  !V_to beta_to'.format(betato)
alwflines[13] = ls+'{0:9.4f} {1:5.3f}  !V_so beta_so'.format(vso_ini[0],betaso)
alwflines[14] = ls+'{0:9.4f} {1:5.3f}  !V_te beta_te'.format(vte_ini[0],betate)
#
alwf = '\n'.join(alwflines)
with open(readfilename,'w') as alwfout:
    alwfout.write(alwf)


"""
    This program is free software: you can redistribute it and/or 
    modify it under the terms of the GNU General Public License as 
    published by the Free Software Foundation, either version 3 of 
    the License, or (at your option) any later version.
"""