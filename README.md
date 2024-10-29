# Optimization Codes
This is a set of python codes used for the minimization and/or optimization of calculated observables to experimental values.
It originaly started with the FBWF code used on my PhD Thesis, and then continue with the post-doc at GANIL.
We will present each code separatly by folders. Remember that the code will act on the previously folder, this is, in "../".

## FBWF-CODE folder

Here we present the codes related to the FBWF code used during the PhD Thesis.

### FBWF_Opt.py
This code optimize the four body energy over a core, by adjusting the potential strengths of the np interaction to get the experimental value. For this matter, it use the least_squares minimization package from the scipy.optimize library.

It can be set to optimize only the T=0 strengths, or the full Vnp.  The necessary data are the initial values, the experimental energy, the angular momentum J and parity, and the initial chi value, all read from the file "Read/alwf*.in".

You have to change the lines 29 and 30 to set the corresponding  filename and spaces for the running version.

## GSMCC-CODE folder

Here we present the codes related to the GSMCC code.

### BasisTest.py

Created on August 2024 by Alan D.K. for 11C_Project.

This code run multiples GSM and GSMCC calculations in order to calculate the GSMCC spectrum with different Basis.parameters.

Comment on October 29: Outdated code, can be improved with an independent input file called BasisTest.in when needed!!

### ChannelTest.py

Created on August 2024 by Alan D.K. for 11C_Project

This code increase the number of channels in the CC calculations.

Comment on October 29: Outdated code, can be improved with an independent input file called ChannelTest.in when needed!!
