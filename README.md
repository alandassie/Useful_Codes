# Useful Codes
This is a set of python codes used for the calculation, minimization and/or optimization of results comming from differents many-body shell model codes.
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

### GSM+GSMCC_run.py

Created on August 2024 by Alan D.K. for 11C_Project

This code run GSM and GSMC. The running information must be in a file called GSM+GSMCC_run.in located in the main folder of the GSMCC code (i.e. CC_dir/GSM+GSMCC_run.py).
An example of the input file can be founded at the end of the code file.

### GSM+GSMCC_ICF.py

Created on October 2024 by Alan D.K. for 11C_Project

This code run GSM (if needed) and GSMCC and then
optimize the complex (or not) interaction corrective factor
to reproduce the experimental data. 
You have to define the input and output files and put
the experimental values (energy and width).

An example of the input file "GMS+GSMCC_ICF.in" can be found
at the end of the code file.

### BasisTest.py

Created on August 2024 by Alan D.K. for 11C_Project.

This code run multiples GSM and GSMCC calculations in order to calculate the GSMCC spectrum with different Basis.parameters.

Comment on October 29 2024: Outdated code, can be improved with an independent input file called BasisTest.in when needed!!

### ChannelTest.py

Created on August 2024 by Alan D.K. for 11C_Project

This code increase the number of channels in the CC calculations.

Comment on October 29 2024: Outdated code, can be improved with an independent input file called ChannelTest.in when needed!!

### ContourMod.py

Created on June 2024 by Alan D.K.

This code Test different shapes of the CC contours

Comment on October 29 2024: Outdated code, can be improved with an independent input file called ContourMod.in when needed!!

### ThrTest.py

Created on August 2024 by Alan D.K. for 11C_Project

This code run GSM and GSMC. You have to define the input and output files.
This code is used to test the behavoir of the GSMCC eigenvalues respect to the position of the thresholds.

Comment on October 29 2024: Outdated code, can be improved with an independent input file called ThrTest.in when needed!!




