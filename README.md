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

### GSMCC_DCS_ene-ang.py

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

### EdithThresholds.py

Created on July 2024 by Alan D.K. for 11C_Project.

This code will edit each one of the files created by
GSM of the form "eigenvector_E_averaged_n_scat_Za_Nb_JPi_i"
where a, b are the number of protons and neutrons, respectively,
JPi is the state of the Target and i is the index.

You have to define each tuple (a,b,JPi,i,ExpEner) with the
experimental energy measured respect to the core in a file
called "EdithThresholds.in" that must be placed in the workspace.
An example can be found at the end of the code.

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

### DATA_ANALYSIS/DataAnalysis_complete.py

Created on June 2024 by Alan D.K.

This is a multitask code used for the analysis of the GSMCC outputs.
The analysis done by the code are the followings:
+ sum up the partial widths for each mass partition.
+ sum up the occ prob for each mass partition 
and shown those greater than 5% of the total.
+ shows the calculated spectroscopic factors in Table format

The name of the file to be analyzed must be defined in the 
file "Data_Analysis.in". An example of the input file is at the end of the code.

