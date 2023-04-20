# Minimization Codes
This is a set of python codes used for the minimization and/or optimization of calculated observables to experimental values.

## This is the full list of codes

### FBWF_Opt.py
This code optimize the four body energy over a core, by adjusting the potential strengths of the np interaction to get the experimental value. For this matter, it use the least_squares minimization package from the scipy.optimize library.

It can be set to optimize only the T=0 strengths, or the full Vnp.  The necessary data are the initial values, the experimental energy, the angular momentum J and parity, and the initial chi value, all read from the file "Read/alwf*.in".

You have to change the lines 29 and 30 to set the corresponding  filename and spaces for the running version.

