"""
    Created on June 2024 by Alan D.K.

    This code sum up the width from a GSMCC output in order to distinguish 
    the mass partition which is more important.
    
    The name of the input file must be defined in the file "Data_Analysis.in"
"""

from datetime import datetime

# Output file name
# now = datetime.now()
# optout = 'Sum%s.out'% (now.strftime("%y.%m.%d_%H:%M"))
# # Erease output file
# with open(optout,'w') as output:
#     pass

# Declaration of funcitons
def print_twice(*args,**kwargs): # Allow to print in file and in terminal at the same line
    print(*args,**kwargs)
    with open(optout,"a") as f:  # appends to file and closes it when finished
        print(file=f,*args,**kwargs)
        
# Input file
inputfile = 'Data_Analysis.in' 
with open(inputfile, 'r') as readfile:
    dataaux = readfile.read().split('\n')
namefile = dataaux[0]
# Reading from the file
with open(readfilename, 'r') as readfile:
    data = readfile.read().split('\n')
    #
# Build two arrays
projectile_width = []
current_width = []
partial_current_width = []
current_width_aux = 0
altcurrent_width = []
partial_altcurrent_width = []
altcurrent_width_aux = 0
j = 0
i = -1
projectile_old = " "
for x in data:
    if x == '':
        continue
    projectile = x.split(' : ')[2].split(' ')[0]
    aux1 = x.split(' : ')[3].split(' , ')
    aux2 = float(aux1[0].split(' keV ')[0])
    aux3 = float(aux1[1].split(' keV ')[0])
    if projectile == projectile_old :
        current_width_aux += aux2
        altcurrent_width_aux += aux3
        partial_current_width[i].append(aux2)
        partial_altcurrent_width[i].append(aux3)
    else:
        if j == 0:
            i += 1
            current_width_aux += aux2
            altcurrent_width_aux += aux3
            projectile_old = projectile
            partial_current_width.append([aux2])
            partial_altcurrent_width.append([aux3])
            projectile_width.append(projectile)
            j +=1
        else:
            current_width.append(current_width_aux)
            altcurrent_width.append(altcurrent_width_aux)
            #
            i+=1
            projectile_old = projectile
            projectile_width.append(projectile)
            current_width_aux = aux2
            altcurrent_width_aux = aux3
            partial_current_width.append([aux2])
            partial_altcurrent_width.append([aux3])
current_width.append(current_width_aux)
altcurrent_width.append(altcurrent_width_aux)

            
for i in range(0,len(projectile_width)):
    print('Projectile mass partition : {0:s}'.format(projectile_width[i]))
    print('Total current width (keV) : {0:.5e}'.format(current_width[i]))
    print('Max partial current width (keV) : {0:.5e}'.format(max(partial_current_width[i])))
    print('Total alt current width (keV) : {0:.5e}'.format(altcurrent_width[i]))
    print('Max partial alt current width (keV) : {0:.5e}'.format(max(partial_altcurrent_width[i])))
    print('\n')






"""
    This program is free software: you can redistribute it and/or 
    modify it under the terms of the GNU General Public License as 
    published by the Free Software Foundation, either version 3 of 
    the License, or (at your option) any later version.
"""
