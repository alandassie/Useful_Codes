"""
    Created on July 2024 by Alan D.K. for 11C_Project.

    This code will edit each one of the files created by
    GSM of the form "eigenvector_E_averaged_n_scat_Za_Nb_JPi_i"
    where a, b are the number of protons and neutrons, respectively,
    JPi is the state of the Target and i is the index.
    
    You have to define each tuple (a,b,JPi,i,ExpEner) with the
    experimental energy measured respect to the core in a file
    called "EdithThresholds.in" that must be placed in the workspace.
    An example can be found at the end of the code.
"""

# Declaration of funcitons
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
#.-

# import numpy as np
# print("\n...The desired target/projectile's states to modify MUST be defined in the python file...\n")
# defined = input("Are the states and experimental energies defined?\n YES or NO: ")
# if defined.lower() == "no":
#     print("Do it!")
#     exit()

# INPUT FILE
readfilename = "EditThresholds.in"
with open(readfilename, 'r') as readfile:
    data = readfile.read().split('\n')
data = data[:-1]
# Define the array of tuples for each state
# Define one array for each target and projectile, and then sum up in a final array
# .-
# TARGETS
theline = searchline(readfilename,"TARGETS")
targets = []
i = 0
while i < len(data)-3:
    if theline+2+i >= len(data):
        break
    aux = data[theline+2+i]
    if aux == "PROJECTILES":
        break
    targets.append( ( int(aux.split(",")[0]),
                    int(aux.split(",")[1]),
                    str(aux.split(",")[2]),
                    int(aux.split(",")[3]),
                    float(aux.split(",")[4]) ) )
    i += 1

# PROJECTILES
theline = searchline(readfilename,"PROJECTILES")
if theline == None:  
    projectiles = []
else:
    projectiles = []
    i = 0
    while i < len(data)-3:
        try:
            aux = data[theline+2+i]
        except:
            break
        projectiles.append( ( int(aux.split(",")[0]),
                            int(aux.split(",")[1]),
                            str(aux.split(",")[2]),
                            int(aux.split(",")[3]),
                            float(aux.split(",")[4]) ) )
        i += 1

keyname = "eigenvector_E_averaged_n_scat_"

# print("...BackUp the GSM calculated version...\n")
for state in targets + projectiles:
    readfile = keyname+"Z%s_N%s_%s_%s"% state[0:4]
    savefile = keyname+"Z%s_N%s_%s_%s_backup"% state[0:4]
    # Check if the state was calculated before
    try:
        with open(readfile,'r') as f:
            line = f.read().split('\n')
    except FileNotFoundError:
        continue
    #
    # If file exist, then append instead of overwriting
    with open(savefile,'a') as fl:
        fl.write(f'\n{line[0]}')
#
# print("...Change threshold files with experimental values...")
# print("\n...Change projectile binding energy...")
for proj in projectiles:
    readfile = keyname+"Z%s_N%s_%s_%s"% proj[0:4]
    try:
        with open(readfile,'r') as f:
            line = f.read().split('\n')
        # print("File: "+readfile)
    except FileNotFoundError:
        # print("The state %s_%s of the Projectile Z:%s,N:%s was not calculated with GSM!"% ( state[2],state[3],state[0],state[1] ))
        continue
    newene = "(%s,0) "% proj[4]
    line[0] = newene + line[0].split(" ")[1]
    with open(readfile,'w') as f:
        f.write(line[0])
# print("\n...Change target binding energy...\n")
for state in targets:
    readfile = keyname+"Z%s_N%s_%s_%s"% state[0:4]
    try:
        with open(readfile,'r') as f:
            line = f.read().split('\n')
        # print("File: "+readfile)
    except FileNotFoundError:
        # print("The state %s_%s of the Target Z:%s,N:%s was not calculated with GSM!"% ( state[2],state[3],state[0],state[1] ))
        continue
    newene = "(%s,0) "% state[4]
    line[0] = newene + line[0].split(" ")[1]
    with open(readfile,'w') as f:
        f.write(line[0])
# print("\n...Finish file editing...")


"""
    Example of EditThresholds.in file:
    _________________________________
    TARGETS
    Z,N,JPi,i,ENE
    4,3,3I2-,0,-9.305
    4,3,3I2-,1,0.595
    4,3,1I2-,0,-8.876
    4,3,1I2-,1,7.695
    4,3,7I2-,0,-4.735
    4,3,7I2-,1,-0.035
    4,3,5I2-,0,-2.575
    4,3,5I2-,1,-2.095
    4,4,0+,0,-28.204
    4,4,2+,0,-25.174
    4,4,4+,0,-16.854
    5,5,3+,0,-36.455
    5,5,3+,1,-31.681
    5,5,3+,2,-29.451
    5,5,1+,0,-35.737
    5,5,1+,1,-34.301
    5,5,1+,2,-31.273
    5,5,1+,3,-28.789
    5,5,0+,0,-34.715
    5,5,0+,1,-28.895
    5,5,2+,0,-32.868
    5,5,2+,1,-31.291
    5,5,2+,2,-30.536
    5,5,2+,3,-28.985
    5,5,2+,4,-28.385
    5,5,2+,5,-27.560
    5,5,2-,0,-31.345
    5,5,2-,1,-28.976
    5,5,2-,2,-28.976
    5,5,2-,3,-28.705
    5,5,4+,0,-30.430
    5,5,4-,0,-29.895
    5,5,1-,0,-29.580
    5,5,1-,1,-29.027
    5,5,1-,2,-28.644
    5,5,1-,3,-16.355
    PROJECTILES
    Z,N,JPi,i,ENE
    2,2,"0+",0,-28.296
    _________________________________
"""

