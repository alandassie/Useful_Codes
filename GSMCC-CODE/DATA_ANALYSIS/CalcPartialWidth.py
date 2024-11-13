"""
    Created on June 2024 by Alan D.K.

    This is a multitask code used for the analysis of the GSMCC outputs.
    The analysis done by the code are the followings:
    + sum up the partial widths for each mass partition 
    and shown those greater than 5% of the total.
    + sum up the occ prob for each mass partition 
    and shown those greater than 5% of the total.
    
    The name of the file to be analyzed must be defined in the 
    file "Data_Analysis.in". An example of the input file is at the end of the code.
"""

# Input file
inputfile = 'Data_Analysis.in' 
with open(inputfile, 'r') as readfile:
    dataaux = readfile.read().split('\n')
namefile = dataaux[1]
jpi_state = dataaux[3]
index_state = dataaux[5]
# Output file
jpi_state_write = jpi_state.replace('/','I')
outputfile = namefile[:-4] + '.DataAnalysis-' + jpi_state_write + '_' + index_state
#
# Declaration of funcitons
# def erease_output_file():
#     with open(logfile,'w') as output:
#         pass
def print_twice(*args,**kwargs): # Allow to print in file and in terminal at the same line
    print(*args,**kwargs)
    with open(outputfile,"a") as f:  # appends to file and closes it when finished
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
# Reading data
with open(readfilename, 'r') as readfile:
    data = readfile.read().split('\n')
    #
#
# Defining the lines of data
line = searchline(namefile,"%s (%s) resonant state"% (jpi_state,index_state))
# Search first occ probabilities
i = 0
while i < 1000:
    aux = data[line+i]
    finder = aux.find('occupation probability :')
    if finder == -1:
        i += 1
        continue
    else:
        line_occprob_i = line+i
        break
# Search last occ probabilities
i = 0
while i < 1000:
    aux = data[line_occprob_i+i]
    finder = aux.find('occupation probability :')
    if finder != -1:
        i += 1
        continue
    else:
        line_occprob_f = line_occprob_i+i
        break
# Search first partial width
i = 0
while i < 1000:
    aux = data[line_occprob_f+i]
    finder = aux.find('(alternative current formula)')
    if finder == -1:
        i += 1
        continue
    else:
        line_partwidth_i = line_occprob_f+i
        break
# Search last partial width
i = 0
while i < 1000:
    aux = data[line_occprob_f+i]
    finder = aux.find('(alternative current formula)')
    if finder != -1:
        i += 1
        continue
    else:
        line_partwidth_f = line_partwidth_i+i
        break

    
    
    
    
    
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
    Example of GSMCC_DCS_ene-ang.in file:
    _________________________________
    GSMCC-FILE
    IN2P3_11C_CC_GSMOpt-24.08.26-11.00_Basis-24.10.30-15.00_ICF-24.11.04-16.30_SF-24.11.13-13.00-3I2-.out
    JPi-STATE
    3/2-
    INDEX
    3
    _________________________________
    
"""