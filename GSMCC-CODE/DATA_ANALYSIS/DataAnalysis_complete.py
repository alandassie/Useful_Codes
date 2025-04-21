"""
    Created on June 2024 by Alan D.K.

    This is a multitask code used for the analysis of the GSMCC outputs.
    The analysis done by the code are the followings:
    + sum up the partial widths for each mass partition.
    + sum up the occ prob for each mass partition 
    and shown those greater than 5% of the total.
    
    The name of the file to be analyzed must be defined in the 
    file "Data_Analysis.in". An example of the input file is at the end of the code.
"""
import ast

# Input file
inputfile = 'Data_Analysis.in' 
with open(inputfile, 'r') as readfile:
    dataaux = readfile.read().split('\n')
namefile = dataaux[1]
jpi_state = dataaux[3]
index_state = dataaux[5]
calc_sf = dataaux[7]
print_info = dataaux[9]
# Output file
jpi_state_write = jpi_state.replace('/','I')
outputfile = namefile[:-4] + '.DataAnalysis-' + jpi_state_write
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
with open(namefile, 'r') as readfile:
    data = readfile.read().split('\n')
    #
#
# Defining the lines of data
if calc_sf == 'YES':
    if index_state == 'ALL':
        print('ERROR! INDEX CANNOT BE ALL FOR SF CALCULATIONS')
        exit()
    else:
        index_values = ast.literal_eval(index_state)
        if type(index_values) == list:
            line_all = []
            for i in index_values:
                aux = searchline(namefile,"%s (%s) resonant state"% (jpi_state,i))
                line_all += aux
        else:
            aux = searchline(namefile,"%s (%s) resonant state"% (jpi_state,index_state))
            line_all = [aux]
else:
    if index_state == 'ALL':
        aux = searchline_all(namefile,"resonant state")
        line_all = aux
    else:
        index_values = ast.literal_eval(index_state)
        if type(index_values) == list:
            line_all = []
            for i in index_values:
                aux = searchline_all(namefile,"%s (%s) resonant state"% (jpi_state,i))
                if type(line_all) == list:
                    line_all += aux
                else:
                    line_all.append(aux)
        else:
            line_all = []
            aux = searchline_all(namefile,"%s (%s) resonant state"% (jpi_state,index_state))
            if type(aux) == list:
                line_all += aux
            else:
                line_all = [aux]
# Identify index, energy and width of the state:
indexes = []
energies = []
widths = []
for line in line_all:
    indexes.append( data[line].split('(')[1].split(')')[0] )
    energies.append( data[line+10].split(' : ')[1] )
    widths.append( data[line+11].split(' : ')[1] )
#
line_occprob_i = []
line_occprob_f = []
line_partwidth_rd_i = []
line_partwidth_rd_f = []
line_partwidth_ri_i = []
line_partwidth_ri_f = []
k = -1
for line in line_all:
    k += 1
    # Search first occ probabilities
    i = 0
    while i < 1000:
        aux = data[line+i]
        finder = aux.find('occupation probability :')
        if finder == -1:
            i += 1
            continue
        else:
            line_occprob_i_aux = line+i
            break
    # Search last occ probabilities
    i = 0
    while i < 1000:
        aux = data[line_occprob_i_aux+i]
        finder = aux.find('occupation probability :')
        if finder != -1:
            i += 1
            continue
        else:
            line_occprob_f_aux = line_occprob_i_aux+i
            break
    # Search first partial width (r-dependent)
    i = 0
    while i < 1000:
        aux = data[line_occprob_f_aux+i]
        finder = aux.find('(r-dependent)')
        if finder == -1:
            i += 1
            continue
        else:
            line_partwidth_rd_i_aux = line_occprob_f_aux+i+3
            break
    # Search last partial width (r-dependent)
    i = 0
    while i < 1000:
        aux = data[line_partwidth_rd_i_aux+i]
        finder = aux.find(' Gamma : ')
        if finder != -1:
            i += 1
            continue
        else:
            line_partwidth_rd_f_aux = line_partwidth_rd_i_aux+i
            break
    # Search first partial width (r-independent)
    i = 0
    while i < 1000:
        aux = data[line_occprob_f_aux+i]
        finder = aux.find('(r-independent)')
        if finder == -1:
            i += 1
            continue
        else:
            line_partwidth_ri_i_aux = line_occprob_f_aux+i+3
            break
    # Search last partial width (r-independent)
    i = 0
    while i < 1000:
        aux = data[line_partwidth_ri_i_aux+i]
        finder = aux.find(' Gamma : ')
        if finder != -1:
            i += 1
            continue
        else:
            line_partwidth_ri_f_aux = line_partwidth_ri_i_aux+i
            break
    #
    line_occprob_i.append(line_occprob_i_aux)
    line_occprob_f.append(line_occprob_f_aux)
    line_partwidth_rd_i.append(line_partwidth_rd_i_aux)
    line_partwidth_rd_f.append(line_partwidth_rd_f_aux)
    line_partwidth_ri_i.append(line_partwidth_ri_i_aux)
    line_partwidth_ri_f.append(line_partwidth_ri_f_aux)
#
if calc_sf == 'YES':
    # Search for spectroscopic factor lines
    sf_lines = searchline_all(namefile,'Spectroscopic factor ')
    # Pick only the lines with the desired JPi_index
    sf_lines_jpi = []
    for i in sf_lines:
        aux = data[i]
        read_jpi = aux.split('--> ')[1].split()[0].split('(')[0]
        read_index = aux.split('--> ')[1].split()[0].split('(')[1].split(')')[0]
        if read_jpi == jpi_state and read_index == index_state:
            sf_lines_jpi.append(aux)
# .-
# .-
#
# Starting data analysis
print_twice('DATA ANALYSIS CODE BUILT BY ALAN D.K.')
print_twice(37*'-')
print_twice('Data in the output file:')
print_twice(namefile+'\n\n')
print_twice(20*'-')
print_twice('Start the data analysis for the state %s with index %s'% (jpi_state,index_state))
print_twice('Number of calculations for this state %s'% (len(energies)))
k = -1
for index, energy, width in zip(indexes,energies,widths):
    k += 1
    if print_info == 'ALL':
        print_twice(160*'-')
        print_twice('{0:2d}  Energy : {1:<.10f} \n    Width  : {2:<.10f}'.format(int(index),float(energy.split('MeV')[0]),float(width.split('keV')[0])))
    else:
        print_twice('{0:2d}   {1:2d}  {2:>14.8f} {3:>14.8f}'.format(k+1,int(index),float(energy.split('MeV')[0]),float(width.split('keV')[0])))
    # .-
    # Occupation probabilities
    if print_info == 'ALL': print_twice('Occupation probabilities:\n')
    #
    # Build two arrays
    channels = []
    projectile_array = []
    occprob_real = []
    occprob_imag = []
    partial_occprob_real = []
    partial_occprob_imag = []
    occprob_real_aux = 0
    occprob_imag_aux = 0
    j = 0
    i = -1
    projectile_old = " "
    line_occprob_i_aux = line_occprob_i[k]
    line_occprob_f_aux = line_occprob_f[k]
    for x in data[line_occprob_i_aux:line_occprob_f_aux]:
        if x == '':
            continue
        projectile = x.split(' : ')[2].split(' ')[0]
        channel = x.split('^(%s'% jpi_state)[0].split('channel ')[1]
        aux1 = x.split(' : ')[3].split(',')
        aux_real = float(aux1[0].split('(')[1])
        aux_imag = float(aux1[1].split(')')[0])
        if projectile == projectile_old :
            occprob_real_aux += aux_real
            occprob_imag_aux += aux_imag
            partial_occprob_real[i].append(aux_real)
            partial_occprob_imag[i].append(aux_imag)
            channels[i].append(channel)
        else:
            if j == 0:
                i += 1
                occprob_real_aux += aux_real
                occprob_imag_aux += aux_imag
                projectile_old = projectile
                partial_occprob_real.append([aux_real])
                partial_occprob_imag.append([aux_imag])
                channels.append([channel])
                projectile_array.append(projectile)
                j +=1
            else:
                occprob_real.append(occprob_real_aux)
                occprob_imag.append(occprob_imag_aux)
                #
                i+=1
                projectile_old = projectile
                projectile_array.append(projectile)
                occprob_real_aux = aux_real
                occprob_imag_aux = aux_imag
                partial_occprob_real.append([aux_real])
                partial_occprob_imag.append([aux_imag])
                channels.append([channel])
    occprob_real.append(occprob_real_aux)
    occprob_imag.append(occprob_imag_aux)
    
    if print_info == 'ALL':
        for i in range(0,len(projectile_array)):
            print_twice('Projectile mass partition : {0:s}'.format(projectile_array[i])) 
            print_twice('Total Occupation Probability : ({0:.5f},{1:.5f})'.format(occprob_real[i],occprob_imag[i]))
            print_twice('Max real Occupation Probability : {0:.5f}'.format(max(partial_occprob_real[i])))
            index_max = max(range(len(partial_occprob_real[i])), key=partial_occprob_real[i].__getitem__)
            print_twice('With the channel : %s'% channels[i][index_max])
            print_twice('  +Occupation probabilities greater than 5%:')
            for j in range(0,len(partial_occprob_real[i])):
                if abs(partial_occprob_real[i][j]) > 0.05:
                    print_twice('    {0:s} -> ({1:.5f},{2:.5f})'.format(channels[i][j],partial_occprob_real[i][j],partial_occprob_imag[i][j]))
            print_twice(10*'-')
            print_twice(' ')
    # .-
    if print_info == 'ALL': print_twice(20*'-')
    # .-
    # Partial Widths (r-dependent)
    if print_info == 'ALL': print_twice('Partial Widths (r-dependent):\n')
    #
    # Build arrays
    channels = []
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
    line_partwidth_i_aux = line_partwidth_rd_i[k]
    line_partwidth_f_aux = line_partwidth_rd_f[k]
    for x in data[line_partwidth_i_aux:line_partwidth_f_aux]:
        if x == '':
            continue
        projectile = x.split(' : ')[2].split(' ')[0]
        channel = x.split('^(%s'% jpi_state)[0]
        aux1 = x.split(' : ')[3].split(' , ')
        aux2 = float(aux1[0].split()[0])
        aux3 = float(aux1[0].split()[3])
        if projectile == projectile_old :
            current_width_aux += aux2
            altcurrent_width_aux += aux3
            partial_current_width[i].append(aux2)
            partial_altcurrent_width[i].append(aux3)
            channels[i].append(channel)
        else:
            if j == 0:
                i += 1
                current_width_aux += aux2
                altcurrent_width_aux += aux3
                projectile_old = projectile
                partial_current_width.append([aux2])
                partial_altcurrent_width.append([aux3])
                channels.append([channel])
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
                channels.append([channel])
    current_width.append(current_width_aux)
    altcurrent_width.append(altcurrent_width_aux)

    if print_info == 'ALL':     
        for i in range(0,len(projectile_width)):
            print_twice('Projectile mass partition : {0:s}'.format(projectile_width[i]))
            print_twice('Total current width (keV) : {0:.5g}'.format(current_width[i]))
            print_twice('Max partial current width (keV) : {0:.5g}'.format(max(partial_current_width[i])))
            index_max = max(range(len(partial_current_width[i])), key=partial_current_width[i].__getitem__)
            print_twice('  With the channel : %s'% channels[i][index_max])
            print_twice('Total alt current width (keV) : {0:.5g}'.format(altcurrent_width[i]))
            print_twice('Max partial alt current width (keV) : {0:.5g}'.format(max(partial_altcurrent_width[i])))
            index_max = max(range(len(partial_altcurrent_width[i])), key=partial_altcurrent_width[i].__getitem__)
            print_twice('  With the channel : %s'% channels[i][index_max])
            print_twice(10*'-')
            print_twice(' ')
    # .-
    # Partial Widths (r-independent)
    if print_info == 'ALL': print_twice('Partial Widths (r-independent):\n')
    #
    # Build arrays
    channels = []
    projectile_width = []
    current_width = []
    current_lower_width = []
    current_upper_width = []
    partial_current_width = []
    partial_current_lower_width = []
    partial_current_upper_width = []
    current_width_aux = 0
    current_width_lower_aux = 0
    current_width_upper_aux = 0
    j = 0
    i = -1
    projectile_old = " "
    line_partwidth_i_aux = line_partwidth_ri_i[k]
    line_partwidth_f_aux = line_partwidth_ri_f[k]
    for x in data[line_partwidth_i_aux:line_partwidth_f_aux]:
        if x == '':
            continue
        projectile = x.split(' : ')[2].split(' ')[0]
        channel = x.split('^(%s'% jpi_state)[0]
        aux1 = x.split(' : ')[3]
        aux2 = float(aux1.split()[0])
        aux3 = float(aux1.split()[3])
        aux4 = float(aux1.split()[9])
        if projectile == projectile_old :
            current_width_aux += aux2
            current_width_lower_aux += aux3
            current_width_upper_aux += aux4
            partial_current_width[i].append(aux2)
            partial_current_lower_width[i].append(aux3)
            partial_current_upper_width[i].append(aux4)
            channels[i].append(channel)
        else:
            if j == 0:
                i += 1
                current_width_aux += aux2
                current_width_lower_aux += aux3
                current_width_upper_aux += aux4
                projectile_old = projectile
                partial_current_width.append([aux2])
                partial_current_lower_width.append([aux3])
                partial_current_upper_width.append([aux4])
                channels.append([channel])
                projectile_width.append(projectile)
                j +=1
            else:
                current_width.append(current_width_aux)
                current_lower_width.append(current_width_lower_aux)
                current_upper_width.append(current_width_upper_aux)
                #
                i+=1
                projectile_old = projectile
                projectile_width.append(projectile)
                current_width_aux = aux2
                current_width_lower_aux = aux3
                current_width_upper_aux = aux4
                partial_current_width.append([aux2])
                partial_current_lower_width.append([aux3])
                partial_current_upper_width.append([aux4])
                channels.append([channel])
    current_width.append(current_width_aux)
    current_lower_width.append(current_width_lower_aux)
    current_upper_width.append(current_width_upper_aux)

    if print_info == 'ALL':     
        for i in range(0,len(projectile_width)):
            print_twice('Projectile mass partition : {0:s}'.format(projectile_width[i]))
            print_twice('Total app. current width (keV) : {0:.5g}'.format(current_width[i]))
            index_max = max(range(len(partial_current_width[i])), key=partial_current_width[i].__getitem__)
            print_twice('Max partial app. current width (keV) : {0:.5g} ({1:.5g},{2:.5g})'.format( 
                        max(partial_current_width[i]), partial_current_lower_width[i][index_max], partial_current_upper_width[i][index_max] ))
            print_twice('  With the channel : %s'% channels[i][index_max])
            print_twice(10*'-')
            print_twice(' ')
    # .-
    # .-
    # Spectroscopic factors:
    if calc_sf == 'YES':
        if print_info == 'ALL': print_twice(20*'-')
        # .-
        # Partial Widths
        if print_info == 'ALL': print_twice('Spectroscopic factors:\n')
        #
        headers = ['Projectile', 'P state', 'T state', 'Non-Anti S',' ','Anti S',' ']
        projectile = []
        p_state = []
        t_state = []
        nonanti_s_r = []
        nonanti_s_i = []
        anti_s_r = []
        anti_s_i = []
        for text in sf_lines_jpi:
            projectile.append(text.split(' -->')[0].split()[-1])
            p_state.append(text.split(' -->')[0].split()[-2])
            t_state.append(text.split(' -->')[0].split()[-4])
            nonanti_s_r.append(float(text.split('Non-antisymmetrized ')[1].split(',')[0].split('(')[1]))
            nonanti_s_i.append(float(text.split('Non-antisymmetrized ')[1].split(',')[1].split(')')[0]))
            anti_s_r.append(float(text.split('antisymmetrized ')[-1].split(',')[0].split('(')[1]))
            anti_s_i.append(float(text.split('antisymmetrized ')[-1].split(',')[1].split(')')[0]))
        # Compare lengths
        max_p = len(max(max(projectile,key=len),headers[0],key=len))
        max_p_state = len(max(max(p_state,key=len),headers[1],key=len))
        max_t_state = len(max(max(t_state,key=len),headers[2],key=len))
        lista = [max_p,max_p_state,max_t_state,7,7,7,7]
        if print_info == 'ALL':
            for i in range(0,len(headers)):
                col = headers[i]
                print_twice(col.ljust(lista[i]), end="")
                print_twice('  ', end = "")
            print_twice()
            for i in range(0,len(projectile)):
                print_twice(projectile[i].ljust(lista[0]), end="")
                print_twice('  ', end = "")
                print_twice(p_state[i].ljust(lista[1]), end="")
                print_twice('  ', end = "")
                print_twice(t_state[i].ljust(lista[2]), end="")
                print_twice('  ', end = "")
                print_twice('{0:7.5f} '.format(nonanti_s_r[i]), end="")
                print_twice('{0:8.5f} '.format(nonanti_s_i[i]), end="")
                print_twice('    ', end = "")
                print_twice('{0:7.5f} '.format(anti_s_r[i]), end="")
                print_twice('{0:8.5f} '.format(anti_s_i[i]))
        
print_twice(' ')
print_twice(30*'-')
"""
    Example of Data_Analysis.in file:
    _________________________________
    GSMCC-FILE:
    IN2P3_11C_CC_GSMOpt-24.08.26-11.00_Basis-24.10.30-15.00_ICF-24.11.04-16.30.out
    JPi-STATE:
    3/2-
    INDEX: # SHOULD BE A NUMBER (ex. "2"), A LIST OF NUMBERS (ex. "[0,2,3]") OR "ALL"
    3
    SF-CALC:
    YES
    PRINT: ALL or ENERGIES
    ENERGIES
    _________________________________
    
"""