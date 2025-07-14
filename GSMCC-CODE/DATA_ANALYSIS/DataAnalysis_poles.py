"""
    Created on June 2024 by Alan D.K.

    This code search for energies of the poles from GSMCC outputs.
    The code search on each files for the desired projectile
    and look for pole.energies and k values
    and then save them in a one single
    file ordered by input file and quantum numbers.
    
    The name of the file to be analyzed must be defined in the 
    file "input.Data_Analysis_Poles". 
    An example of the input file is at the end of the code.
"""

#
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
# FILES
# Input file
inputfile = 'input.Data_Analysis_Poles' 
with open(inputfile, 'r') as readfile:
    dataaux = readfile.read().split('\n')
aux_l = searchline(inputfile, "GSMCC-NFILES:")
numberfiles = int(dataaux[aux_l+1])
aux_l = searchline(inputfile, "GSMCC-FILE:")
namefiles = []
for j in range(aux_l+1,aux_l+1+numberfiles):
    namefiles.append( dataaux[j] )
aux_l = searchline(inputfile, "PROJECTILE:")
all_projectile = dataaux[aux_l+1]
aux_pw = searchline(inputfile, "PARTIAL_WAVES:")
if aux_pw is not None:
    all_partial_waves = dataaux[aux_pw+1].split(',')
aux_l = searchline(inputfile, "OUTPUT-FILE:")
outputfile = dataaux[aux_l+1]+'_'+all_projectile
if aux_pw is not None:
    outputfile += '-'+'.'.join(all_partial_waves)
#
# CONTINUE FUNDEF
def print_twice(*args,**kwargs): # Allow to print in file and in terminal at the same line
    print(*args,**kwargs)
    with open(outputfile,"a") as f:  # appends to file and closes it when finished
        print(file=f,*args,**kwargs)
# .-
def erease_output_file():
    with open(outputfile,'w') as output:
        pass
# .-
erease_output_file()
# .-
# Write information on output file
print_twice('#  INPUT FILES:  ')
for i in range(0, numberfiles):
    print_twice('#    {0:2d}    {1:<s}'.format(i, namefiles[i]))
print_twice('########################################\n')
# Start data aquisition
n_poles = 0
data_full = []
for i in range(0,numberfiles):
    namefile = namefiles[i]
    print(namefile)
    # Reading data
    with open(namefile, 'r') as readfile:
        data = readfile.read().split('\n')
        #
    #
    line_all = searchline_all(namefile,"  k : ")
    # Analyzing the specific projectile
    # Identify Jpi, index, energies and withs of the states:
    aux_save = []
    for line in line_all:
        projectile_aux = data[line].split()[0]
        if projectile_aux == all_projectile:
            if i == 0:
                n_poles += 1
            aux_jpi    = data[line].split()[1]
            aux_k_real = float(data[line].split(':')[1].split('(')[1].split(')')[0].split(',')[0])
            aux_k_imag = float(data[line].split(':')[1].split('(')[1].split(')')[0].split(',')[1])
            aux_energy = float(data[line].split(':')[2].split()[0])
            aux_widths = float(data[line].split(':')[3].split()[0])
            aux_jostff = float(data[line].split(':')[-1])
            aux_save.append( [aux_jpi, aux_k_real, aux_k_imag, aux_energy, aux_widths, aux_jostff] )
    #
    data_full.append(aux_save)
    # # Writing process
    # for j in range(0, len(all_jpi)):
    #     aux_print = data_complete[j]
    #     aux_print_sort = sorted(aux_print, key=lambda x: x[1])
    #     aux_print_sort_1D = ''
    #     for item in aux_print_sort:
    #         aux_print_sort_1D += ' {0:2d}  {1:s}'.format(i,all_jpi[j])
    #         for num in item:
    #             aux_print_sort_1D += ' ' + str(num)
    #         aux_print_sort_1D += '\n'
    #     aux_print_sort_1D += '\n\n'
    #     with open(outputfile,"a") as f:  # appends to file and closes it when finished
    #         f.write(aux_print_sort_1D)

# Writing process
with open(outputfile,"a") as f:  # appends to file and closes it when finished
    f.write('# F   JPi  {0:>16} {1:>16} {2:>16} {3:>16}\n'.format('k.real(fm^-1)','k.imag(fm^-1)','Ene(MeV)','Gamma(keV)','|Jost|'))
aux_print_sort_1D = ''
for i in range(0, numberfiles):
    aux_old = '0'
    for j in range(0,n_poles):
        aux_data = data_full[i][j]
        jpi = aux_data[0]
        for k in jpi:
            if k.isdigit() is not True:
                l = k
                break
        if aux_pw is not None:
            if l not in all_partial_waves:
                continue
        else:
            if l != aux_old:
                aux_print_sort_1D += '\n'
        aux_print_sort_1D += ' {0:2d}  {1:s} {2:16.8f} {3:16.8f} {4:16.8f} {5:16.8f} \n'.format(i,jpi,aux_data[1],aux_data[2],aux_data[3],aux_data[4],aux_data[5])
        aux_old = aux_data[0][-1]
    aux_print_sort_1D += '\n'
with open(outputfile,"a") as f:  # appends to file and closes it when finished
    f.write(aux_print_sort_1D)
        

"""
    Example of input.Data_Analysis_Poles file:
    _________________________________
    GSMCC-NFILES: # Number of files
    2
    
    GSMCC-FILE: # Name of each file
    IN2P3_11C_CC_GSMOpt-24.08.26-11.00_Basis-24.10.30-15.00_ICF-24.11.04-16.30.out
    IN2P3_11C_CC_GSMOpt-24.08.26-11.00_Basis-24.10.30-15.00_ICF-24.11.04-19.30.out
    
    OUTPUT-FILE:
    IN2P3_11C_CC_GSMOpt-24.08.26-11.00_Basis-24.10.30-15.00.Energies
    
    PROJECTILE:
    proton
    
    PARTIAL_WAVES: # Optional, if you want to analyze only some partial waves, comma separated
    S
    _________________________________
    
"""