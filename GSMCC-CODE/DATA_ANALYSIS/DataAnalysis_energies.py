"""
    Created on June 2024 by Alan D.K.

    This code search for energies of the GSMCC outputs.
    The code search on each files for the JPi, pole.energies
    and complete.energies and then save them in a one single
    file ordered by input file and JPi.
    
    The name of the file to be analyzed must be defined in the 
    file "input.Data_Analysis_Energies". 
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
inputfile = 'input.Data_Analysis_Energies' 
with open(inputfile, 'r') as readfile:
    dataaux = readfile.read().split('\n')
aux_l = searchline(inputfile, "GSMCC-NFILES:")
numberfiles = int(dataaux[aux_l+1])
aux_l = searchline(inputfile, "GSMCC-FILE:")
namefiles = []
for j in range(aux_l+1,aux_l+1+numberfiles):
    namefiles.append( dataaux[j] )
aux_l = searchline(inputfile, "JPI-STATES:")
all_jpi = dataaux[aux_l+1].split(',')
aux_l = searchline(inputfile, "OUTPUT-FILE:")
if len(all_jpi) == 1:
    outputfile = dataaux[aux_l+1]+'_'+all_jpi[0].replace('/','I')
else:
    outputfile = dataaux[aux_l+1]
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
data_full = []
for i in range(0,numberfiles):
    namefile = namefiles[i]
    print(namefile)
    # Reading data
    with open(namefile, 'r') as readfile:
        data = readfile.read().split('\n')
        #
    #
    data_complete = []
    line_all = searchline_all(namefile,"resonant state")
    # Analyzing each jpi
    for j in range(0, len(all_jpi)):
        # Identify Jpi, index, energies and withs of the states:
        aux_save = []
        for line in line_all:
            jpi_aux = data[line].split()[0]
            if jpi_aux == all_jpi[j]:
                aux_index = int(data[line].split('(')[1].split(')')[0])
                aux_energies = float(data[line+10].split(' : ')[1].split()[0])
                aux_widths = float(data[line+11].split(' : ')[1].split()[0])
                aux_energies_pole = float(data[line+6].split(' : ')[1].split()[0])
                aux_widths_pole = float(data[line+7].split(' : ')[1].split()[0])
                aux_save.append( [aux_index, aux_energies, aux_energies_pole, aux_widths, aux_widths_pole] )
        data_complete.append(aux_save)
    #
    data_full.append(data_complete)
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
    f.write('# F   JPi  i {0:>16} {1:>16} {2:>16} {3:>16}\n'.format('ENE(MeV)','ENE_pole(MeV)','GAMMA(keV)','GAMMA_pole(keV)'))
for j in range(0, len(all_jpi)):
    aux_print_sort_1D = ''
    for i in range(0, numberfiles):
        aux_data = data_full[i]
        aux_print = aux_data[j]
        aux_print_sort = sorted(aux_print, key=lambda x: x[1])
        for item in aux_print_sort:
            aux_print_sort_1D += ' {0:2d}  {1:s}'.format(i,all_jpi[j])
            for k in range(0,len(item)): 
                num = item[k]
                if k == 0:
                    aux_print_sort_1D += ' {0:2d}'.format(num)
                else:
                    aux_print_sort_1D += ' {0:16.8f}'.format(num)
            aux_print_sort_1D += '\n'
        aux_print_sort_1D += '\n'
    aux_print_sort_1D += '\n\n'
    with open(outputfile,"a") as f:  # appends to file and closes it when finished
        f.write(aux_print_sort_1D)
        

"""
    Example of input.Data_Analysis_Energies file:
    _________________________________
    GSMCC-NFILES: # Number of files
    2
    
    GSMCC-FILE: # Name of each file
    IN2P3_11C_CC_GSMOpt-24.08.26-11.00_Basis-24.10.30-15.00_ICF-24.11.04-16.30.out
    IN2P3_11C_CC_GSMOpt-24.08.26-11.00_Basis-24.10.30-15.00_ICF-24.11.04-19.30.out
    
    OUTPUT-FILE:
    IN2P3_11C_CC_GSMOpt-24.08.26-11.00_Basis-24.10.30-15.00.Energies
    
    JPI-STATES: # Comma separated
    3/2-,5/2-,5/2+,7/2+
    _________________________________
    
"""