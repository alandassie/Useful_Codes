"""
    Created on February 2025 by Alan D.K.

    This code search for cross-sections of the GSMCC outputs.
    The code search on all the files of a particular folder
    the defined line name with the keyboard.
"""

import os

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
# Folder
folder = input("Enter the name of the folder:\n ")
print("Looking for files in %s"% folder)
outputfile = folder + '.crosssections'
# Looking for all files:
items = os.listdir(folder)
sorted_items = sorted(items)
print('List of files in the selected folder: ')
print(sorted_items)
# Name line
nameinline = input("Enter the line name (for ex. -total radiative capture cross section : antisymmetrized-):\n ")
print("Looking for lines with -%s-"% nameinline)
# .-
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
# Start data aquisition
data_full = []
for i in range(0,len(sorted_items)):
    namefile = sorted_items[i]
    print('File %s'% namefile)
    # Reading data
    with open(folder + '/' + namefile, 'r') as readfile:
        data = readfile.read().split('\n')
        #
    #
    theline = searchline(folder + '/' + namefile,'MeV(E.kinetic.total.system.CM)')
    energy = data[theline].split()[0]
    print('Energy = %s MeV'% energy)
    #
    theline = searchline(folder + '/' + namefile, nameinline)
    cs = data[theline].split('=')[1]
    #
    print_twice('  %s  %s '% (energy,cs))
    data_complete = []
    