"""
    Created on August 2024 by Alan D.K. for 11C_Project

    This code increase the number of channels in the CC calculations.
    
    Comment on October 29: Outdated code, can be improved with an independent
    input file called ChannelTest.in when needed!!
"""

import numpy as np
import subprocess as sp


# Output file name
# now = datetime.now()
# optout = 'conttest%s.out'% (now.strftime("%y.%m.%d-%H:%M"))
# # Erease output file
# with open(optout,'w') as output:
#     pass
# Read file name
readfilename = 'IN2P3_11C_CC_GSMOpt-24.07.24-10.00_Basis-24.08.07-19.00_3I2-.in'

# Declaration of funcitons
def print_twice(*args,**kwargs): # Allow to print in file and in terminal at the same line
    print(*args,**kwargs)
    with open(optout,"a") as f:  # appends to file and closes it when finished
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
# Idea from https://stackoverflow.com/questions/176918/how-to-find-the-index-for-a-given-item-in-a-list
def indexof (obj, elem, offset=0):
    if elem in obj[offset:]:
        return offset + obj[offset:].index(elem)
    return -1
# .-

# Define the arrays with the channels (ordered by energy)
protonchannels_exp = ["  pole.target(yes)      3+  0          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      1+  0          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      0+  0          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      1+  1          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      2+  0          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      3+  1          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      2-  0          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      2+  1          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      1+  2          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      2+  2          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      4+  0          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      3-  0          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      4-  0          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      1-  0          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      3+  2          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      1-  1          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      2-  1          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      0+  1          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      2-  2          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      1-  2          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      3-  1          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 "]
protonchannels_gsm = ["  pole.target(yes)      0+  2          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      3-  2          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      4+  1          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      4+  2          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      4-  1          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      4-  2          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      5+  0          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      5+  1          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      5+  2          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      5-  0          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      5-  1          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 ",
                      "  pole.target(yes)      5-  2          proton             entrance.channel(no)       9                                p1/2 p3/2 d5/2 s1/2 d3/2 f7/2 f5/2 g9/2 g7/2 "]
protonchannels = protonchannels_exp + protonchannels_gsm

alphachannels_exp = ["  pole.target(yes)      3/2- 0         alpha              entrance.channel(no)       5                                S P D F G ",
                     "  pole.target(yes)      1/2- 0         alpha              entrance.channel(no)       5                                S P D F G ",
                     "  pole.target(yes)      7/2- 0         alpha              entrance.channel(no)       5                                S P D F G ",
                     "  pole.target(yes)      5/2- 0         alpha              entrance.channel(no)       5                                S P D F G ",
                     "  pole.target(yes)      5/2- 1         alpha              entrance.channel(no)       5                                S P D F G ",
                     "  pole.target(yes)      7/2- 1         alpha              entrance.channel(no)       5                                S P D F G "]
alphachannels_gsm = ["  pole.target(yes)      3/2- 1         alpha              entrance.channel(no)       5                                S P D F G ",
                     "  pole.target(yes)      1/2- 1         alpha              entrance.channel(no)       5                                S P D F G ",
                     "  pole.target(yes)      5/2+ 0         alpha              entrance.channel(no)       5                                S P D F G ",
                     "  pole.target(yes)      5/2+ 1         alpha              entrance.channel(no)       5                                S P D F G ",
                     "  pole.target(yes)      3/2+ 0         alpha              entrance.channel(no)       5                                S P D F G ",
                     "  pole.target(yes)      3/2+ 1         alpha              entrance.channel(no)       5                                S P D F G"]
alphachannels = alphachannels_exp + alphachannels_gsm

# Define limits and step
min_protonchannels = 6
min_alphachannels = 4
step_protonchannels = 4
step_alphachannels = 2

# Start calculation loop
loop1 = np.arange(min_protonchannels, len(protonchannels)+step_protonchannels/2, step_protonchannels,dtype=int)
loop2 = np.arange(min_alphachannels, len(alphachannels)+step_alphachannels/2, step_alphachannels,dtype=int)
# Proton loop
for n1 in loop1:
    for datain_a in [alphachannels_exp,alphachannels]:
        with open(readfilename, 'r') as readfile:
            data = readfile.read().split('\n')
            #
        # Search the channels line
        theline = searchline(readfilename,"target.projectile.state(s)")
        thelineend = searchline(readfilename,"CC.cluster.Berggren.basis")
        # Prepare data
        print("Number of proton channels : %s"% len(protonchannels[:n1]))
        print("Number of alpha channels : %s"% len(datain_a))
        datain = protonchannels[:n1] + datain_a
        data[theline] = "%s target.projectile.state(s)"% len(datain)
        data[theline+2:thelineend-1] = datain
        inputfile_aux = '\n'.join(data)
        with open(readfilename,'w') as fileout:
            fileout.write(inputfile_aux)
        #
        sp.run(['mpirun -np 4 ./CC_exe < '+readfilename+' >> IN2P3_11C_CC_GSMOpt-24.07.24-10.00_Basis-24.08.07-19.00_3I2-_24.08.08-11.00-ChannelsTest.out'], shell=True)
        #
