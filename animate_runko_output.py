import numpy as np
import sys, os
import time
from numpy import pi, log, log10, exp, sin, cos, tanh

import re

RElap = re.compile(r'------ lap: (\S+)')

#-------------------------------------------------- 
if __name__ == "__main__":

    n_start = 45
    n_dt    = 68


    with open('pcap.txt') as file:
        sim = file.readlines()

    
    ns = []
    for i in range(len(sim)):
        m2 = RElap.match(sim[i])
        if m2: ns.append(i)

    UP = "\x1B[3A"
    CLR = "\x1B[0K"

    for i in range(len(ns)-1):
        i1 = ns[i]
        i2 = ns[i+1]


        block = sim[i1:i2]
        #print(block[0])

        #os.system("clear && printf '\e[3J'")

        #sys.stdout.write("\x1b[1A\x1b[2K")
        for b in block: sys.stdout.write(b)
        #print(block)

        time.sleep(0.2)
        #os.system('clear')
        print(CLR)
        #input()




    #i1 = n_start
    #i2 = n_start + n_dt

    #lap = 0
    #while True:
    #    block = sim[i1:i2]

    #    #sys.stdout.write("\x1b[1A\x1b[2K")

    #    print('----------line ', i1, i2)
    #    print(block[0])

    #    #for b in block: sys.stdout.write(b)

    #    i1 += n_dt
    #    i2 += n_dt
    #    lap += 1

    #    time.sleep(0.1)

    #    if i2 > len(sim): break

    #    if lap > 5: break














