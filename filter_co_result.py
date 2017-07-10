#!/usr/bin/python

#usage: python filter_co_result.py result rho> filtered_result

import sys 
cut = float(sys.argv[2])
for line in open(sys.argv[1],'r'):
    spl = line.strip().split('\t')
    rho = float(spl[3])
    if rho >= cut or rho <= -cut:
        print line,
