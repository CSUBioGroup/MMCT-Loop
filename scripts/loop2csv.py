#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import sys
import csv

id = sys.argv[1]

data = open('results_%s.loop'%(id),'r')
heads = data.readline().strip().split('\t')
sum_sig = 0
sum_es = 0
with open('results_%s_loop.csv'%(id), 'w') as f:
    writer = csv.writer(f)
    writer.writerow(heads)
    for line in data:
        r = line.strip().split('\t')
        if r[5] != 'inf':
            if float(r[-1]) == 1 :
                sum_sig += 1
                sum_es += float(r[5])
        writer.writerow(r)
    print('significant = %d\tES = %lf\n'%(sum_sig, sum_es/sum_sig))
