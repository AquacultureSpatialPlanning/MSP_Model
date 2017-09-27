#!/usr/bin/env python
import os
import numpy as np
from scipy.io import loadmat  # this is the SciPy module that loads mat-files
import matplotlib.pyplot as plt
from datetime import datetime, date, time
import pandas as pd

print os.path.abspath(os.curdir)
os.chdir("..")
fid_mat = loadmat('Input/Data/Aqua_Dev_Indices.mat')['Aqua_Dev_Indices']
fid = fid_mat.tolist()

columns = []
for i in range(0,5):
    number = i + 1
    strg = "Scenario{}".format(number)
    columns.append(strg)

print columns


mat = loadmat('Input/Data/bc_set_policies.mat')  # load mat-file
df = pd.DataFrame(mat['bc_set_policies'].reshape(1061, 5)\
                       , index=fid, columns = columns)

df.to_csv('Output/Data/Revised_Seeding_Policies.csv', columns=columns, index = True, index_label='FID')
# df.to_csv()
