import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from os import listdir
from os.path import isfile, join , basename, exists
import sys

filepath = sys.argv[1]

onlyfiles = [f for f in listdir(filepath) if isfile(join(filepath, f))]

onlyfiles.sort()

for f in onlyfiles : 
    name = basename(f)
    name = name[:-4]
    if f.endswith('Grimaces.csv') and (exists(filepath + name + "_static.csv")) :
        print(f)
        dataframe = pd.read_csv(filepath + f)
        dataframe_static = pd.read_csv(filepath + name + "_static.csv")

        df_columns = [col for col in dataframe.columns if (col.endswith("_r"))]

        for i in range(len(df_columns)) : 
            x = dataframe.index.to_numpy()
            y = dataframe[df_columns[i]].to_numpy()

            x_static = dataframe_static.index.to_numpy()
            y_static = dataframe_static[df_columns[i]].to_numpy()

            plt.ylim(-0.5, 5)
            plt.plot(x, y, label=f'{df_columns[i]} dynamic')
            plt.plot(x_static, y_static, label=f'{df_columns[i]} static')
            plt.xlabel('Frame')
            plt.ylabel('Intensity detected by OpenFace')
            plt.title(f'{df_columns[i]} Intensities dynamic and static')
            plt.legend()
            plt.savefig(filepath + name + ' ' + df_columns[i] + '.png')
            plt.close() 