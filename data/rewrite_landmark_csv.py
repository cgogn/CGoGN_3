import pandas as pd
import os
from os import listdir
from os.path import isfile, join
import sys

filepath = sys.argv[1]

# filepath = os.path.dirname(os.path.abspath(__file__)) + filepath

print(filepath)

onlyfiles = [f for f in listdir(filepath) if isfile(join(filepath, f))]

onlyfiles.sort()

for f in onlyfiles : 
    if f.endswith('.csv') :
        print("PYTHON")
        print(f)
        dataframe = pd.read_csv(filepath + f , sep=',')
        if(dataframe.columns.size == 1) : 
            dataframe = pd.read_csv(filepath + f , sep=';')

        csv = False
        mv_columns = [col for col in dataframe.columns if (col.startswith("eye") or col.startswith("gaze") or col.startswith("p") or col.startswith("AU") or col.startswith("Z_") or col.startswith("Y_") or col.startswith("X_"))]
        print(mv_columns)
        if mv_columns :
            csv = True
            print("LA")
            dataframe = dataframe.drop(columns=mv_columns)
            print(dataframe.columns)
            
        if(csv): 
            dataframe.to_csv(join(filepath,f) , index=None , sep=',')