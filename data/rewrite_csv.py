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
        print(f)
        dataframe = pd.read_csv(filepath + f , sep=',')
        if(dataframe.columns.size == 1) : 
            dataframe = pd.read_csv(filepath + f , sep=';')
        csv = False
        if(' AU28_c' in dataframe.columns): 
            dataframe = dataframe.drop(columns=' AU28_c')
            csv = True
        if('AU28_c' in dataframe.columns): 
            dataframe = dataframe.drop(columns='AU28_c')
            csv = True
        if('face_id' in dataframe.columns): 
            dataframe = dataframe.drop(columns='face_id')
            csv = True
        if('confidence' in dataframe.columns): 
            dataframe = dataframe.drop(columns='confidence')
            csv = True
        if('success' in dataframe.columns): 
            dataframe = dataframe.drop(columns='success')
            csv = True

        if(csv) : 
            dataframe.to_csv(join(filepath,f) , index=None , sep=',')
