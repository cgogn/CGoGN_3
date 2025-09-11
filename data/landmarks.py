from matplotlib import image
from matplotlib import pyplot as plt
import pandas as pd
import os
from os import listdir
from os.path import isfile, join
import sys

filepath = sys.argv[1]

# filepath = os.path.dirname(os.path.abspath(__file__)) + filepath

print(filepath)

data = image.imread('../build/stage/bin/Screenshot_0000.jpg')

dataframe = pd.read_csv(filepath , sep=',')

x_columns = [col for col in dataframe.columns if col.startswith("x_")]
y_columns = [col for col in dataframe.columns if col.startswith("y_")]

x_list = []
y_list = []

for i in range(len(dataframe.columns)):
    if(dataframe.columns[i].startswith("x_")) : 
        x_list.append(dataframe.at[0,dataframe.columns[i]])
    if(dataframe.columns[i].startswith("y_")) : 
        y_list.append(dataframe.at[0,dataframe.columns[i]])

print(x_list)
print(len(x_list))
print(len(y_list))

for i in range(len(x_list)) :     
    plt.plot(x_list[i], y_list[i] , marker='o' , color='red')

plt.imshow(data)
plt.show()
