import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from os import listdir
from os.path import isfile, join , basename

filepath = "data/meshes/AUs/buste_no_brows_clean/CSV/"

onlyfiles = [f for f in listdir(filepath) if isfile(join(filepath, f))]

onlyfiles.sort()

file = open("data/meshes/AUs/buste_no_brows/CSV_VIDEO/weights.txt", "w")   # 'r' for reading and 'w' for writing

for f in onlyfiles : 
    if f.endswith('.csv') :
        print(f)
        name_au = basename(f)
        name_au = name_au[:-4] + "_r"
        dataframe = pd.read_csv(filepath + f)
        fig = plt.figure()
        x = dataframe.index.to_numpy()
        y = dataframe[name_au].to_numpy()

        tmp = 0
        poids = 0.
        new_x = []
        for i in range(dataframe.index.size) :
            if (dataframe[name_au][i] >= 1.000) and (tmp <= 2) :
                print(name_au + " time : " + str(i) + " poids : " + str(poids))
                if(tmp == 1): 
                    file.write(name_au + " time : " + str(i) + " poids : " + str(poids) + "\n")
                tmp += 1
            new_x.append(poids)
            poids += 0.003

        new_x = pd.DataFrame(new_x , columns=['weights'])
        new_x = new_x['weights'].to_numpy()

        tmp = 0
        poids = 0.
        for i in range(dataframe.index.size) :
            if (dataframe[name_au][i] <= 1.000) and (tmp <= 2) :
                print("Reverse " + name_au + " time : " + str(i) + " poids : " + str(poids))
                if(tmp == 1): 
                    file.write("Reverse " + name_au + " time : " + str(i) + " poids : " + str(poids) + "\n")
                tmp += 1
            poids += 0.003

        x_smooth = np.linspace(new_x.min(), new_x.max(), 1500)
        spline = make_interp_spline(new_x, y, k=3)
        y_smooth = spline(x_smooth)

        plt.plot(x_smooth, y_smooth, label=f'{name_au} (smoothed)')
        plt.xlabel('Weight CGoGN')
        plt.ylabel('Intensity detected by OpenFace')
        plt.title(f'{name_au} Intensity over Time')
        plt.legend()
        #plt.show()
        plt.savefig(name_au + '.png')

file.close()  