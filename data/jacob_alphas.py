import pandas as pd



import numpy as np


import matplotlib.pyplot as plt


from scipy.interpolate import make_interp_spline


from os import listdir


from os.path import isfile, join, basename, exists


import sys





filepath = sys.argv[1]





onlyfiles = [f for f in listdir(filepath) if isfile(join(filepath, f))]





onlyfiles.sort()





file = open(sys.argv[2] + "slopes.txt", "w")





nb_files = [files for files in onlyfiles if (files.endswith('.csv') and files.startswith("AU"))]





file.write(str(len(nb_files)))


file.write("\n")





alphas = [[0 for x in range(len(nb_files))] for y in range(len(nb_files))]


betas = [[0 for x in range(len(nb_files))] for y in range(len(nb_files))]





row_jacob = 1





for f in nb_files :


    if(not (f.endswith('_static.csv'))) : 


        name = basename(f)


        name = name[:-4]


        if not((exists(filepath + name + "_static.csv"))) : 


            if f.endswith('.csv') :


                print(f)


                name_au = basename(f)


                name_au = name_au[:-4] + "_r"


                dataframe = pd.read_csv(filepath + f)


                columns = [col for col in dataframe.columns if (col.startswith("AU") and col.endswith("_r") and not(col.startswith("AU28")))]





                for i in range(len(columns)) :


                    x = dataframe.index.to_numpy()


                    y = dataframe[columns[i]].to_numpy()


                    poids = 0.





                    x0 = 0


                    y0 = 0


                    y1 = 0


                    x1 = 0


                    draw = True


                    


                    for j in range(dataframe.index.size) :                





                        if(j == 0) :


                            dataframe.loc[j,columns[i]] = dataframe[columns[i]][j+2]


                            dataframe.loc[j+1,columns[i]] = dataframe[columns[i]][j+2]


                            x0 = poids


                            y0 = dataframe.loc[0,columns[i]]


                        if(dataframe[columns[i]][j] == 0) and (dataframe[columns[i]][0] == 0) and (poids < 1) :


                            x0 = poids


                        if(poids >= x0 + 1.5) and (draw): 


                            draw = False


                            y1 = (dataframe.loc[j,columns[i]] + dataframe.loc[j+1,columns[i]] + dataframe.loc[j+2,columns[i]] + dataframe.loc[j+3,columns[i]] + dataframe.loc[j-1,columns[i]] + dataframe.loc[j-2,columns[i]] + dataframe.loc[j-3,columns[i]]) / 7.


                            x1 = poids


                        poids += 0.003





                    alpha = (y1 - y0) / (x1 - x0)


                    beta = y1 - alpha*x1





                    if(f.endswith("AU45.csv")) :


                        alpha = 0





                    if(abs(alpha) < 0.3) :


                        alpha = 0





                    alphas[row_jacob - 1][i] = alpha


                    betas[row_jacob - 1][i] = beta


                


                row_jacob += 1


                





for i in range(len(alphas)) :


    for j in range(len(alphas[i])) :

        if(i == j) : 
            file.write(str(round(alphas[i][j] , 2)) + " ")
        else :
            file.write(str(0) + " ")


    file.write("\n")


for i in range(len(betas)) : 


    for j in range(len(betas[i])) :


        file.write(str(round(betas[i][j] , 2)) + " ")


    file.write("\n")








file.close()  