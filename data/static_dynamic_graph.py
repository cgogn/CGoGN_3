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

for f in onlyfiles :
    if(not (f.endswith('_static.csv'))) : 
        name = basename(f)
        name = name[:-4]
        if not((exists(filepath + name + "_static.csv"))) : 
            if f.endswith('.csv') :
                print(f)
                name_au = basename(f)
                name_au = name_au[:-4] + "_r"
                dataframe = pd.read_csv(filepath + f)
                x = dataframe.index.to_numpy()
                y = dataframe[name_au].to_numpy()
                poids = 0.

                x0 = 0
                y0 = 0
                y1 = 0
                x1 = 0
                draw = True
                
                new_x = []
                for i in range(dataframe.index.size) :                
                    # if(dataframe[name_au][0] != 0):
                    #     dataframe.loc[i,name_au] = dataframe[name_au][i] - dataframe[name_au][0]
                    if(i == 0) :
                        dataframe.loc[i,name_au] = dataframe[name_au][i+2]
                        dataframe.loc[i+1,name_au] = dataframe[name_au][i+2]
                        x0 = poids
                        y0 = dataframe.loc[0,name_au]
                    if(dataframe[name_au][i] == 0) and (dataframe[name_au][0] == 0) and (poids < 1) :
                        x0 = poids
                    if(poids >= x0 + 1.5) and (draw): 
                        draw = False
                        y1 = (dataframe.loc[i,name_au] + dataframe.loc[i+1,name_au] + dataframe.loc[i+2,name_au] + dataframe.loc[i+3,name_au] + dataframe.loc[i-1,name_au] + dataframe.loc[i-2,name_au] + dataframe.loc[i-3,name_au]) / 7.
                        x1 = poids
                    new_x.append(poids)
                    poids += 0.003

                new_x = pd.DataFrame(new_x , columns=['weights'])
                new_x = new_x['weights'].to_numpy()

                x_smooth = np.linspace(new_x.min(), new_x.max(), 1500)
                spline = make_interp_spline(new_x, y, k=3)
                y_smooth = spline(x_smooth)

                plt.plot(x_smooth, y_smooth, label=f'{name_au} (smoothed)')
                if(not(draw)) : 
                    plt.ylim(-1, 5.5)
                    plt.axline((x0,y0),(x1,y1) , color="red" , linewidth=2 , label="Approximation of linear weight")
                plt.xlabel('Weight CGoGN')
                plt.ylabel('Intensity detected by OpenFace')
                plt.title(f'{name_au[:-2]} Intensity over Time')
                plt.legend()

                #plt.show()
                plt.savefig(filepath + name_au + '.png')
                plt.close()
        else :
            if f.endswith('.csv') :
                name_au = basename(f)
                name_au = name_au[:-4] + "_r"
                dataframe = pd.read_csv(filepath + f)
                dataframe_static = pd.read_csv(filepath + name + "_static.csv")

                x = dataframe.index.to_numpy()
                y = dataframe[name_au].to_numpy()

                x_static = dataframe_static.index.to_numpy()
                y_static = dataframe_static[name_au].to_numpy()

                poids = 0.
                poids_static = 0.

                x0 = 0
                y0 = 0
                y1 = 0
                x1 = 0

                x0_static = 0
                y0_static = 0
                y1_static = 0
                x1_static = 0
                draw = True
                draw_static = True
                
                new_x = []
                for i in range(dataframe.index.size) :                
                    # if(dataframe[name_au][0] != 0):
                    #     dataframe.loc[i,name_au] = dataframe[name_au][i] - dataframe[name_au][0]
                    if(i == 0) :
                        dataframe.loc[i,name_au] = dataframe[name_au][i+2]
                        dataframe.loc[i+1,name_au] = dataframe[name_au][i+2]
                        x0 = poids
                        y0 = dataframe.loc[0,name_au]
                    if(dataframe[name_au][i] == 0) and (dataframe[name_au][0] == 0) and (poids < 1) :
                        x0 = poids
                    if(poids >= x0 + 1.5) and (draw): 
                        draw = False
                        y1 = (dataframe.loc[i,name_au] + dataframe.loc[i+1,name_au] + dataframe.loc[i+2,name_au] + dataframe.loc[i+3,name_au] + dataframe.loc[i-1,name_au] + dataframe.loc[i-2,name_au] + dataframe.loc[i-3,name_au]) / 7.
                        x1 = poids
                    new_x.append(poids)
                    poids += 0.003

                new_x_static = []
                for i in range(dataframe_static.index.size) :                
                    # if(dataframe_static[name_au][0] != 0):
                    #     dataframe_static.loc[i,name_au] = dataframe_static[name_au][i] - dataframe_static[name_au][0]
                    if(i == 0) :
                        dataframe_static.loc[i,name_au] = dataframe_static[name_au][i+2]
                        dataframe_static.loc[i+1,name_au] = dataframe_static[name_au][i+2]
                        x0_static = poids_static
                        y0_static = dataframe_static.loc[0,name_au]
                    if(dataframe_static[name_au][i] == 0) and (dataframe_static[name_au][0] == 0) and (poids_static < 1) :
                        x0_static = poids_static
                    if(poids_static >= x0_static + 1.5) and (draw_static): 
                        draw_static = False
                        y1_static = (dataframe_static.loc[i,name_au] + dataframe_static.loc[i+1,name_au] + dataframe_static.loc[i+2,name_au] + dataframe_static.loc[i+3,name_au] + dataframe_static.loc[i-1,name_au] + dataframe_static.loc[i-2,name_au] + dataframe_static.loc[i-3,name_au]) / 7.
                        x1_static = poids_static
                    new_x_static.append(poids_static)
                    poids_static += 0.003

                new_x_static = pd.DataFrame(new_x_static , columns=['weights'])
                new_x_static = new_x_static['weights'].to_numpy()

                plt.plot(new_x, y, label=f'{name_au} (dynamic)')
                plt.plot(new_x_static, y_static, label=f'{name_au} (static)')
                if(not(draw) and not(draw_static)) : 
                    plt.ylim(-1, 5.5)
                    plt.axline((x0_static,y0_static),(x1_static,y1_static) , color="red" , linewidth=2 , label="Approximation of linear weight (static)")
                    plt.axline((x0,y0),(x1,y1) , color="blue" , linewidth=2 , label="Approximation of linear weight (dynamic)")
                plt.xlabel('Weight CGoGN')
                plt.ylabel('Intensity detected by OpenFace')
                plt.title(f'{name_au[:-2]} Intensity over Time')
                plt.legend()

                #plt.show()
                plt.savefig(filepath + name_au + '.png')
                plt.close()