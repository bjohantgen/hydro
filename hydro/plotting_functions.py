#!/usr/bin/env python
#
#
#
# Author: Brayden JoHantgen
# Last Update: 7/12/2026

#############
# Importing
#############
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import math_functions as mfun
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

#############
# Functions
#############
def plotparams(ax, labelsize=10):
    '''
    Basic plot params
        :param ax: axes to modify
        :type ax: matplotlib axes object
        :returns: modified matplotlib axes object
    '''
    ax.minorticks_on()
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(direction='in', which='both', labelsize=labelsize)
    ax.tick_params('both', length=8, width=1.8, which='major')
    ax.tick_params('both', length=4, width=1, which='minor')
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)
    return ax

def line_num(num, directory = "time_step_files"):
    '''
    Input:
        num: the number of the txt file.
        directory: the location of the txt files.
    Output:
        line_count: the number of lines in the given txt file.
    Description:
        This function opens the given text file and reads the number
        of lines in the file. 
    '''
    file_path = directory+"/"+str(num)+".txt"
    with open(file_path, 'r') as f:
        line_count = len(f.readlines())
    return line_count

def read_txt_files(num, directory = "time_step_files"):
    '''
    Input:
    Output:
    Description:
    '''
    file_path = directory+"/"+str(num)+".txt"

    file_data = {}
    line_count = line_num(num, directory = directory)

    with open(file_path, 'r') as f:
        line_list = f.readlines()
        for i in range(line_count):
            line_dict_split = line_list[i].split(": ")
            float_list = line_dict_split[1].split(" ")
            float_data = []
            for j in range(len(float_list)-1):
                float_data.append(float(float_list[j]))
            line_dict = {line_dict_split[0]: float_data}
            file_data.update(line_dict)
    df = pd.DataFrame(list(file_data.items()), columns = ['Physics', 'Values'])
    return df

def one_dimension_list(od_list, xnum):
    '''
    Input:
    Output:
    Description:
    '''
    one_d = od_list[1:-1]
    return one_d

def two_dimension_list(td_list, xnum, ynum):
    '''
    Input:
    Output:
    Description:
    '''
    two_d = []
    for i in range(xnum):
        fill = []
        for j in range(ynum):
            fill.append(td_list[(i * ynum + j)])
        two_d.append(fill)
    transpose = [list(row) for row in zip(*two_d)]
    return transpose

def plot_frame_1d(ax, df, param, directory = "time_step_files"):
    '''
    Input:
    Output:
    Description:
    '''
    cells = int(df["Values"][0][0])
    rnames = np.linspace(0.0, 1.0, (cells+1))
    r = 0.5 * (rnames[:-1] + rnames[1:])

    frame, = ax.plot([], [], alpha=0.5, color='black')

    if param == "vx":
        fill_list = one_dimension_list(df["Values"][9], cells)
        
    elif param == "vy":
        fill_list = one_dimension_list(df["Values"][10], cells)

    elif param == "vz":
        fill_list = one_dimension_list(df["Values"][11], cells)

    elif param == "Bx":
        fill_list = one_dimension_list(df["Values"][12], cells)

    elif param == "By":
        fill_list = one_dimension_list(df["Values"][13], cells)

    elif param == "Bz":
        fill_list = one_dimension_list(df["Values"][14], cells)

    elif param == "pressure":
        fill_list = one_dimension_list(df["Values"][7], cells)

    else:
        fill_list = one_dimension_list(df["Values"][8], cells)

    frame.set_data(r, fill_list)
    return frame

def animate_1d(i, animated_plot, param, directory = "time_step_files"):
    '''
    '''
    frame_df = read_txt_files(i, directory = directory)

    plt.title('t='+str(round(frame_df["Values"][6][0], 4)), fontsize=20)

    cells = frame_df["Values"][0][0]
    rnames = np.linspace(0.0, 1.0, (int(cells)+1))
    r = 0.5 * (rnames[:-1] + rnames[1:])

    if param == "vx":
        filllist = one_dimension_list(frame_df["Values"][9], cells)

    elif param == "vy":
        filllist = one_dimension_list(frame_df["Values"][10], cells)

    elif param == "vz":
        filllist = one_dimension_list(frame_df["Values"][11], cells)

    elif param == "Bx":
        filllist = one_dimension_list(frame_df["Values"][12], cells)

    elif param == "By":
        filllist = one_dimension_list(frame_df["Values"][13], cells)

    elif param == "Bz":
        filllist = one_dimension_list(frame_df["Values"][14], cells)

    elif param == "pressure":
        filllist = one_dimension_list(frame_df["Values"][7], cells)

    else:
        filllist = one_dimension_list(frame_df["Values"][8], cells)

    animated_plot.set_data(r, filllist)

def plot_frame_2d(ax, df, param, directory = "time_step_files"):
    '''
    Input:
    Output:
    Description:
    '''
    xnum = int(df["Values"][0][0])
    ynum = int(df["Values"][1][0])

    x_left = df["Values"][2][0]
    x_right = df["Values"][3][0]

    y_left = df["Values"][4][0]
    y_right = df["Values"][5][0]

    if param == "vx":
        A = ax.imshow(np.zeros((xnum, ynum)), cmap='plasma', interpolation='nearest',
                    vmax=1.0, vmin=-1.0, extent = [x_left, x_right, y_left, y_right])
        fill_list = two_dimension_list(df["Values"][9], xnum, ynum)

    elif param == "vy":
        A = ax.imshow(np.zeros((xnum, ynum)), cmap='plasma', interpolation='nearest',
                    vmax=1.0, vmin=-1.0, extent = [x_left, x_right, y_left, y_right])
        fill_list = two_dimension_list(df["Values"][10], xnum, ynum)

    elif param == "vz":
        A = ax.imshow(np.zeros((xnum, ynum)), cmap='plasma', interpolation='nearest',
                    vmax=1.0, vmin=-1.0, extent = [x_left, x_right, y_left, y_right])
        fill_list = two_dimension_list(df["Values"][11], xnum, ynum)

    elif param == "Bx":
        A = ax.imshow(np.zeros((xnum, ynum)), cmap='plasma', interpolation='nearest',
                    vmax=3.0, vmin=-3.0, extent = [x_left, x_right, y_left, y_right])
        fill_list = two_dimension_list(df["Values"][12], xnum, ynum)

    elif param == "By":
        A = ax.imshow(np.zeros((xnum, ynum)), cmap='plasma', interpolation='nearest',
                    vmax=3.0, vmin=-3.0, extent = [x_left, x_right, y_left, y_right])
        fill_list = two_dimension_list(df["Values"][13], xnum, ynum)

    elif param == "Bz":
        A = ax.imshow(np.zeros((xnum, ynum)), cmap='plasma', interpolation='nearest',
                    vmax=3.0, vmin=-3.0, extent = [x_left, x_right, y_left, y_right])
        fill_list = two_dimension_list(df["Values"][14], xnum, ynum)

    elif param == "pressure":
        A = ax.imshow(np.zeros((xnum, ynum)), cmap='plasma', interpolation='nearest',
                    vmax=5.0, extent = [x_left, x_right, y_left, y_right])
        fill_list = two_dimension_list(df["Values"][7], xnum, ynum)

    else:
        A = ax.imshow(np.zeros((xnum, ynum)), cmap='plasma', interpolation='nearest',
                    vmax=3.0, extent = [x_left, x_right, y_left, y_right])
        fill_list = two_dimension_list(df["Values"][8], xnum, ynum)

    A.set_data(fill_list)
    return A

def animate_2d(i, A, param, directory = "time_step_files"):
    '''
    '''
    frame_df = read_txt_files(i, directory = directory)

    xnum = int(frame_df["Values"][0][0])
    ynum = int(frame_df["Values"][1][0])

    plt.title("t="+str(round(frame_df["Values"][6][0], 4)), fontsize=20)

    if param == "vx":
        fill_list = two_dimension_list(frame_df["Values"][9], xnum, ynum)

    elif param == "vy":
        fill_list = two_dimension_list(frame_df["Values"][10], xnum, ynum)

    elif param == "vz":
        fill_list = two_dimension_list(frame_df["Values"][11], xnum, ynum)

    elif param == "Bx":
        fill_list = two_dimension_list(frame_df["Values"][12], xnum, ynum)

    elif param == "By":
        fill_list = two_dimension_list(frame_df["Values"][13], xnum, ynum)

    elif param == "Bz":
        fill_list = two_dimension_list(frame_df["Values"][14], xnum, ynum)

    elif param == "pressure":
        fill_list = two_dimension_list(frame_df["Values"][7], xnum, ynum)

    else:
        fill_list = two_dimension_list(frame_df["Values"][8], xnum, ynum)
    
    A.set_data(fill_list)
