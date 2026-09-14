#!/usr/bin/env python
#
# This python file is used to animate the checkpoint file data from
# the RUST simulation of MHD.
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
import plotting_functions as pfun
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
from functools import partial

#####################
# Loading directory
#####################
directory = "time_step_files"
lst = os.listdir(directory)
filenum = len(lst) - 1

################
# Text Prompts
################
print("Please take the time to answer the following questions.")
print("What parameter would you like to plot? (density, pressure, vx, vy, vz, Bx, By, Bz)")
input_param = input()
print("Would you like an animation?")
print("Please answer with yes or no.")
input_anim = input()

if input_anim == "no":
    animation = False
    print("There are",filenum-1,"frames. Which frame number would you like to plot?")
    input_frame = input()
else: 
    animation = True


print("Would you like to plot any of the following?")
print("Magnetic Pressure (y/n)")
input_mag_press = input()
#print("Kinetic Energy (x direction): yes or no")
#input_kex = input()
#print("Kinetic Energy (y direction): yes or no")
#input_key = input()
#print("Kinetic Energy (z direction): yes or no")
#input_kez = input()
#print("Avg Mach Number: yes or no")
#input_mach = input()
#print("Magnetic Energy: yes or no")
#input_mag = input()
if input_mag_press == "y":
    print("There are",filenum-1,"frames. Which frame number would you like to plot?")
    input_extra_frame = input()

############
# Plotting
############
fig, ax = plt.subplots(1,1)
ax = pfun.plotparams(ax)

if animation == False:
    df = pfun.read_txt_files(int(input_frame))
    xnum = int(df["Values"][0][0])
    ynum = int(df["Values"][1][0])

    plt.title("t="+str(round(df["Values"][6][0], 4)), fontsize=20)

    if ynum == 0:
        pfun.plot_frame_1d(ax, df, input_param)

        if input_param == "vx":
            ax.set_ylim(-0.4, 0.8)
            ax.set_ylabel(r'$v_{x}$', fontsize=15)

        elif input_param == "vy":
            ax.set_ylim(-1.8, 0.2)
            ax.set_ylabel(r'$v_{y}$', fontsize=15)

        elif input_param == "vz":
            ax.set_ylim(-1.1, 1.1)
            ax.set_ylabel(r'$v_{z}$', fontsize=15)

        elif input_param == "Bx":
            ax.set_ylim(-1.25, 1.25)
            ax.set_ylabel(r'$B_{x}$', fontsize=15)

        elif input_param == "By":    
            ax.set_ylim(-1.25, 1.25)
            ax.set_ylabel(r'$B_{y}$', fontsize=15)
        
        elif input_param == "Bz":
            ax.set_ylim(-1.25, 1.25)
            ax.set_ylabel(r'$B_{z}$', fontsize=15)

        elif input_param == "pressure":
            ax.set_ylim(-0.05, 1.15)
            ax.set_ylabel(r'$P$', fontsize=15)

        else:
            ax.set_ylim(0.0, 1.05)
            ax.set_ylabel(r'$\rho$', fontsize=15)

        ax.set_xlim(-0.1, 1.1)
        ax.set_xlabel(r'$x$', fontsize=15)
        
    else:
        if input_param == "vx":
            cdata = pfun.plot_frame_2d(ax, df, input_param)
            cbar = fig.colorbar(cdata, ax = ax)
            cbar.set_label(r'$v_x$', size=15)

        elif input_param == "vy":
            cdata = pfun.plot_frame_2d(ax, df, input_param)
            cbar = fig.colorbar(cdata, ax = ax)
            cbar.set_label(r'$v_y$', size=15)

        elif input_param == "vz":
            cdata = pfun.plot_frame_2d(ax, df, input_param)
            cbar = fig.colorbar(cdata, ax = ax)
            cbar.set_label(r'$v_z$', size=15)

        elif input_param == "Bx":
            cdata = pfun.plot_frame_2d(ax, df, input_param)
            cbar = fig.colorbar(cdata, ax = ax)
            cbar.set_label(r'$B_x$', size=15)

        elif input_param == "By":
            cdata = pfun.plot_frame_2d(ax, df, input_param)
            cbar = fig.colorbar(cdata, ax = ax)
            cbar.set_label(r'$B_y$', size=15)

        elif input_param == "Bz":
            cdata = pfun.plot_frame_2d(ax, df, input_param)
            cbar = fig.colorbar(cdata, ax = ax)
            cbar.set_label(r'$B_z$', size=15)

        elif input_param == "pressure":
            cdata = pfun.plot_frame_2d(ax, df, input_param)
            cbar = fig.colorbar(cdata, ax = ax)
            cbar.set_label(r'$P$', size=15)

        else:
            cdata = pfun.plot_frame_2d(ax, df, input_param)
            cbar = fig.colorbar(cdata, ax = ax)
            cbar.set_label(r'$\rho$', size=15)

        ax.set_xlabel(r'$x$', fontsize=15)
        ax.set_ylabel(r'$y$', fontsize=15)
        cbar.ax.tick_params(labelsize=10)

else:
    df0 = pfun.read_txt_files(0)
    xnum = int(df0["Values"][0][0])
    ynum = int(df0["Values"][1][0])
    
    if ynum == 0:
        animated_plot, = ax.plot([], [], alpha=0.5, color='black')
        animation = FuncAnimation(fig=fig, 
                    func=partial(pfun.animate_1d, animated_plot=animated_plot, param=input_param), 
                    frames=(filenum-1),
                    interval=100,
                    repeat=True
                    )

        if input_param == "vx":
            ax.set_ylim(-0.4, 0.8)
            ax.set_ylabel(r'$v_{x}$', fontsize=15)

        elif input_param == "vy":
            ax.set_ylim(-1.8, 0.2)
            ax.set_ylabel(r'$v_{y}$', fontsize=15)

        elif input_param == "vz":
            ax.set_ylim(-1.1, 1.1)
            ax.set_ylabel(r'$v_{z}$', fontsize=15)

        elif input_param == "Bx":
            ax.set_ylim(-1.25, 1.25)
            ax.set_ylabel(r'$B_{x}$', fontsize=15)

        elif input_param == "By":
            ax.set_ylim(-1.25, 1.25)
            ax.set_ylabel(r'$B_{y}$', fontsize=15)

        elif input_param == "Bz":
            ax.set_ylim(-1.25, 1.25)
            ax.set_ylabel(r'$B_{z}$', fontsize=15)

        elif input_param == "pressure":
            ax.set_ylim(-0.05, 1.15)
            ax.set_ylabel(r'$P$', fontsize=15)

        else:
            ax.set_ylim(0.0, 1.05)
            ax.set_ylabel(r'$\rho$', fontsize=15)
    
        ax.set_xlim(-0.1, 1.1)
        ax.set_xlabel(r'$x$', fontsize=15)

    else: 
        x_left = df0["Values"][2][0]
        x_right = df0["Values"][3][0]

        y_left = df0["Values"][4][0]
        y_right = df0["Values"][5][0]

        if input_param == "vx":
            A = ax.imshow(np.zeros((int(df0["Values"][0][0]), int(df0["Values"][1][0]))), cmap='plasma', interpolation='nearest',
                    vmax=1.0, vmin=-1.0, extent = [-0.5, 0.5, -0.5, 0.5])
            animation = FuncAnimation(fig=fig, 
                    func=partial(pfun.animate_2d, A=A, param=input_param), 
                    frames=filenum,
                    interval=30,
                    repeat=True,
                    )
            cbar = fig.colorbar(A, ax = ax)
            cbar.set_label(r'$v_{x}$', size=15)

        elif input_param == "vy":
            A = ax.imshow(np.zeros((int(df0["Values"][0][0]), int(df0["Values"][1][0]))), cmap='plasma', interpolation='nearest',
                    vmax=1.0, vmin=-1.0, extent = [x_left, x_right, y_left, y_right])
            animation = FuncAnimation(fig=fig, 
                    func=partial(pfun.animate_2d, A=A, param=input_param), 
                    frames=filenum,
                    interval=30,
                    repeat=True,
                    )
            cbar = fig.colorbar(A, ax = ax)
            cbar.set_label(r'$v_{y}$', size=15)

        elif input_param == "vz":
            A = ax.imshow(np.zeros((int(df0["Values"][0][0]), int(df0["Values"][1][0]))), cmap='plasma', interpolation='nearest',
                    vmax=1.0, vmin=-1.0, extent = [x_left, x_right, y_left, y_right])
            animation = FuncAnimation(fig=fig, 
                    func=partial(pfun.animate_2d, A=A, param=input_param), 
                    frames=filenum,
                    interval=30,
                    repeat=True,
                    )
            cbar = fig.colorbar(A, ax = ax)
            cbar.set_label(r'$v_{z}$', size=15)

        elif input_param == "Bx":
            A = ax.imshow(np.zeros((int(df0["Values"][0][0]), int(df0["Values"][1][0]))), cmap='plasma', interpolation='nearest',
                    vmax=3.0, vmin=-3.0, extent = [x_left, x_right, y_left, y_right])
            animation = FuncAnimation(fig=fig, 
                    func=partial(pfun.animate_2d, A=A, param=input_param), 
                    frames=filenum,
                    interval=30,
                    repeat=True,
                    )
            cbar = fig.colorbar(A, ax = ax)
            cbar.set_label(r'$B_{x}$', size=15)

        elif input_param == "By":
            A = ax.imshow(np.zeros((int(df0["Values"][0][0]), int(df0["Values"][1][0]))), cmap='plasma', interpolation='nearest',
                    vmax=3.0, vmin=-3.0, extent = [x_left, x_right, y_left, y_right])
            animation = FuncAnimation(fig=fig, 
                    func=partial(pfun.animate_2d, A=A, param=input_param), 
                    frames=filenum,
                    interval=30,
                    repeat=True,
                    )
            cbar = fig.colorbar(A, ax = ax)
            cbar.set_label(r'$B_{y}$', size=15)

        elif input_param == "Bz":
            A = ax.imshow(np.zeros((int(df0["Values"][0][0]), int(df0["Values"][1][0]))), cmap='plasma', interpolation='nearest',
                    vmax=3.0, vmin=-3.0, extent = [x_left, x_right, y_left, y_right])
            animation = FuncAnimation(fig=fig, 
                    func=partial(pfun.animate_2d, A=A, param=input_param), 
                    frames=filenum,
                    interval=30,
                    repeat=True,
                    )
            cbar = fig.colorbar(A, ax = ax)
            cbar.set_label(r'$B_{z}$', size=15)

        elif input_param == "pressure":
            A = ax.imshow(np.zeros((int(df0["Values"][0][0]), int(df0["Values"][1][0]))), cmap='plasma', interpolation='nearest',
                    vmax=5.0, extent = [x_left, x_right, y_left, y_right])
            animation = FuncAnimation(fig=fig, 
                    func=partial(pfun.animate_2d, A=A, param=input_param), 
                    frames=filenum,
                    interval=30,
                    repeat=True,
                    )
            cbar = fig.colorbar(A, ax = ax)
            cbar.set_label(r'$P$', size=15)

        else:
            A = ax.imshow(np.zeros((int(df0["Values"][0][0]), int(df0["Values"][1][0]))), cmap='plasma', interpolation='nearest',
                    vmax=3.0, extent = [x_left, x_right, y_left, y_right])
            animation = FuncAnimation(fig=fig, 
                    func=partial(pfun.animate_2d, A=A, param=input_param), 
                    frames=filenum,
                    interval=30,
                    repeat=True,
                    )
            cbar = fig.colorbar(A, ax = ax)
            cbar.set_label(r'$\rho$', size=15)

        ax.set_xlabel(r'$x$', fontsize=15)
        ax.set_ylabel(r'$y$', fontsize=15)
        cbar.ax.tick_params(labelsize=10)

plt.show()

if input_mag_press == "y":
    fig, ax = plt.subplots(1,1)
    ax = pfun.plotparams(ax)
    df = pfun.read_txt_files(input_extra_frame)

    xnum = int(df["Values"][0][0])
    ynum = int(df["Values"][1][0])

    x_left = df["Values"][2][0]
    x_right = df["Values"][3][0]

    y_left = df["Values"][4][0]
    y_right = df["Values"][5][0]

    mp_list = mfun.mag_press(df["Values"][12], df["Values"][13])

    A = ax.imshow(np.zeros((xnum, ynum)), cmap='plasma', interpolation='nearest',
                 extent = [x_left, x_right, y_left, y_right], aspect='auto')
    fill_list = pfun.two_dimension_list(mp_list, xnum, ynum)

    A.set_data(fill_list)

    cbar = fig.colorbar(A, ax = ax)
    cbar.set_label(r'$P_M$', size=15)

    ax.set_xlabel(r'$x$', fontsize=15)
    ax.set_ylabel(r'$y$', fontsize=15)
    cbar.ax.tick_params(labelsize=10)

    plt.show()
    

#if input_key == "yes":
#    fig, ax = plt.subplots(1,1)
#    ax = plotparams(ax)
#    time, key_evo = time_evolution_ke(filenum-1, "key")
#    ax.plot(time, key_evo, alpha = 0.5, color = "black")
#    ax.set_yscale('log')
#    ax.set_xlabel(r'$t$', fontsize=15)
#    ax.set_ylabel(r'$KE_{y}$', fontsize=15)
#    plt.show()

#if input_mag == "yes":
#    fig, ax = plt.subplots(1,1)
#    ax = plotparams(ax)
#    time, mag_evo = time_evolution_mag(filenum-1)
#    ax.plot(time, mag_evo, alpha = 0.5, color = "black")
#    ax.set_yscale('log')
#    ax.set_xlabel(r'$t$', fontsize=15)
#    ax.set_ylabel(r'$E_{mag}$', fontsize=15)
#    plt.show()