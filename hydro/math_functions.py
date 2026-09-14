#!/usr/bin/env python
#
# This python file is used to perform some of the math for
# the plotting of the RUST simulation data.
#
# Author: Brayden JoHantgen
# Last Update: 7/11/2026

#############
# Importing
#############
import numpy as np

#############
# Functions
#############
def list_avg(l):
    '''
    '''
    avg = sum(l) / len(l)
    return avg

def total_pressure(p, bx, by, bz):
    '''
    '''
    b_sq = bx**2 + by**2 + bz**2
    p_t = p + 0.5 * b_sq
    return p_t

def total_energy(p, rho, vx, vy, vz, bx, by, bz, gamma):
    '''
    '''
    b_sq = bx**2 + by**2 + bz**2
    v_sq = vx**2 + vy**2 + vz**2
    e = p / (gamma - 1) + 0.5 * rho * v_sq + 0.5 * b_sq
    return e

def sound_speed(p, rho, gamma):
    '''
    '''
    c = np.sqrt(gamma * p / rho)
    return c

#def fast_magsonic_speed

#def slow_magsonic_speed

#def alfven_speed

#def ke_x

#def ke_y

#def ke_z

#def mag_energy

#def temp

#def entropy

def mach_number(p, rho, vx, vy, vz, gamma):
    '''
    '''
    v = np.sqrt(vx**2 + vy**2 + vz**2)
    c = sound_speed(p, rho, gamma)
    n = v / c
    return n

#def emf



def list_avg(l):
    '''
    Input:
    Output:
    Description:
    '''
    avg = sum(l) / len(l)
    return avg

def avg_ke_y(rho_list, vy_list):
    '''
    Input:
    Output:
    Description:
    '''
    key_list = []
    for i in range(len(rho_list)):
        key_list.append(0.5 * rho_list[i] * (vy_list[i])**2)

    avg_key = list_avg(key_list)
    return avg_key

def avg_mag(bx_list, by_list, bz_list):
    '''
    Input:
    Output:
    Description:
    '''
    mag_list = []
    for i in range(len(bx_list)):
        val = bx_list[i]**2 + by_list[i]**2 + bz_list[i]**2
        mag_list.append(val)

    avg_m = list_avg(mag_list)
    return avg_m

def mag_press(bx_list, by_list):
    '''
    '''
    mp_list = []
    for i in range(len(bx_list)):
        val = bx_list[i]**2 + by_list[i]**2
        mp_list.append(val)
    return mp_list

#def time_evolution_mag(filenum):
#    '''
#    Input:
#    Output:
#    Description:
#    '''
#    time = []
#    t_evo = []

#    df0 = read_txt_files(0)
#    b0_sq = avg_mag(df0["Values"][8], df0["Values"][9], df0["Values"][10])

#    for i in range(filenum):
#        df = read_txt_files(i)
#        time.append(df["Values"][2][0])
#        b_sq = avg_mag(df["Values"][8], df["Values"][9], df["Values"][10])
#        fill = (b_sq - b0_sq) / (8 * np.pi)
#        t_evo.append(fill)

#    return time, t_evo

#def time_evolution_ke(filenum, param):
#    '''
#    Input:
#    Output:
#    Description:
#    '''
#    time = []
#    t_evo = []
#    for i in range(filenum):
#        df = read_txt_files(i)
#        time.append(df["Values"][2][0])

#        if param == "key":
#            fill = avg_ke_y(df["Values"][4], df["Values"][6])
        
#        t_evo.append(fill)
#
#    return time, t_evo