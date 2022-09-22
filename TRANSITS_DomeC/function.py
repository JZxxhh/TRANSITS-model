# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 01:28:48 2022

@author: lenovo
"""
# Embedded file name: function.py
import numpy as np
from pylab import *
import datetime as dt
import operator
import matplotlib
import matplotlib.cbook as cbook
import matplotlib.dates as dates
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math as m

def Read_inputs(i, p, file_name):
    """Read inputs from file and store data to object"""
    file = open(file_name, 'r')
    time_series = np.genfromtxt(file, skip_header=2, delimiter='')
    file.close()
    time_series = np.transpose(time_series)
    i.accu_mass_rep = time_series[1, :]
    i.AT_conc = time_series[2, :]
    i.AT_height = time_series[3, :]
    i.FS_mass_rep = time_series[4, :]
    i.FS_d15N = time_series[5, :]
    i.FS_D17O = time_series[6, :]
    i.FT_mass_rep = time_series[7, :]
    i.FT_d15N = time_series[8, :]
    i.FT_D17O = time_series[9, :]
    i.f_exp = time_series[10, :]
    i.O3colDU = time_series[11, :]
    i.alpha = time_series[12, :]
    i.D17O_O3 = time_series[13, :]
    i.x_fact = time_series[14, :]
    i.SZA_mean = time_series[15, :]
    i.air_temp = time_series[16, :]
    i.air_pres = time_series[17, :]
    i.calc_AT_mass(p)


def accu_push(len):
    """Creates the array to use to push all layers one level down"""
    v = np.ones(len)
    return np.diag(v, -1)


def indexMin(tt):
    """Returns the index of the smallest element in list"""
    i, minv = (0, 1000000)
    while i < len(tt):
        if tt[i] < minv:
            minv, imin = tt[i], i
        i = i + 1

    return imin


class Ozone_column(object):
    """Define the ozone column class"""

    def __init__(self):
        self.DU100 = 0
        self.DU125 = 0
        self.DU150 = 0
        self.DU175 = 0
        self.DU200 = 0
        self.DU225 = 0
        self.DU250 = 0
        self.DU275 = 0
        self.DU300 = 0
        self.DU325 = 0
        self.DU350 = 0
        self.DU375 = 0
        self.DU400 = 0
        self.DU425 = 0
        self.DU450 = 0
        self.DU475 = 0
        self.DU500 = 0
        self.DU750 = 0
        self.DU1000 = 0
        self.DU1500 = 0
        self.DU2000 = 0
        self.DU300 = 0


class Options(object):
    """Define the options class"""
    pass


class Params(object):
    """Define the parameters class"""

    def __init__(self):
        self.Dt = 0
        self.nb_yr = 0
        self.nb_ts = 0

    def calc_cp_dur(self):
        self.cp_dur = self.nb_yr * self.nb_ts

    def calc_rho(self):
        self.SN_rho = self.SN_d * 1000

    def calc_M_NO3(self):
        self.M_NO3 = 1 * self.M_N + 3 * self.M_O

    def calc_layers_nb(self):
        self.SN_final_nb = int(self.SN_set_height / self.SN_set_thick)


class State(object):
    """Define the state class"""

    def __init__(self, p):
        number = p.SN_final_nb + 1
        self.SN_depth = np.zeros((number, 1))
        self.SN_thick = np.zeros((number, 1))
        self.SN_date = np.zeros((number, 1))
        self.SN_mass = np.zeros((number, 1))
        self.SN_d15N = np.zeros((number, 1))
        self.SN_D17O = np.zeros((number, 1))
        self.SN_J14 = np.zeros((number, 1))
        self.SN_eps15 = np.zeros((number, 1))
        self.SN_massR = np.zeros((number, 1))
        self.SN_d15NR = np.zeros((number, 1))
        self.SN_D17OR = np.zeros((number, 1))
        self.SN_massE = np.zeros((number, 1))
        self.SN_d15NE = np.zeros((number, 1))
        self.SN_D17OE = np.zeros((number, 1))
        self.SN_massC = np.zeros((number, 1))
        self.SN_d15NC = np.zeros((number, 1))
        self.SN_D17OC = np.zeros((number, 1))
        self.SN_massD = np.zeros((number, 1))
        self.SN_fracD = np.zeros((number, 1))
        self.SN_d15ND = np.zeros((number, 1))
        self.SN_D17OD = np.zeros((number, 1))
        self.SN_JO3 = np.zeros((number, 1))
        self.SN_d15N_mass = np.zeros((number, 1))
        self.SN_D17O_mass = np.zeros((number, 1))
        self.SN_date_mass = np.zeros((number, 1))
        self.SN_JO3_mass = np.zeros((number, 1))
        self.FP_mass = 0.0
        self.FP_d15N = 0.0
        self.FP_D17O = 0.0
        self.FD_d15N = 0.0
        self.FD_D17O = 0.0
        self.AT_D17O_addO = 0.0
        self.AT_D17O_NO2 = 0.0
        for i in range(p.SN_final_nb + 1):
            self.SN_thick[i][0] = p.SN_set_height / p.SN_final_nb
            self.SN_depth[i][0] = (i + 1) * (p.SN_set_height / p.SN_final_nb)
            self.SN_date[i][0] = -i * (p.SN_set_height / p.SN_final_nb) / (p.accu_mass_yr / p.SN_rho)
            conc = 100.0
            self.SN_mass[i][0] = conc * (p.M_N / p.M_NO3) * 1e-12 * p.SN_rho * (p.SN_set_height / p.SN_final_nb * p.cp_S) * 1000
            self.SN_d15N[i][0] = 50.0
            self.SN_D17O[i][0] = 30.0

        self.AT_conc = 5.0
        ABL_height = 100.0
        self.AT_mass = self.AT_conc * 1e-12 * (p.M_N / p.M_NO3) * (ABL_height * p.cp_S)
        self.AT_d15N = 5.0
        self.AT_D17O = 30.0
        self.FD_d15N = self.AT_d15N
        self.FD_D17O = self.AT_D17O
        self.FD_mass = 0.0


class Inputs(object):
    """Define the inputs class"""

    def calc_accu_mass(self, p):
        self.accu_mass = p.accu_mass_yr * self.accu_mass_rep

    def calc_FS_mass(self, p):
        self.FS_mass = p.FS_FPI_frac * p.FPI_mass_yr * self.FS_mass_rep / p.Dt

    def calc_FT_mass(self, p):
        self.FT_mass = p.FT_FPI_frac * p.FPI_mass_yr * self.FT_mass_rep / p.Dt

    def calc_AT_mass(self, p):
        self.AT_mass = p.M_N / p.M_NO3 * self.AT_conc * p.cp_S * self.AT_height * pow(10, -12)


class Outputs(object):
    """Define the output class"""

    def __init__(self, p):
        self.SN_depth = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_thick = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_mass = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_date = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_mass_int_TS = np.zeros((1, p.cp_dur))
        self.SN_mass_int_YR = np.zeros((1, p.nb_yr))
        self.SN_JO3 = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_J14 = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_eps15 = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_d15N = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_D17O = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_FP_mass_contri = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_FP_d15N_contri = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_FP_D17O_contri = np.zeros((p.SN_final_nb, p.nb_ts))
        self.FP_d15N_NOx = np.zeros((1, p.nb_ts))
        self.FP_D17O_NOx = np.zeros((1, p.nb_ts))
        self.FP_mass = np.zeros((1, p.nb_ts))
        self.FP_d15N = np.zeros((1, p.nb_ts))
        self.FP_D17O = np.zeros((1, p.nb_ts))
        self.FA_mass = np.zeros((1, p.cp_dur))
        self.FA_conc = np.zeros((1, p.cp_dur))
        self.FA_thick = np.zeros((1, p.cp_dur))
        self.FA_date = np.zeros((1, p.cp_dur))
        self.FA_JO3 = np.zeros((1, p.cp_dur))
        self.FA_d15N = np.zeros((1, p.cp_dur))
        self.FA_D17O = np.zeros((1, p.cp_dur))
        self.AT_mass_int_TS = np.zeros((1, p.cp_dur))
        self.AT_mass_int_YR = np.zeros((1, p.nb_yr))
        self.FD_mass = np.zeros((1, p.nb_ts))
        self.FD_d15N = np.zeros((1, p.nb_ts))
        self.FD_D17O = np.zeros((1, p.nb_ts))
        self.FE_mass = np.zeros((1, p.nb_ts))
        self.FE_d15N = np.zeros((1, p.nb_ts))
        self.FE_D17O = np.zeros((1, p.nb_ts))
        self.AT_mass = np.zeros((1, p.nb_ts))
        self.AT_conc = np.zeros((1, p.nb_ts))
        self.AT_d15N = np.zeros((1, p.nb_ts))
        self.AT_D17O = np.zeros((1, p.nb_ts))
        self.AT_D17O_addO = np.zeros((1, p.nb_ts))
        self.AT_D17O_NO2_PSS = np.zeros((1, p.nb_ts))
        self.FDall_mass = np.zeros((1, p.cp_dur))
        self.FDall_d15N = np.zeros((1, p.cp_dur))
        self.FDall_D17O = np.zeros((1, p.cp_dur))
        self.FEall_mass = np.zeros((1, p.cp_dur))
        self.FEall_d15N = np.zeros((1, p.cp_dur))
        self.FEall_D17O = np.zeros((1, p.cp_dur))
        self.FPall_mass = np.zeros((1, p.cp_dur))
        self.FPall_d15N = np.zeros((1, p.cp_dur))
        self.FPall_D17O = np.zeros((1, p.cp_dur))
        self.frem = np.zeros((p.SN_final_nb + 1, p.cp_dur))
        self.depth_re = np.zeros((p.SN_final_nb + 1, p.cp_dur))
        self.thick_re = np.zeros((p.SN_final_nb + 1, p.cp_dur))

    def calc_mid_depth(self):
        """Calculates the mid depth in each layer"""
        self.SN_depth_mid = self.SN_depth - 0.5 * self.SN_thick

    def calc_conc(self, p):
        """Calculates nitrate concentrations in snow"""
        self.SN_conc = p.M_NO3 / p.M_N * self.SN_mass / (self.SN_thick * p.cp_S * p.SN_rho) * pow(10, 9)

    def calc_NOx_flux_molecule(self, p):
        self.FP_mass_mol = self.FP_mass / p.M_N * p.N_A * pow(10, 3)

    def calc_NOx_flux_nmol(self, p):
        self.FP_mass_amount = self.FP_mass / p.M_N * pow(10, 12) * 3600

    def calc_skinlayer(self, p):
        """Calculated skinlayer according to a given thickness"""
        self.SL_depth = np.zeros((1, p.nb_ts))
        self.SL_thick = np.zeros((1, p.nb_ts))
        self.SL_mass = np.zeros((1, p.nb_ts))
        self.SL_conc = np.zeros((1, p.nb_ts))
        self.SL_d15N = np.zeros((1, p.nb_ts))
        self.SL_D17O = np.zeros((1, p.nb_ts))
        rows, cols = self.SN_depth.shape
        T = np.zeros((1, rows))
        for t in range(cols):
            if self.SN_depth[:, t].max() > p.SL_thick:
                T = np.zeros((1, rows))
                n = 0
                while self.SN_depth[n, t] < p.SL_thick:
                    T[0, n] = 1
                    n = n + 1

                T[0, n] = 1 - (self.SN_depth[n, t] - p.SL_thick) / self.SN_thick[n, t]
                self.SL_thick[:, t] = np.dot(T, self.SN_thick[:, t]).sum()
                self.SL_depth[:, t] = np.dot(T, self.SN_thick[:, t]).sum()
                self.SL_mass[:, t] = np.dot(T, self.SN_mass[:, t])
                self.SL_conc[:, t] = p.M_NO3 / p.M_N * self.SL_mass[:, t] / (self.SL_thick[:, t] * p.cp_S * p.SN_rho) * pow(10, 9)
                self.SL_d15N[:, t] = np.dot(T, self.SN_mass[:, t] * self.SN_d15N[:, t]) / self.SL_mass[:, t]
                self.SL_D17O[:, t] = np.dot(T, self.SN_mass[:, t] * self.SN_D17O[:, t]) / self.SL_mass[:, t]


def resample_profile(vect_init_bot, vect_init_thi, vect_sub_bot):
    """Resample initial vector to subsampled vector and returns the conversion array"""
    TRANS = np.zeros((len(vect_sub_bot), len(vect_init_bot)))
    vect_init_top = vect_init_bot - vect_init_thi
    for line in range(len(vect_sub_bot)):
        col = 0
        while vect_init_bot[col] < vect_sub_bot[line]:
            TRANS[line, col] = 1
            col = col + 1

        TRANS[line, col] = (vect_sub_bot[line] - vect_init_top[col]) / (vect_init_bot[col] - vect_init_top[col])

    for line in range(1, len(vect_sub_bot)):
        for col in range(len(vect_init_bot)):
            sum_col = sum(TRANS[0:line, col])
            TRANS[line, col] = TRANS[line, col] - sum_col

    vect_sub_thi = np.dot(np.eye(len(vect_sub_bot)) - np.diag(np.ones(len(vect_sub_bot) - 1), -1), vect_sub_bot)
    if abs((vect_sub_thi - np.dot(TRANS, vect_init_thi)).sum()) > 1e-06:
        print(abs((vect_sub_thi - np.dot(TRANS, vect_init_thi)).sum()))
        print('Problem with profile resampling')
    return TRANS


def resample_after_deposition(vect_init, accu, thick):
    """Resample snowpack after deposition of a new layer of a given accumulation"""
    dim = len(vect_init)
    TRANS = np.zeros((dim, dim))
    A = np.zeros((dim, dim))
    b = 1 + int(accu / thick) - accu / thick
    for i in range(dim - 1):
        A[i, i] = b
        A[i + 1, i] = 1 - b

    a = thick / accu
    TRANS[int(accu / thick):, 1:] = A[0:dim - int(accu / thick), 0:dim - 1]
    TRANS[0:int(accu / thick), 0] = a
    TRANS[int(accu / thick), 0] = 1 - int(accu / thick) * a
    for i in range(dim):
        TRANS[dim - 1, i] = 1 - TRANS[0:dim - 1, i].sum()

    return TRANS


def DgeffHNO3(SSA, rhosnow, tau, Temp, Pressure):
    """
    Provided by Mickael Kerbrat
    SSA: Snow Specific Surface Area in m2/kg        DC surf ~ 40-60m2/kg
    
    rhosnow : density of snow in kg/m3
    
    tau: tortuosity         fresh snow ~ 0.7
                            aged snow ~ 0.5
    
    Temp : Temperature in Kelvin
    
    Pressure : Pressure in mbar
    """
    KlinCHNO3 = 7.5e-07 * exp(4585.0 / Temp)
    DgHNO30 = 0.118
    DgHNO3 = DgHNO30 * 0.0001 * (1013.25 / Pressure) * (Temp / 298.15) ** 1.81
    rhoice = 917.0
    Phi = 1.0 - rhosnow / rhoice
    a_on_v_exp = SSA * rhoice * (1.0 - Phi) / Phi
    DgeffHNO3 = 1.0 / (tau ** 2 / (DgHNO3 * Phi) * (a_on_v_exp * KlinCHNO3 + 1.0))
    return DgeffHNO3


def int_erfc_tau(alpha, tau, Pi):
    """Returns the integral of erfc function erfc(x/tau) inbetween 0 and alpha"""
    return tau * (alpha / tau * m.erfc(alpha / tau) + 1 / Pi ** 0.5 * (1 - e ** (-(alpha / tau) ** 2)))


def calc_diffusion(s, Diff, Dt, Dt_short, Dx):
    """Calculate nitrate diffusion in the snowpack"""
    dim = len(s.SN_mass) - 1
    d15N_mass = s.SN_mass[:] * s.SN_d15N[:]
    D17O_mass = s.SN_mass[:] * s.SN_D17O[:]
    matrix_small = -2 * np.eye(dim) + np.diag(np.ones(dim - 1), -1) + np.diag(np.ones(dim - 1), 1)
    matrix_small[(0, 0)] = -1
    matrix_small[dim - 1, dim - 1] = -1
    matrix_wide = np.eye(dim + 1)
    matrix_wide[:dim, :dim] = np.eye(dim) + Diff * Dt_short / (Dx * Dx) * matrix_small
    time = 0
    for t in range(int(Dt / Dt_short)):
        time = time + Dt_short
        s.SN_mass[:] = np.dot(matrix_wide, s.SN_mass[:])
        d15N_mass[:] = np.dot(matrix_wide, d15N_mass[:])
        D17O_mass[:] = np.dot(matrix_wide, D17O_mass[:])

    Dt_last = Dt - time
    matrix_wide = np.eye(dim + 1)
    matrix_wide[:dim, :dim] = np.eye(dim) + Diff * Dt_last / (Dx * Dx) * matrix_small
    s.SN_mass[:] = np.dot(matrix_wide, s.SN_mass[:])
    d15N_mass[:] = np.dot(matrix_wide, d15N_mass[:])
    D17O_mass[:] = np.dot(matrix_wide, D17O_mass[:])
    s.SN_d15N[:] = d15N_mass[:] / s.SN_mass[:]
    s.SN_D17O[:] = D17O_mass[:] / s.SN_mass[:]


def plot_2D(p, o, what):
    range_year = range(p.cp_dur - p.nb_ts, p.cp_dur)
    if what == 'd15N':
        data = o.SPL_d15N
        vmin_val = int(data[:, range_year].min() / 10) * 10
        vmax_val = int(data[:, range_year].max() / 10) * 10 + 10
        inter = 50
        legend = '10$^3$ * $\\delta$$^{15}$N(NO$_3$$^{-})$'
    elif what == 'D17O':
        data = o.SPL_D17O
        vmin_val = int(data[:, range_year].min())
        vmax_val = int(data[:, range_year].max()) + 1
        inter = 5
        legend = '10$^3$ * $\\Delta$$^{17}$O(NO$_3$$^{-}$)'
    elif what == 'conc':
        data = o.SPL_conc
        vmin_val = int(data[:, range_year].min() / 100) * 100
        vmax_val = int(data[:, range_year].max() / 100) * 100 + 100
        inter = 250
        legend = '[NO$_3$$^{-}$], ng.g$^{-}$'
    Cmap = matplotlib.cm.jet
    Norm = matplotlib.colors.Normalize(vmin=vmin_val, vmax=vmax_val)
    scmap = matplotlib.cm.ScalarMappable(norm=Norm, cmap=Cmap)
    bounds = range(vmin_val, vmax_val, inter)
    fig = figure()
    fig.set_size_inches(6.0, 3.0)
    axs = axes()
    xlim(p.cp_dur - p.nb_ts, p.cp_dur)
    ylim(1.0, 0.0)
    xlabel('Week from June 21st in the last computed year')
    ylabel('Snow depth (m)')
    for i in range(p.nb_ts):
        ts = p.cp_dur - p.nb_ts + i
        for j in range(len(o.SPL_depth[:, ts])):
            coord = (ts - 0.5, o.SPL_depth[j, ts] - o.SPL_thick[j, ts])
            height = o.SPL_thick[j, ts]
            width = 1
            pat = patches.Rectangle(coord, width, height, edgecolor='none', fill=True, facecolor=scmap.to_rgba(data[j, ts]))
            axs.add_patch(pat)

    ax2 = fig.add_axes([10,
     0.1,
     0.03,
     0.8])
    cbar = matplotlib.colorbar.ColorbarBase(ax2, cmap=Cmap, norm=Norm, ticks=bounds, orientation='vertical')
    cbar.set_label(legend)
    name2 = 'snowpack_%s' % what
    plt.savefig(name2, bbox_inches='tight', dpi=300)
    plt.close()


def plot_conc_t(o, t):
    """Plots the nitrate concentration vs snow depth at time step t"""
    plot(o.SN_conc[0:t - 1, t], o.SN_depth_mid[0:t - 1, t], 'r-')
    name0 = 'plot_NO3_conc_week_' + str(t)
    title('Nitrate concentration profile, time step = ' + str(t))
    plt.xlabel('[NO$_3$$^{-}$], ng/g', fontsize=16)
    plt.ylabel('Snow depth, m', fontsize=16)
    plt.axis((0, 1000, 1, 0))
    plt.savefig(name0)
    plt.close()


def plot_d15N_t(o, t):
    """Plots the d15N in nitrate vs snow depth at time step t"""
    plot(o.SN_d15N[0:t - 1, t], o.SN_depth_mid[0:t - 1, t], 'r-')
    name0 = 'plot_NO3_d15N_week_' + str(t)
    plt.xlabel('10$^3$ * $\\delta$$^{15}$N(NO$_3$$^{-})$', fontsize=16)
    plt.ylabel('Snow depth, m', fontsize=16)
    plt.axis((0, 400, 1, 0))
    plt.savefig(name0)
    plt.close()


def plot_D17O_t(o, t):
    """Plots the D17O in nitrate vs snow depth at time step t"""
    plot(o.SN_D17O[0:t - 1, t], o.SN_depth_mid[0:t - 1, t], 'r-')
    name0 = 'plot_NO3_D17O_week_' + str(t)
    plt.xlabel('10$^3$ * $\\Delta$$^{17}$O(NO$_3$$^{-})$', fontsize=16)
    plt.ylabel('Snow depth, m', fontsize=16)
    plt.axis((0, 42, 1, 0))
    plt.savefig(name0)
    plt.close()


def Return_list_of_runs(case):
    """Returns a list of sensitivity tests"""
    L = []
    if case == 'sensitivity':
        L.append(['_sensi_AT_eps15_dep/',
         'AT_eps15_dep',
         0,
         '',
         0,
         '',
         0])
        L.append(['_sensi_AT_eps15_dep/',
         'AT_eps15_dep',
         -10,
         '',
         0,
         '',
         0])
        L.append(['_sensi_AT_height/',
         'AT_height',
         25,
         '',
         0,
         '',
         0])
        L.append(['_sensi_AT_height/',
         'AT_height',
         100,
         '',
         0,
         '',
         0])
        L.append(['_sensi_Diff/',
         'Diff',
         1.33e-10,
         '',
         0,
         '',
         0])
        L.append(['_sensi_Diff/',
         'Diff',
         1.33e-12,
         '',
         0,
         '',
         0])
    elif case == 'test':
        L.append(['_sensi_d15N_FT/',
         'FT_d15N',
         9999,
         '',
         0,
         '',
         0])
        L.append(['_sensi_D17O_FT/',
         'FT_D17O',
         8888,
         '',
         0,
         '',
         0])
        L.append(['_sensi_D17O_FS/',
         'FS_D17O',
         7777,
         '',
         0,
         '',
         0])
    elif case == 'East_Antarctica':
        L.append(['_sensi_ACCU/',
         'accu_mass',
         300,
         'f_exp',
         0.2,
         '',
         0])
        L.append(['_sensi_ACCU/',
         'accu_mass',
         25,
         'f_exp',
         0.2,
         '',
         0])
        L.append(['_sensi_ACCU/',
         'accu_mass',
         30,
         'f_exp',
         0.2,
         '',
         0])
        L.append(['_sensi_ACCU/',
         'accu_mass',
         40,
         'f_exp',
         0.2,
         '',
         0])
        L.append(['_sensi_ACCU/',
         'accu_mass',
         50,
         'f_exp',
         0.2,
         '',
         0])
        L.append(['_sensi_ACCU/',
         'accu_mass',
         75,
         'f_exp',
         0.2,
         '',
         0])
        L.append(['_sensi_ACCU/',
         'accu_mass',
         100,
         'f_exp',
         0.2,
         '',
         0])
        L.append(['_sensi_ACCU/',
         'accu_mass',
         150,
         'f_exp',
         0.2,
         '',
         0])
        L.append(['_sensi_ACCU/',
         'accu_mass',
         200,
         'f_exp',
         0.2,
         '',
         0])
    return L