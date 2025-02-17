
# -*- coding: utf-8 -*-
# import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math as m
import pandas as pd
from matplotlib.ticker import MaxNLocator
from netCDF4 import Dataset
from statistics import mean
import time

def Read_inputs(i, p, filename):
    "Read inputs from file and store data to object"
    # Read data from compiled Excel file
    # Open Excel file
    df = pd.read_excel(filename,sheet_name='inputs')

    i.accu_mass_rep = df[df.columns[1]][1:].values.astype(float)
    i.AT_conc       = df[df.columns[2]][1:].values.astype(float)
    i.AT_height     = df[df.columns[3]][1:].values.astype(float)
    i.AT_Temp       = df[df.columns[4]][1:].values.astype(float)
    i.AT_Pres       = df[df.columns[5]][1:].values.astype(float)
    i.AT_O3_conc    = df[df.columns[6]][1:].values.astype(float)
    i.AT_D17O_O3b   = df[df.columns[7]][1:].values.astype(float)
    i.AT_f_exp      = df[df.columns[8]][1:].values.astype(float)
    i.FS_mass_rep   = df[df.columns[9]][1:].values.astype(float)
    i.FS_d15N       = df[df.columns[10]][1:].values.astype(float)
    i.FS_D17O       = df[df.columns[11]][1:].values.astype(float)
    i.FT_mass_rep   = df[df.columns[12]][1:].values.astype(float)
    i.FT_d15N       = df[df.columns[13]][1:].values.astype(float)
    i.FT_D17O       = df[df.columns[14]][1:].values.astype(float)
    i.O3colDU       = df[df.columns[15]][1:].values.astype(float)
    i.SZA_mean      = df[df.columns[16]][1:].values.astype(float)

    i.calc_AT_mass(p)


def accu_push(len):
    "Creates the array to use to push all layers one level down"
    v = np.ones(len)
    return np.diag(v,-1)

def indexMin(tt):
    "Returns the index of the smallest element in list"
    i, minv = 0, 1000000
    while i < len(tt):
        if tt[i] < minv :
            minv, imin = tt[i], i
        i = i + 1
    return imin


class Ozone_column(object):
    "Define the ozone column class"

    def __init__(self):
        self.DU200  = 0
        self.DU225  = 0
        self.DU250  = 0
        self.DU275  = 0
        self.DU300  = 0
        self.DU325  = 0
        self.DU350  = 0
        self.DU375  = 0
        self.DU400  = 0
        self.DU425  = 0
        self.DU450  = 0
        self.DU475  = 0
        self.DU500  = 0

class Options(object):
    "Define the options class"
    pass


class Params(object):
    "Define the parameters class"

    def __init__(self):
        self.Dt     = 0
        self.nb_yr  = 0
        self.nb_ts  = 0

    def calc_cp_dur(self):
        self.cp_dur = self.nb_yr * self.nb_ts + self.nb_end

    def calc_rho(self):
        self.SN_rho = self.SN_d * 1000

    def calc_M_NO3(self):
        self.M_NO3 = 1 * self.M_N + 3 * self.M_O

    def calc_layers_nb(self):
        self.SN_final_nb = int(self.SN_set_height/self.SN_set_thick)
        
    def cal_SN_SSA(self):
        self.SN_SSA = -308.2*np.log(self.SN_d)-205.96  # Domine et al., 2008


class State(object):
    "Define the state class"

    def __init__(self, p, i):
        load = open('initial_condition.txt','r')
        inco = np.genfromtxt(load,skip_header=1)
        load.close()
        index = [int(a/p.SN_set_thick*1e-2) for a in list(np.append(0,np.cumsum(inco[:,0])))] 
        number = index[-1] #total numbers of column
        self.SN_depth      = np.zeros((number))
        self.SN_thick      = np.zeros((number))
        self.SN_date       = np.zeros((number))
        self.SN_mass       = np.zeros((number)) #self.SN_thick*(p.M_N/p.M_NO3) * p.SN_rho * inco[:,3] * 1000 *1e-12 ,unit: kg(NO3-)
        self.SN_d15N       = np.zeros((number))
        self.SN_D17O       = np.zeros((number))
        self.SN_year       = np.zeros((number))
        for nbl in range(len(index)-1):
            index1 = index[nbl] #top layer
            index2 = index[nbl+1] #bottom layer
            for nb in range(index1,index2):
                self.SN_depth[nb] = (nb-0.5)*p.SN_set_thick
                self.SN_thick[nb] = p.SN_set_thick
                self.SN_date[nb] = -1
                self.SN_year[nb] = -1
                self.SN_mass[nb] = 0 # p.SN_set_thick*(p.M_N/p.M_NO3) * p.SN_rho * inco[nbl,1] * 1000 *1e-12
                self.SN_d15N[nb] = 0 # inco[nbl,2]
                self.SN_D17O[nb] = 0 # inco[nbl,3]        
        print('=== Initial conditions ===')
        print( 'Initial NO3- mass in snowpack ' + str(self.SN_mass.sum())+' kgN')
        print( 'Initial NO3- d15N in snowpack ' + str((self.SN_d15N*self.SN_mass).sum()/self.SN_mass.sum()) + ' permil')
        print( 'Initial NO3- D17O in snowpack ' + str((self.SN_D17O*self.SN_mass).sum()/self.SN_mass.sum()) + ' permil')
        
        self.AT_conc       = 2 # ngNO3.m-3 
        ABL_height         = 120. # m
        self.AT_mass       = self.AT_conc * (1e-12) * (p.M_N/p.M_NO3) * (ABL_height* p.cp_S)
        print('Initial NO3- mass in AT ' + str(self.AT_mass) + ' kgN')
        self.AT_d15N       = 0  
        self.AT_D17O       = 30.    
        self.FD_d15N       = self.AT_d15N
        self.FD_D17O       = self.AT_D17O
        print( 'Initial NO3- d15N in AT ' + str(self.AT_d15N) + ' permil')
        print('Initial NO3- D17O in AT ' + str(self.AT_D17O) + ' permil')
        # Set initial deposition flux to zero
        self.FD_mass = i.FT_mass[-1]*(1-i.FT_mass_rep[-1]) # kgN.m-2.s-1  
        print( '==========================')


class Inputs(object):
    "Define the inputs class"

    def calc_accu_mass(self, p):
        self.accu_mass = p.accu_mass_yr * self.accu_mass_rep

    def calc_FS_mass(self, p):
        # * cpS
        self.FS_mass = p.FS_FPI_frac * p.FPI_mass_yr * self.FS_mass_rep / p.Dt   # stratospheric denitrification flux at each time step (kgN.m^(-2).s^(-1))

    def calc_FT_mass(self, p):
        # * cpS
        self.FT_mass = p.FT_FPI_frac * p.FPI_mass_yr * self.FT_mass_rep / p.Dt   # stratospheric denitrification flux at each time step (kgN.m^(-2).s^(-1))

    def calc_AT_mass(self, p):
        self.AT_mass = (p.M_N / p.M_NO3) * self.AT_conc * p.cp_S * self.AT_height * pow(10, -12)     # nitrate mass in the atmospheric box (kgN.m^(-2))
        
    
class Outputs(object):
    "Define the output class"
    
    def __init__(self, p):
        self.FP_d15N_NO2     = np.zeros((1, p.cp_dur))
        self.FP_D17O_NO2     = np.zeros((1, p.cp_dur))
        self.FP_mass         = np.zeros((1, p.cp_dur))
        self.FP_d15N         = np.zeros((1, p.cp_dur))
        self.FP_D17O         = np.zeros((1, p.cp_dur))
        self.FD_mass         = np.zeros((1, p.cp_dur))
        self.FD_d15N         = np.zeros((1, p.cp_dur))
        self.FD_D17O         = np.zeros((1, p.cp_dur))
        self.AT_mass         = np.zeros((1, p.cp_dur))
        self.AT_conc         = np.zeros((1, p.cp_dur))
        self.AT_d15N         = np.zeros((1, p.cp_dur))
        self.AT_D17O         = np.zeros((1, p.cp_dur))
        self.AT_D17O_addO    = np.zeros((1, p.cp_dur))
        self.OH_frac         = np.zeros((1, p.cp_dur))
        self.D17O_NO2ini     = np.zeros((1, p.cp_dur))
        self.AT_D17O_addO    = np.zeros((1, p.cp_dur))
        self.AT_OH_conc      = np.zeros((1, p.cp_dur))
        self.AT_HO2_conc     = np.zeros((1, p.cp_dur))
        self.AT_CH3O2_conc   = np.zeros((1, p.cp_dur))
        self.AT_D17O_NO2_PSS = np.zeros((1, p.cp_dur))
        self.AT_alpha        = np.zeros((1, p.cp_dur))
        self.AT_D17O_OH      = np.zeros((1, p.cp_dur))
        self.AT_JNO2         = np.zeros((1, p.cp_dur))

        self.FDall_mass         = np.zeros((1, p.cp_dur))
        self.FDall_d15N         = np.zeros((1, p.cp_dur))
        self.FDall_D17O         = np.zeros((1, p.cp_dur))
        self.FEall_mass         = np.zeros((1, p.cp_dur))
        self.FEall_d15N         = np.zeros((1, p.cp_dur))
        self.FEall_D17O         = np.zeros((1, p.cp_dur))
        self.FPall_mass         = np.zeros((1, p.cp_dur))
        self.FPall_d15N         = np.zeros((1, p.cp_dur))
        self.FPall_D17O         = np.zeros((1, p.cp_dur))
        #add output items
        
        self.FD_massall         = np.zeros((1,p.cp_dur))
        self.FP_massall         = np.zeros((1,p.cp_dur))
        self.J_NO3              = np.zeros((1,p.cp_dur))  #mean rate at surface of reaction NO3- + hv -->  NO2 + OH
        self.eps_surface        = np.zeros((1,p.cp_dur))
        self.frem               = np.zeros((p.cp_dur,p.cp_dur))
        self.surface_conc       = np.zeros((1, p.cp_dur))
        self.surface_d15N       = np.zeros((1, p.cp_dur))
        self.surface_D17O       = np.zeros((1, p.cp_dur))

    def calc_mid_depth(self):
        "Calculates the mid depth in each layer"
        self.SN_depth_mid = self.SN_depth - 0.5 * self.SN_thick

    def calc_conc(self, p):
        "Calculates nitrate concentrations in snow"
        self.SN_conc = (p.M_NO3 / p.M_N) * self.SN_mass / (self.SN_thick * p.cp_S * p.SN_rho) * pow(10, 9)      # in ngNO3/gSNOW

    def calc_NO2_flux_molecule(self, p):
        # Converts the photolytically produced NO2 flux (originally in kgN.m^(-2).s^(-1)) in molecules.cm^(-2).s^(-1)
        self.FP_mass_mol = (self.FP_mass / p.M_N) * p.N_A * pow(10, 3)*1e-4

    def calc_NO2_flux_nmol(self, p):
        # Converts the photolytically produced NO2 flux (originally in kgN.m^(-2).s^(-1)) in nmol.m^(-2).h^(-1)
        self.FP_mass_amount = (self.FP_mass / p.M_N) * pow(10, 12) * 3600

    def calc_skinlayer(self, p):
        "Calculated skinlayer according to a given thickness"

        # Define outputs
        self.SL_thick       =np.zeros((1, p.nb_ts))
        self.SL_depth       = np.zeros((1, p.nb_ts))
        self.SL_mass        = np.zeros((1, p.nb_ts))
        self.SL_conc        = np.zeros((1, p.nb_ts))
        self.SL_d15N        = np.zeros((1, p.nb_ts))
        self.SL_D17O        = np.zeros((1, p.nb_ts))

        # Gets size of the object to use to calculate the skinlayer
        (rows, cols) = self.SN_depth[-1].shape
        # Prepare transfer array
        T = np.zeros((1, rows))

        # Loop
        for t in range(cols):
            # if not enough depth then we keep zeros in arrays
            if self.SN_depth[:, t].max() > p.SL_thick:
                T = np.zeros((1, rows))
                n = 0
                while self.SN_depth[n, t] < p.SL_thick:
                    T[0, n] = 1
                    n = n+1
                T[0, n] = 1 - (self.SN_depth[n, t] - p.SL_thick)/self.SN_thick[n, t]

                self.SL_thick[:, t] = (np.dot(T, self.SN_thick[:, t])).sum()
                self.SL_depth[:, t] = (np.dot(T, self.SN_thick[:, t])).sum()
                self.SL_mass[:, t]  = np.dot(T, self.SN_mass[:, t])
                self.SL_conc[:, t]  = (p.M_NO3 / p.M_N) * self.SL_mass[:, t]/(self.SL_thick[:, t] * p.cp_S * p.SN_rho) * pow(10, 9)  # in ngNO3/gSNOW
                self.SL_d15N[:, t]  = np.dot(T, self.SN_mass[:, t] * self.SN_d15N[:, t]) / self.SL_mass[:, t]
                self.SL_D17O[:, t]  = np.dot(T, self.SN_mass[:, t] * self.SN_D17O[:, t]) / self.SL_mass[:, t]


def resample_profile(vect_init_bot, vect_init_thi, vect_sub_bot):
    "Resample initial vector to subsampled vector and returns the conversion array"

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

    # Calculate thickness in vect_sub_bot
    vect_sub_thi = np.dot(np.eye(len(vect_sub_bot)) - np.diag(np.ones(len(vect_sub_bot)-1),-1), vect_sub_bot)

    if abs((vect_sub_thi - np.dot(TRANS, vect_init_thi)).sum()) > 1e-6:
        print( abs((vect_sub_thi - np.dot(TRANS, vect_init_thi)).sum()))
        print( 'Problem with profile resampling')

    return TRANS



def resample_after_deposition(vect_init, accu, thick):
    "Resample snowpack after deposition of a new layer of a given accumulation"

    dim = len(vect_init)

    # Define transformation and intermediate arrays
    TRANS = np.zeros((dim, dim))
    A     = np.zeros((dim, dim))    # Intermediate array

    # Write data in intermediate array
    b = (1 + int(accu/thick)) - accu/thick
    
    for i in range(dim-1):
        A[i, i]   = b
        A[i+1, i] = 1 - b

    # Write data in transformation array
    a = thick / accu
    TRANS[int(accu/thick):, 1:] = A[0:(dim-int(accu/thick)), 0:dim-1]
    TRANS[0:int(accu/thick), 0] = a
    TRANS[int(accu/thick), 0] = 1 - int(accu/thick)*a
    for i in range(dim):
        TRANS[dim-1, i] = 1 - TRANS[0:dim-1, i].sum()

    return TRANS


def DgeffHNO3(SSA, rhosnow, tau, Temp, Pressure):
	"""
        Provided by Mickael Kerbrat
	SSA: Snow Specific Surface Area in m2/kg	DC surf ~ 40-60m2/kg
	
	rhosnow : density of snow in kg/m3
	
	tau: tortuosity 	fresh snow ~ 0.7
				aged snow ~ 0.5
	
	Temp : Temperature in Kelvin
	
	Pressure : Pressure in mbar
	"""
	
	# Crowley et al. IUPAC, ACP, 2010
	KlinCHNO3	= 7.5e-7*np.exp(4585.0/Temp) # [m]
	
	# Durham et al. Atmospheric Environment, Volume 20, Issue 3, 1986, Pages 559-563 
	DgHNO30		= 0.118 # [cm2/s]
	
	# Tempearture Dependance Dg selon Massmann, Atmospheric Environment, 1998
	DgHNO3		= DgHNO30*(1e-2)**2*(1013.25/Pressure)*(Temp/298.15)**1.81 # [m2/s]
	
	#Conversion porosity
	rhoice		= 917.0 # [kg/m3]
	Phi		= 1. - rhosnow/rhoice # [-]
	
	# conversion aexp/vexp
	a_on_v_exp	= SSA*rhoice*(1.-Phi)/Phi # [m]
	
	
	DgeffHNO3	= 1./((tau**2/(DgHNO3*Phi))*(a_on_v_exp*KlinCHNO3+1.)) # [m2/s]
	
	return DgeffHNO3


def int_erfc_tau(alpha, tau, Pi):
    "Returns the integral of erfc function erfc(x/tau) inbetween 0 and alpha"

    return tau * ((alpha/tau) * m.erfc(alpha/tau) + (1 / (Pi**0.5)) * (1 - np.e**(-(alpha/tau)**2)))



def calc_diffusion(s, p, Diff, Dt, Dt_short, Dx):
    "Calculate nitrate diffusion in the snowpack"

    # Calculate number of layers which must undergo diffusion
    nb_layers = int(p.SN_set_height / p.photic_zone / p.SN_set_thick) + 1
    if nb_layers == 1001:
        nb_layers = 1000

    dim_small = nb_layers
    dim_wide  = len(s.SN_mass)

    # Fill these arrays with initial data
    d15N_mass = s.SN_mass[:] * s.SN_d15N[:]
    D17O_mass = s.SN_mass[:] * s.SN_D17O[:]

    # Create diffusion matrix
    # First create a small one (one col/line less)
    # Indeed, the last snow layer (the archived one, below 1m) must remain unchanged
    matrix_small = -2*np.eye(dim_small) + np.diag(np.ones(dim_small-1),-1) + np.diag(np.ones(dim_small-1),1)
    matrix_small[0, 0]         = -1
    matrix_small[dim_small-1, dim_small-1] = -1

    # Create a wider diffusion array
    matrix_wide = np.eye(dim_wide)
    matrix_wide[:dim_small, :dim_small] = np.eye(dim_small) + Diff*Dt_short/(Dx*Dx)*matrix_small

##    # Transform this array to remove rows and cols if photic zone is compressed
##    if p.photic_zone != 1:
##        # Calculate number of layers which must undergo diffusion
##        nb_layers = int(p.SN_set_height / p.photic_zone / p.SN_set_thick) + 1
##        if nb_layers == 1001:
##            nb_layers = 1000
##
##        # Remove rows and cols in the matrix_wide array which will be defined
##        # For this, calculate the corresponding indexes of the rows and cols
##        list_indexes = range(dim + 1)
##        for i in range(int((dim + 1)/2) - int((dim - nb_layers)/2), int((dim + 1)/2) - int((dim - nb_layers)/2) + (dim - nb_layers)):
##            list_indexes.remove(i)
##
##        temp = matrix_wide[list_indexes, :]
##        temp =        temp[:, list_indexes]
##        matrix_wide                             = np.eye(dim+1)
##        matrix_wide[:nb_layers+1, :nb_layers+1] = temp       

    time = 0
    for t in range(int(Dt/Dt_short)):
        time = time + Dt_short

        s.SN_mass[:] = np.dot(matrix_wide, s.SN_mass[:])
        d15N_mass[:] = np.dot(matrix_wide, d15N_mass[:])
        D17O_mass[:] = np.dot(matrix_wide, D17O_mass[:])
        
    Dt_last = Dt - time
    matrix_wide = np.eye(dim_wide)
    matrix_wide[:dim_small, :dim_small] = np.eye(dim_small) + Diff*Dt_last/(Dx*Dx)*matrix_small

    s.SN_mass[:] = np.dot(matrix_wide, s.SN_mass[:])
    d15N_mass[:] = np.dot(matrix_wide, d15N_mass[:])
    D17O_mass[:] = np.dot(matrix_wide, D17O_mass[:])

    # To finish, calculate the isotope arrays at end of diffusion routine
    s.SN_d15N[:] = d15N_mass[:] / s.SN_mass[:]
    s.SN_D17O[:] = D17O_mass[:] / s.SN_mass[:]


def plot_2D(p, o, what):

    range_year = range(p.cp_dur - p.nb_ts, p.cp_dur)

    if what == 'd15N':
        data = o.SPL_d15N
        vmin_val = int(data[:, range_year].min() / 10) * 10
        vmax_val = int(data[:, range_year].max() / 10) * 10 + 10
        inter    = 50
        legend='10$^3$ * $\delta$$^{15}$N(NO$_3$$^{-})$'
    elif what == 'D17O':
        data = o.SPL_D17O
        vmin_val = int(data[:, range_year].min())
        vmax_val = int(data[:, range_year].max()) + 1
        inter    = 5
        legend='10$^3$ * $\Delta$$^{17}$O(NO$_3$$^{-}$)'
    elif what == 'conc':
        data = o.SPL_conc
        vmin_val = int(data[:, range_year].min() / 100) * 100
        vmax_val = int(data[:, range_year].max() / 100) * 100 + 100
        inter    = 250
        legend='[NO$_3$$^{-}$], ng.g$^{-}$'

    Cmap=matplotlib.cm.jet
    Norm=matplotlib.colors.Normalize(vmin = vmin_val, vmax=vmax_val)
    scmap=matplotlib.cm.ScalarMappable(norm=Norm,cmap=Cmap)

    bounds = range(vmin_val, vmax_val , inter)

    fig = figure()
    # Define the size of the figure
    fig.set_size_inches(6.,3.)
    axs = axes()
##    axs.xaxis_date()
##    xlim(dat2[0],dat2[-1])
    # Define limits of the 2 axis
    xlim(p.cp_dur - p.nb_ts, p.cp_dur)
    ylim(1., 0.)
    # Define the labels of the 2 axis
    xlabel('Week from June 21st in the last computed year')
    ylabel('Snow depth (m)')

    # Double loop to plot each rectangle
    for i in range(p.nb_ts):
        # Below, we prefer to deal with the time step
        ts = p.cp_dur - p.nb_ts + i
        for j in range(len(o.SPL_depth[:, ts])):
            # the following coordinates represent the bottom left corner of each rectangle in the case of a normal y-axis!
            # Be aware that our y-axis is reverse!!
            coord = (ts-0.5, o.SPL_depth[j, ts] - o.SPL_thick[j, ts])
            # Height of the rectangle = thickness of the layer
            height = o.SPL_thick[j, ts]
            # Width of the rectangle = duration of one time step
            width = 1
            pat = patches.Rectangle(coord, width, height, edgecolor = 'none', fill = True, facecolor = scmap.to_rgba(data[j,ts]))
            axs.add_patch(pat)

    # Draw the legend (applies whatever mode : horizontal or vertical)
    # First argument : position of the box from left of the frame
    # Second argument : position of the box from bottom of the frame
    # Third argument : width of the box
    # Fourth argument : height of the box
    ax2=fig.add_axes([10, 0.1, 0.03, 0.8])
    # Color in the legend have the same properties as in the drawing
    cbar=matplotlib.colorbar.ColorbarBase(ax2, cmap=Cmap, norm=Norm, ticks=bounds, orientation='vertical')
    cbar.set_label(legend)

    # Export plot
    name2 = 'snowpack_%s' % (what)
    plt.savefig(name2, bbox_inches='tight', dpi=300)
    plt.close()


def plot_conc_t(o, t):
    "Plots the nitrate concentration vs snow depth at time step t"
    plot(o.SN_conc[0:len(o.SN_depth_mid), t], o.SN_depth_mid[0:len(o.SN_depth_mid), t], 'r-')                                     # plot definition
    name0 = 'plot_NO3_conc_week_' + str(t)                                      # plot's name
    title('Nitrate concentration profile, time step = ' + str(t))               # plot's title
    plt.xlabel('[NO$_3$$^{-}$], ng/g', fontsize=16)                         
    plt.ylabel('Snow depth, m', fontsize=16)
    plt.axis((0,1000,1,0))
    plt.savefig(name0)                                                          # saves figure in PNG
    plt.close()

def plot_d15N_t(o,p,t):
    "Plots the d15N in nitrate vs snow depth at time step t"
    plot(o.SN_d15N[0:p.SN_final, t], o.SN_depth_mid[0:p.SN_final, t], 'r-')                                     # plot definition
    name0 = 'plot_NO3_d15N_week_' + str(t)                                      # plot's name
##    title('Nitrate concentration profile, time step = ' + str(t))               # plot's title
    plt.xlabel('10$^3$ * $\delta$$^{15}$N(NO$_3$$^{-})$', fontsize=16)                         
    plt.ylabel('Snow depth, m', fontsize=16)
    plt.axis((0,400,1,0))
    plt.savefig(name0)                                                          # saves figure in PNG
    plt.close()

def plot_D17O_t(o, t):
    "Plots the D17O in nitrate vs snow depth at time step t"
    plot(o.SN_D17O[0:t-1, t], o.SN_depth_mid[0:t-1, t], 'r-')                                     # plot definition
    name0 = 'plot_NO3_D17O_week_' + str(t)                                      # plot's name
##    title('Nitrate concentration profile, time step = ' + str(t))               # plot's title
    plt.xlabel('10$^3$ * $\Delta$$^{17}$O(NO$_3$$^{-})$', fontsize=16)                         
    plt.ylabel('Snow depth, m', fontsize=16)
    plt.axis((0,42,1,0))
    plt.savefig(name0)                                                          # saves figure in PNG
    plt.close()

def calc_AIR_conc_in_moleccm3(P, T):
    'calculates AIR concentration in molecules.cm-3'
    'P in bar'
    'T in K'

    R   = 8.31              # (J.mol-1.kg-1)  ... pefect gas constant
    N_A = 6.022*10**23      # molecules.mol-1 ... Avogadro number

    return P*pow(10, 5) * N_A / (R*T*pow(10, 6))

def conv_P_in_hPa_to_bar(P):
    'converts P in hPa to bar'

    return P / 1000.

def O3_NO_195_308K(T):
    'rate constant in cm3 molecule-1 s-1'

    # O3 + NO (Atkinson et al., 2004, preferred value)
    # T range = 195-308 K
    # Atkinson et al., 2004, ACP Evaluated kinetic and photochemical data for atmospheric chemistry: Volume I - gas phase reactions of Ox, HOx, NOx and SOx species
    # Reaction # 54
    # k = 1.4*E−12 * exp(-1310/T) cm3 molecule−1 s−1   
    return 1.4*10**-12 * np.exp(-1310/T)

def BrO_NO_220_430K(T):
    'rate constant in cm3 molecule-1 s-1'

    # BrO + NO (Atkinson et al., 2007, preferred value)
    # T range = 220-430 K
    # Atkinson et al., 2007, ACP Evaluated kinetic and photochemical data for atmospheric chemistry: Volume III – gas phase reactions of inorganic halogens
    # Reaction # 76
    # k = 8.7*E−12 * exp(260/T) cm3 molecule−1 s−1   
    return 8.7*10**-12 * np.exp(260/T)

def HO2_NO_200_400K(T):
    'rate constant in cm3 molecule-1 s-1'

    # HO2 + NO (Atkinson et al., 2004, preferred value)
    # T range = 200-400 K
    # Atkinson et al., 2004, ACP Evaluated kinetic and photochemical data for atmospheric chemistry: Volume I - gas phase reactions of Ox, HOx, NOx and SOx species
    # Reaction # 45
    # k = 3.6*E−12 * exp(270/T) cm3 molecule−1 s−1   
    return 3.6*10**-12 * np.exp(270/T)

def CH3O2_NO_200_430K(T):
    'rate constant in cm3 molecule-1 s-1'

    # CH3O2 + NO (Atkinson et al., 2006, preferred value)
    # T range = 200-430 K
    # Atkinson et al., 2006, ACP Evaluated kinetic and photochemical data for atmospheric chemistry: Volume II – gas phase reactions of organic species
    # Reaction # 120
    # k = 2.3*E−12 * exp(360/T) cm3 molecule−1 s−1   
    return 2.3*10**-12 * np.exp(360/T)

def OH_H2O_iso_300_420K(T):
    'rate constant in cm3 molecule-1 s-1'

    # HO + H2O (isotopic exchange)
    # T range = 300-420 K
    # Dubey et al., 1997, J. Phys. Chem.
    # k = 2.3 +/- 1.0 * 10^-13 exp[-(2100 +/- 250)/T] cm3 molecule-1 s-1
    return 2.3*10**-13 * np.exp(-2100/T)

def OH_CO_200_300K(P, T):
    'rate constant in cm3 molecule-1 s-1'

    # HO + CO (Atkinson et al., 2006, preferred value)
    # T range = 200-300 K
    # P range = 0-1 bar of N2
    # Atkinson et al., 2006, ACP Evaluated kinetic and photochemical data for atmospheric chemistry: Volume II - gas phase reactions of organic species
    # k = 1.44*E-13 * (1 + P_N2/4.2 * E+19 molecule cm-3) cm3 molecule-1 s-1

    # Conversion of P (air pressure in bar) to P_N2 (partial pressure of N2 in molecules.cm-3)
    # the air contains 78% of N2
    P_N2 = 0.78 * calc_AIR_conc_in_moleccm3(P, T)
    
    return 1.44*10**-13 * (1 + P_N2/(4.2*10**19))

def OH_CH4_200_300K(T):
    'rate constant in cm3 molecule-1 s-1'

    # HO + CH4 (Atkinson et al., 2006, preferred value)
    # T range = 200-300 K
    # Atkinson et al., 2006, ACP Evaluated kinetic and photochemical data for atmospheric chemistry: Volume II - gas phase reactions of organic species
    # k = 1.85*E−12 * exp(-1690/T) cm3 molecule−1 s−1
    return 1.85*10**-12 * np.exp(-1690/T)

def calc_H2O_mixing_ratio_Bolton(P, T, RH):
    'calculates H2O mixing ratio, given RH, P and T'

    # From Bolton et al., 1980
    P_H2O = RH/100. * 6.112*np.exp(17.67*(T-273)/(T-29.5)) # in hPa
    P_H2O = conv_P_in_hPa_to_bar(P_H2O)                 # in bar
    return P_H2O/P

def calc_x_factor(P, T, RH, CO, CH4):
    """
    Returns the x factor necessary to calculate D17O(OH)
    No dimension
    """

    # CO  is expected in mol/mol
    # CH4 is expected in mol/mol
    # T   is in K
    # P   is in mbar

    # Import the necessary rate constants at the given temperature
    k_OH_H2O = OH_H2O_iso_300_420K(T)
    k_OH_CO  = OH_CO_200_300K(P, T)
    k_OH_CH4 = OH_CH4_200_300K(T)

    # Calculate H2O mixing ratio
    H2O = calc_H2O_mixing_ratio_Bolton(P, T, RH)
##    H20 = 3.5E-4  " Helene Barral, pers. comm. 2014)

    x = (k_OH_CO * CO + k_OH_CH4 * CH4) / (k_OH_CO * CO + k_OH_CH4 * CH4 + k_OH_H2O * H2O)
    return x

def my_resample(arr,accu,thickness,typea):
    #the depth array:first is the deposition layer, lower is the accumulated layers
    if accu <= 0: #remove the top layer
        num_lyr = int(-accu/thickness)+1 #number of layer
        return arr[num_lyr:]
    else:       
        if accu <= thickness:
            return arr
        else:
            num_lyr = len(arr)-1+int(accu/thickness)+1 #
            new_ar = np.zeros((num_lyr))
            new_ar[int(accu/thickness)+1:] = arr[1:]
            if typea == 'depth':
                for i in range(int(accu/thickness)):
                    new_ar[i] = (i+1)*thickness
                    new_ar[int(accu/thickness)] = accu
            elif typea == 'thickness':
                for i in range(int(accu/thickness)):
                    new_ar[i] = thickness
                    new_ar[int(accu/thickness)] = accu-int(accu/thickness)*thickness
            elif typea == 'd15N':
                new_ar[:int(accu/thickness)+1] = arr[0]
            elif typea == 'D17O':
                new_ar[:int(accu/thickness)+1] = arr[0]
            elif typea == 'date':
                new_ar[:int(accu/thickness)+1] = arr[0]
            elif typea == 'mass':
                mass_mean = arr[0]/accu #density
                for i in range(int(accu/thickness)):
                    new_ar[i] = thickness*mass_mean
                    new_ar[int(accu/thickness)] = (accu-int(accu/thickness)*thickness)*mass_mean
            elif typea == 'year':
                new_ar[:int(accu/thickness)+1] = arr[0]
            return new_ar

def plot_weekly(conc,thickness,d15N,D17O,week,SN_accu):
    #plot the weekly average data
    trans = np.zeros((52,3))
    for wek in range(52):
        concr = []
        thicknessr = []
        d15Nr = []
        D17Or = []
        for i in range(len(conc)):
            if week[i] == wek:
                concr.append(conc[i])
                thicknessr.append(thickness[i])
                d15Nr.append(d15N[i])
                D17Or.append(D17O[i])
        trans[wek,0] = (np.array(concr)*np.array(thicknessr)).sum()/np.array(thicknessr).sum()
        trans[wek,1] = (np.array(concr)*np.array(thicknessr)*np.array(d15Nr)).sum()/(np.array(thicknessr)*np.array(concr)).sum()
        trans[wek,2] = (np.array(concr)*np.array(thicknessr)*np.array(D17Or)).sum()/(np.array(thicknessr)*np.array(concr)).sum()
    #below we plot the weekly results
    #first squeeze the date in case of nan value
    conc = trans[:,0]
    d15N = trans[:,1]
    D17O = trans[:,2]
    conc = np.hstack((np.arange(52).reshape((52,1)),conc.reshape((52,1))))
    conc = conc[~np.isnan(conc).any(axis=1)]
    d15N = np.hstack((np.arange(52).reshape((52,1)),d15N.reshape((52,1))))
    d15N = d15N[~np.isnan(d15N).any(axis=1)]
    D17O = np.hstack((np.arange(52).reshape((52,1)),D17O.reshape((52,1))))
    D17O = D17O[~np.isnan(D17O).any(axis=1)]
    ###################plot     
    fig,(ax1,ax2) = plt.subplots(2,1)
    lns1 = ax1.plot(conc[:,0],conc[:,1],'bo-',markersize = 2, label = 'concentration')
    ax4 = ax1.twinx()
    lns2 = ax4.plot(range(52),SN_accu,linestyle = 'dashed',label = 'friction')
    ax4.set_xlim([0,52])
#    ax1.hlines(conc[:,1].mean(),0,51,linestyles=u'dashed')
#    ax1.text(10,ax1.get_ylim()[1]+3,'mean concentration = '+ format(conc[:,1].mean(),'0.2f')+' ng/g',fontsize=12)
    ax4.get_xaxis().set_visible(False)
    ax1.get_xaxis().set_visible(False)
    ax1.set_ylabel('ng/g', fontsize=12)
    ax4.set_ylabel('accumulation_friction', fontsize=12)
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax4.legend(lns, labs, loc=0)
    
    lns1 = ax2.plot(d15N[:,0],d15N[:,1],'ro-',markersize = 2,label = '$\delta$$^{15}$N')
    ax2.set_xlabel('week',fontsize=16)
    ax2.set_ylabel('10$^3$ * $\delta$$^{15}$N', fontsize=12)
    ax3 = ax2.twinx()
    lns2 = ax3.plot(D17O[:,0],D17O[:,1],'ko-',markersize = 2,label = '$\Delta$$^{17}$O' )
    ax3.set_xlabel('week',fontsize=16)
    ax3.set_ylabel('10$^3$ * $\Delta$$^{17}$O', fontsize=12)
    ax3.set_xlim([0,52])
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs, loc=3)
    ax2.yaxis.set_major_locator(MaxNLocator(prune='upper'))
    ax3.yaxis.set_major_locator(MaxNLocator(prune='upper'))
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.show()
    
def write_excel(o,title):
    '''write output in excel file'''
    # creat an array include thickness
    new_ma = np.zeros((len(o.mid_depth),5))
    new_ma[:,0] = o.mid_depth
    new_ma[:,1] = o.SN_conc
    new_ma[:,2] = o.SN_d15N
    new_ma[:,3] = o.SN_D17O
    new_ma[:,4] = o.SN_date
    data1 = pd.DataFrame(new_ma,columns = ['mid_depth','conc','d15N','D17O','week'])
    # another array include each flux
    me = np.zeros((len(o.FP_mass[0,:]),6))
    me[:,0] = o.FD_mass[0,:]
    me[:,1] = o.FD_d15N[0,:]
    me[:,2] = o.FD_D17O[0,:]
    me[:,3] = o.FP_mass[0,:]
    me[:,4] = o.FP_d15N[0,:]
    me[:,5] = o.FP_D17O[0,:]
    data2 = pd.DataFrame(me,columns = ['FD_mass','FD_d15N','FD_D17O','FP_mass','FP_d15N','FP_D17O'])
    # write data
    path = 'output/{}.xlsx'.format(title)
    with pd.ExcelWriter(path) as writer:  # doctest: +SKIP
        data1.to_excel(writer, sheet_name='Depth')
        data2.to_excel(writer, sheet_name='Flux')
    
def find_eps(data,frem,t):
    ary = np.ones((t+1))
    for k_data in range(t+1):
        res = []
        for layer in range(len(frem)):
            if data[layer] == k_data:
                res.append(frem[layer])
        ary[t-k_data] = mean(res)
    ary[np.isnan(ary)] = 1
    return ary

def frem_frac(frem):
    ts = len(frem[:,0])
    res = np.zeros((ts))
    for i in range(ts):
        me = 1
        ary = frem[:ts-i,i:]
        for j in range(ts-i):
            me = me*ary[j,j]
        res[i] = me
    return res

def create_ncfile(filename,data1,data2,paraml):
    nc = Dataset(filename, "w", format="NETCDF4")
    nc.description = 'TRANSITS model output'
    nc.history     = 'Created ' + time.ctime(time.time())
    
    snow_profile = nc.createGroup('outputs1')
    snow_profile.description = 'Group with outputs of snow depth profile'
    nb_layer, _ = data1.shape # number of layers
    snow_profile.createDimension('LAnumber', nb_layer)
    # write values for each variable
    ncSN_depth = snow_profile.createVariable('SN_depth', 'f', ('LAnumber',))
    ncSN_depth.long_name = 'Snow : middle depth of the layer'
    ncSN_depth.units     = 'meters from snowpack surface'
    ncSN_depth[:]        = data1[:,0]
    ncSN_conc = snow_profile.createVariable('SN_conc', 'f', ('LAnumber',))
    ncSN_conc.long_name = 'Snow : middle depth of the layer'
    ncSN_conc.units     = 'ng/g nitrate in snow'
    ncSN_conc[:]        = data1[:,1]
    ncSN_d15N = snow_profile.createVariable('SN_d15N', 'f', ('LAnumber',))
    ncSN_d15N.long_name = 'd15N of nitrate of the layer'
    ncSN_d15N.units     = 'per mil'
    ncSN_d15N[:]        = data1[:,2]
    ncSN_D17O = snow_profile.createVariable('SN_D17O', 'f', ('LAnumber',))
    ncSN_D17O.long_name = 'D17O of nitrate of the layer'
    ncSN_D17O.units     = 'per mil'
    ncSN_D17O[:]        = data1[:,3]
    ncSN_data = snow_profile.createVariable('SN_data', 'i', ('LAnumber',))
    ncSN_data.long_name = 'deposition time of the layer'
    ncSN_data.units     = 'week'
    ncSN_data[:]        = data1[:,4].astype('int')
    ncSN_thickness = snow_profile.createVariable('SN_thickness', 'f', ('LAnumber',))
    ncSN_thickness.long_name = 'thickness of the layer'
    ncSN_thickness.units     = 'meter'
    ncSN_thickness[:]        = data1[:,5]
    # weekly variables
    weekly_profile = nc.createGroup('outputs2')
    weekly_profile.description = 'Group with outputs of weekly flux and eps'
    _, nb_week = data2.shape
    weekly_profile.createDimension('LWnumber', nb_week)
    # write weekly variable results
    ncFP_mass = weekly_profile.createVariable('FP_mass', 'f', ('LWnumber',))
    ncFP_mass.long_name = 'photolysis flux of all weeks'
    ncFP_mass.units     = 'kgN.m-2.s-1'
    ncFP_mass[:] = data2[0,:]
    ncFP_d15N = weekly_profile.createVariable('FP_d15N', 'f', ('LWnumber',))
    ncFP_d15N.long_name = 'd15N of photolysis flux of all weeks'
    ncFP_d15N.units     = 'per mil'
    ncFP_d15N[:] = data2[1,:]
    ncFP_D17O = weekly_profile.createVariable('FP_D17O', 'f', ('LWnumber',))
    ncFP_D17O.long_name = 'D17O of photolysis flux of all weeks'
    ncFP_D17O.units     = 'per mil'
    ncFP_D17O[:] = data2[2,:]
    ncFD_mass = weekly_profile.createVariable('FD_mass', 'f', ('LWnumber',))
    ncFD_mass.long_name = 'deposition flux of all weeks'
    ncFD_mass.units     = 'kgN.m-2.s-1'
    ncFD_mass[:] = data2[3,:]
    ncFD_d15N = weekly_profile.createVariable('FD_d15N', 'f', ('LWnumber',))
    ncFD_d15N.long_name = 'd15N of deposition flux of all weeks'
    ncFD_d15N.units     = 'per mil'
    ncFD_d15N[:] = data2[4,:]
    ncFD_D17O = weekly_profile.createVariable('FD_D17O', 'f', ('LWnumber',))
    ncFD_D17O.long_name = 'D17O of deposition flux of all weeks'
    ncFD_D17O.units     = 'per mil'
    ncFD_D17O[:] = data2[5,:]
    # surface snow eps value
    nceps_surface = weekly_profile.createVariable('eps_surface', 'f', ('LWnumber',))
    nceps_surface.long_name = 'nitrate photolysis fractionation factor'
    nceps_surface.units     = 'per mil'
    nceps_surface[:] = data2[6,:]
    ncfrac_rem = weekly_profile.createVariable('frac_rem', 'f', ('LWnumber',))
    ncfrac_rem.long_name = 'remaing fraction of nitrate after photolysis'
    ncfrac_rem.units     = 'fraction'
    ncfrac_rem[:] = data2[7,:]
    ncsurface_conc = weekly_profile.createVariable('surface_conc', 'f', ('LWnumber',))
    ncsurface_conc.long_name = 'surface snow concentration'
    ncsurface_conc.units     = 'ng,g-1'
    ncsurface_conc[:] = data2[8,:]            
    ncsurface_d15N = weekly_profile.createVariable('surface_d15N', 'f', ('LWnumber',))
    ncsurface_d15N.long_name = 'd15N of surface snow'
    ncsurface_d15N.units     = 'per mil'
    ncsurface_d15N[:] = data2[9,:]
    ncsurface_D17O = weekly_profile.createVariable('surface_D17O', 'f', ('LWnumber',))
    ncsurface_D17O.long_name = 'D17O of surface snow'
    ncsurface_D17O.units     = 'per mil'
    ncsurface_D17O[:] = data2[10,:]                                    
    # major parameters used in model
    param = nc.createGroup('parameters')
    param.description = 'major parameters used in model'
    param.createDimension('Dimone',1)
    # write values
    ncfout = param.createVariable('fout','f',('Dimone',))
    ncfout.long_name = 'export fraction of photolysis NOx flux'
    ncfout.units     = 'fraction'
    ncfout[:] = paraml[0]
    ncfcage = param.createVariable('fcage','f',('Dimone',))
    ncfcage.long_name = 'cage effect of nitrate photolysis'
    ncfcage.units     = 'fraction'
    ncfcage[:] = paraml[1]
    ncphi = param.createVariable('phi','f',('Dimone',))
    ncphi.long_name = 'photolysis quantum yield'
    ncphi.units     = 'fraction'
    ncphi[:] = paraml[2]   
    nc.close()











