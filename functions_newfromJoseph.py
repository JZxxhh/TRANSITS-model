# -*- coding: utf-8 -*-
# import modules
from __future__ import division                             # true division
import numpy as np
from StringIO import StringIO                           
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
import xlrd


def Read_inputs(i, p, filename):
    "Read inputs from file and store data to object"

    # Read data from compiled Excel file
    # Open Excel file
    wb = xlrd.open_workbook(filename)
    # Open sheet in Excel file
    ws = wb.sheet_by_name('inputs')
    # Count number of rows
    nb_rows = ws.nrows

    # Create arrays for all inputs
    i.accu_mass_rep = np.zeros((p.nb_ts))   # fraction              accumulated fractions repartition per time step
    i.AT_conc       = np.zeros((p.nb_ts))   # ngNO3.m-3             idealized atmospheric nitrate concentration
    i.AT_height     = np.zeros((p.nb_ts))   # meters                atmospheric boundary layer height
    i.AT_Temp       = np.zeros((p.nb_ts))   # Kelvin                Air temperature
    i.AT_Pres       = np.zeros((p.nb_ts))   # mbar                  Air pressure
    i.AT_O3_conc    = np.zeros((p.nb_ts))   # ppbv                  Ozone concentration in (lower) atmosphere 
    i.AT_D17O_O3b   = np.zeros((p.nb_ts))   # permil                D17O of local ozone
    i.AT_f_exp      = np.zeros((p.nb_ts))   # fraction              horizontally exported fraction
    i.FS_mass_rep   = np.zeros((p.nb_ts))   # fraction              stratospheric mass fractions repartition per year
    i.FS_d15N       = np.zeros((p.nb_ts))   # permil                stratospheric d15N in nitrate
    i.FS_D17O       = np.zeros((p.nb_ts))   # permil                stratospheric D17O in nitrate
    i.FT_mass_rep   = np.zeros((p.nb_ts))   # fraction              long distance transport  mass fractions repartition per year
    i.FT_d15N       = np.zeros((p.nb_ts))   # permil                long distance transport d15N in nitrate
    i.FT_D17O       = np.zeros((p.nb_ts))   # permil                long distance transport D17O in nitrate
    i.O3colDU       = np.zeros((p.nb_ts))   # Dobson units          mean ozone column
    i.SZA_mean      = np.zeros((p.nb_ts))   # degrees               mean Solar Zenith Angle in time step

    for j in range(2, nb_rows):
        i.accu_mass_rep[j-2] = ws.cell_value(j, 1)
        i.AT_conc[j-2]       = ws.cell_value(j, 2)
        i.AT_height[j-2]     = ws.cell_value(j, 3)
        i.AT_Temp[j-2]       = ws.cell_value(j, 4)
        i.AT_Pres[j-2]       = ws.cell_value(j, 5)
        i.AT_O3_conc[j-2]    = ws.cell_value(j, 6)
        i.AT_D17O_O3b[j-2]   = ws.cell_value(j, 7)
        i.AT_f_exp[j-2]      = ws.cell_value(j, 8)
        i.FS_mass_rep[j-2]   = ws.cell_value(j, 9)
        i.FS_d15N[j-2]       = ws.cell_value(j, 10)
        i.FS_D17O[j-2]       = ws.cell_value(j, 11)
        i.FT_mass_rep[j-2]   = ws.cell_value(j, 12)
        i.FT_d15N[j-2]       = ws.cell_value(j, 13)
        i.FT_D17O[j-2]       = ws.cell_value(j, 14)
        i.O3colDU[j-2]       = ws.cell_value(j, 15)
        i.SZA_mean[j-2]      = ws.cell_value(j, 16)

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
        self.DU100  = 0
        self.DU125  = 0
        self.DU150  = 0
        self.DU175  = 0
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
        self.DU750  = 0
        self.DU1000  = 0
        self.DU1500  = 0
        self.DU2000  = 0
        self.DU300  = 0

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
        self.cp_dur = self.nb_yr * self.nb_ts

    def calc_rho(self):
        self.SN_rho = self.SN_d * 1000

    def calc_M_NO3(self):
        self.M_NO3 = 1 * self.M_N + 3 * self.M_O

    def calc_layers_nb(self):

        self.SN_final_nb = int(self.SN_set_height/self.SN_set_thick)


class State(object):
    "Define the state class"

    def __init__(self, p):
        number = p.SN_final_nb+1
        self.SN_depth      = np.zeros((number, 1))
        self.SN_thick      = np.zeros((number, 1))
        self.SN_date       = np.zeros((number, 1))
        self.SN_mass       = np.zeros((number, 1))
        self.SN_conc       = np.zeros((number, 1))
        self.SN_d15N       = np.zeros((number, 1))
        self.SN_D17O       = np.zeros((number, 1))
        self.SN_J14        = np.zeros((number, 1))
        self.SN_eps15      = np.zeros((number, 1))
        self.SN_massR      = np.zeros((number, 1))
        self.SN_d15NR      = np.zeros((number, 1))
        self.SN_D17OR      = np.zeros((number, 1))
        self.SN_massE      = np.zeros((number, 1))
        self.SN_d15NE      = np.zeros((number, 1))
        self.SN_D17OE      = np.zeros((number, 1))
        self.SN_massC      = np.zeros((number, 1))
        self.SN_d15NC      = np.zeros((number, 1))
        self.SN_D17OC      = np.zeros((number, 1))
        self.SN_massD      = np.zeros((number, 1))
        self.SN_fracD      = np.zeros((number, 1))
        self.SN_d15ND      = np.zeros((number, 1))
        self.SN_D17OD      = np.zeros((number, 1))
        self.SN_JO3        = np.zeros((number, 1))
        self.SN_JNO2       = np.zeros((number, 1))
        self.SN_d15N_mass  = np.zeros((number, 1))
        self.SN_D17O_mass  = np.zeros((number, 1))
        self.SN_date_mass  = np.zeros((number, 1))
        self.SN_JO3_mass   = np.zeros((number, 1))
        self.SN_JNO2_mass  = np.zeros((number, 1))
        self.FP_mass       = 0.
        self.FP_d15N       = 0.
        self.FP_D17O       = 0.
        self.FD_d15N       = 0.
        self.FD_D17O       = 0.
        self.AT_D17O_addO  = 0.
        self.AT_D17O_NO2   = 0.
        self.AT_JNO2       = 0.


        # Read initialization data from compiled Excel file
        filename = 'initial_snow_conditions.xlsx'
        # Open Excel file
        wb = xlrd.open_workbook(filename)
        # Open sheet in Excel file
        ws = wb.sheet_by_name('pit11')              # EXAMPLE for pit11
        # Count number of rows
        nb_rows = ws.nrows

        # Create arrays for all inputs
        init_depth = np.zeros((p.SN_final_nb+1))
        init_conc  = np.zeros((p.SN_final_nb+1))
        init_d15N  = np.zeros((p.SN_final_nb+1))
        init_D17O  = np.zeros((p.SN_final_nb+1))

        # Read the data and place them in the arrays
        for j in range(1, nb_rows):
            init_depth[j-1]      = ws.cell_value(j, 0)
            init_conc[j-1]       = ws.cell_value(j, 1)
            init_d15N[j-1]       = ws.cell_value(j, 2)
            init_D17O[j-1]       = ws.cell_value(j, 3)

        # Now, initialize the snowpack
        for i in range(p.SN_final_nb+1):
            # Give each layer a lower depth and a thickness (constant)
            self.SN_thick[i][0]  = p.SN_set_height/p.SN_final_nb
            self.SN_depth[i][0]  = (i+1) * (p.SN_set_height/p.SN_final_nb)
            # Set a date (of deposition, in years since beginning of simulation) to each layer 
            self.SN_date[i][0]   = -i * (p.SN_set_height/p.SN_final_nb) / (p.accu_mass_yr/p.SN_rho)
            # Set the nitrate concentration in each layer
            self.SN_conc[i][0] = init_conc[i] # ng g-1
            # Set the nitrate mass in each layer
            self.SN_mass[i][0]   = self.SN_conc[i][0] * (p.M_N/p.M_NO3) * 1e-12 * p.SN_rho * ((p.SN_set_height/p.SN_final_nb) * p.cp_S) * 1000
            # Set d15N value in initial nitrate
            self.SN_d15N[i][0] = init_d15N[i] # permil
            # Set D17O value in initial nitrate
            self.SN_D17O[i][0] = init_D17O[i] # permil
        self.AT_conc       = 6293. # ngNO3.m-3
        ABL_height         = 60. # m
        self.AT_mass       = self.AT_conc * (1e-12) * (p.M_N/p.M_NO3) * (ABL_height* p.cp_S)
        self.AT_d15N       = 7.6 #(pit1=6.1, pit2=7.4, pit3=8.3, pit4=5.3, pit5=7.5, pit6=5.2, pit7=11.1, pit8=7.6, pit9=-5.52, pit10=7.6, pit11=7.6, pit12=7.4)
        self.AT_D17O       = 28.
        self.FD_d15N       = self.AT_d15N
        self.FD_D17O       = self.AT_D17O
        # Set initial deposition flux to zero
        self.FD_mass = 0. # kgN.m-2.s-1  



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
        self.SN_depth        = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_thick        = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_mass         = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_date         = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_mass_int_TS  = np.zeros((1, p.cp_dur))
        self.SN_mass_int_YR  = np.zeros((1, p.nb_yr))
        self.SN_JO3          = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_JNO2         = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_J14          = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_eps15        = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_d15N         = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_D17O         = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_FP_mass_contri = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_FP_d15N_contri = np.zeros((p.SN_final_nb, p.nb_ts))
        self.SN_FP_D17O_contri = np.zeros((p.SN_final_nb, p.nb_ts))
        self.FP_d15N_NO2     = np.zeros((1, p.nb_ts))
        self.FP_D17O_NO2     = np.zeros((1, p.nb_ts))
        self.FP_mass         = np.zeros((1, p.nb_ts))
        self.FP_d15N         = np.zeros((1, p.nb_ts))
        self.FP_D17O         = np.zeros((1, p.nb_ts))
        self.FA_mass         = np.zeros((1, p.cp_dur))
        self.FA_conc         = np.zeros((1, p.cp_dur))
        self.FA_thick        = np.zeros((1, p.cp_dur))
        self.FA_date         = np.zeros((1, p.cp_dur))
        self.FA_JO3          = np.zeros((1, p.cp_dur))
        self.FA_JNO2         = np.zeros((1, p.cp_dur))
        self.FA_d15N         = np.zeros((1, p.cp_dur))
        self.FA_D17O         = np.zeros((1, p.cp_dur))
        self.AT_mass_int_TS  = np.zeros((1, p.cp_dur))
        self.AT_mass_int_YR  = np.zeros((1, p.nb_yr))
        self.FD_mass         = np.zeros((1, p.nb_ts))
        self.FD_d15N         = np.zeros((1, p.nb_ts))
        self.FD_D17O         = np.zeros((1, p.nb_ts))
        self.FE_mass         = np.zeros((1, p.nb_ts))
        self.FE_d15N         = np.zeros((1, p.nb_ts))
        self.FE_D17O         = np.zeros((1, p.nb_ts))
        self.AT_mass         = np.zeros((1, p.nb_ts))
        self.AT_conc         = np.zeros((1, p.nb_ts))
        self.AT_d15N         = np.zeros((1, p.nb_ts))
        self.AT_D17O         = np.zeros((1, p.nb_ts))
        self.AT_D17O_addO    = np.zeros((1, p.nb_ts))
        self.AT_frac_OH_oxidation = np.zeros((1, p.nb_ts))
        self.AT_OH_conc      = np.zeros((1, p.nb_ts))
        self.AT_HO2_conc     = np.zeros((1, p.nb_ts))
        self.AT_CH3O2_conc   = np.zeros((1, p.nb_ts))
        self.AT_D17O_NO2_PSS = np.zeros((1, p.nb_ts))
        self.AT_alpha        = np.zeros((1, p.nb_ts))
        self.AT_D17O_OH      = np.zeros((1, p.nb_ts))
        self.AT_JNO2         = np.zeros((1, p.nb_ts))
        self.SN_mass_reservoir     = np.zeros((1, p.nb_ts))
        self.SN_D17O_reservoir     = np.zeros((1, p.nb_ts))
        self.SN_d15N_reservoir     = np.zeros((1, p.nb_ts))
        self.FDall_mass         = np.zeros((1, p.cp_dur))
        self.FDall_d15N         = np.zeros((1, p.cp_dur))
        self.FDall_D17O         = np.zeros((1, p.cp_dur))
        self.FEall_mass         = np.zeros((1, p.cp_dur))
        self.FEall_d15N         = np.zeros((1, p.cp_dur))
        self.FEall_D17O         = np.zeros((1, p.cp_dur))
        self.FPall_mass         = np.zeros((1, p.cp_dur))
        self.FPall_d15N         = np.zeros((1, p.cp_dur))
        self.FPall_D17O         = np.zeros((1, p.cp_dur))



    def calc_mid_depth(self):
        "Calculates the mid depth in each layer"
        self.SN_depth_mid = self.SN_depth - 0.5 * self.SN_thick

    def calc_conc(self, p):
        "Calculates nitrate concentrations in snow"
        self.SN_conc = (p.M_NO3 / p.M_N) * self.SN_mass / (self.SN_thick * p.cp_S * p.SN_rho) * pow(10, 9)      # in ngNO3/gSNOW

    def calc_NO2_flux_molecule(self, p):
        # Converts the photolytically produced NO2 flux (originally in kgN.m^(-2).s^(-1)) in molecules.m^(-2).s^(-1)
        self.FP_mass_mol = (self.FP_mass / p.M_N) * p.N_A * pow(10, 3)

    def calc_NO2_flux_nmol(self, p):
        # Converts the photolytically produced NO2 flux (originally in kgN.m^(-2).s^(-1)) in nmol.m^(-2).h^(-1)
        self.FP_mass_amount = (self.FP_mass / p.M_N) * pow(10, 12) * 3600

    def calc_skinlayer(self, p):
        "Calculated skinlayer according to a given thickness"

        # Define outputs
        self.SL_depth       = np.zeros((1, p.nb_ts))
        self.SL_thick       = np.zeros((1, p.nb_ts))
        self.SL_mass        = np.zeros((1, p.nb_ts))
        self.SL_conc        = np.zeros((1, p.nb_ts))
        self.SL_d15N        = np.zeros((1, p.nb_ts))
        self.SL_D17O        = np.zeros((1, p.nb_ts))

        # Gets size of the object to use to calculate the skinlayer
        (rows, cols) = self.SN_depth.shape
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
        print abs((vect_sub_thi - np.dot(TRANS, vect_init_thi)).sum())
        print 'Problem with profile resampling'

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
	KlinCHNO3	= 7.5e-7*exp(4585.0/Temp) # [m]
	
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

    return tau * ((alpha/tau) * m.erfc(alpha/tau) + (1 / (Pi**0.5)) * (1 - e**(-(alpha/tau)**2)))



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
    plot(o.SN_conc[0:t-1, t], o.SN_depth_mid[0:t-1, t], 'r-')                                     # plot definition
    name0 = 'plot_NO3_conc_week_' + str(t)                                      # plot's name
    title('Nitrate concentration profile, time step = ' + str(t))               # plot's title
    plt.xlabel('[NO$_3$$^{-}$], ng/g', fontsize=16)                         
    plt.ylabel('Snow depth, m', fontsize=16)
    plt.axis((0,1000,1,0))
    plt.savefig(name0)                                                          # saves figure in PNG
    plt.close()

def plot_d15N_t(o, t):
    "Plots the d15N in nitrate vs snow depth at time step t"
    plot(o.SN_d15N[0:t-1, t], o.SN_depth_mid[0:t-1, t], 'r-')                                     # plot definition
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

def Return_list_runs_sensitivity_tests(L, scale, folder_name):
    "Returns a list of sensitivity tests"

    # Normal run
    L.append([folder_name,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'FS_d15N',  119,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'FT_d15N',  100,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'AT_eps15_dep',  0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'FPI_scaled',  10,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'FS_FPI_frac', scale*0.5,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'f_cage',  scale*0.15,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'AT_f_exp',  scale*0.2,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'accu_mass',  scale*28,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'SN_rho',  scale*300,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'photic_zone',  scale*1,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'actinic_flux',  scale*1,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'Phi',  scale*0.026,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'Diff',  scale*(1.3E-11),
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'AT_height',  500,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'accu_repartition_time_series',  'W=2xS',
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'accu_repartition_time_series',  'S=2xW',
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'O3_col_time_series',  'flat_100DU',
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'O3_col_time_series',  'flat_300DU',
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'O3_col_time_series',  'flat_500DU',
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'O3_col_time_series',  'flat_300DU+hole_100DU',
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'AT_conc_time_series',  'real_idealized_x10',
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])


    # Change D17O only
    
    L.append([folder_name,
                   'FS_D17O',    0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])


    L.append([folder_name,
                   'FT_D17O',    0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'AT_D17O_O3b', 0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'AT_O3_conc_scale', 2,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'AT_HO2_conc_scale', 2,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'AT_CH3O2_conc_scale', 2,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'AT_BrO_conc',    2.5E-12 * 2,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    L.append([folder_name,
                   'AT_Temp_shift',  -10,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0,
                   '',           0
                   ])

    return L




def Return_list_runs_East_Antarctica(L, folder_name, list_A):
    "Returns a list of sensitivity tests"

    for A in list_A:
        if A < 30:
            nb_yr = 35
        else:
            nb_yr = 25

        L.append([folder_name,
                  'accu_mass',   A,
                  'nb_yr',   nb_yr,
                  '',            0,
                  '',            0,
                  '',            0,
                  '',            0
                  ])        

    return L



def Return_list_runs_sensitivity_for_curves(L, folder_name, list_Phi):
    "Returns a list of sensitivity tests"

    for val_Phi in list_Phi:


        L_runs.append([folder_name,
                       'Phi',            val_Phi,
                       'O3_col_time_series',  'flat_100DU',
                       '',            0,
                       '',            0,
                       '',            0,
                       '',            0
                       ])

        L_runs.append([folder_name,
                       'Phi',            val_Phi,
                       'O3_col_time_series',  'flat_300DU',
                       '',            0,
                       '',            0,
                       '',            0,
                       '',            0
                       ])

        L_runs.append([folder_name,
                       'Phi',            val_Phi,
                       'O3_col_time_series',  'flat_500DU',
                       '',            0,
                       '',            0,
                       '',            0,
                       '',            0
                       ])

        L_runs.append([folder_name,
                       'Phi',            val_Phi,
                       'O3_col_time_series',  'flat_300DU+hole_100DU',
                       '',            0,
                       '',            0,
                       '',            0,
                       '',            0
                       ])

        L_runs.append([folder_name,
                       'Phi',            val_Phi,
                       '',            0,
                       '',            0,
                       '',            0,
                       '',            0,
                       '',            0
                       ])

        L_runs.append([folder_name,
                       'Phi',            val_Phi,
                       'FT_d15N',   100,
                       '',            0,
                       '',            0,
                       '',            0,
                       '',            0
                       ])

        L_runs.append([folder_name,
                       'Phi',            val_Phi,
                       'FS_d15N',   119,
                       '',            0,
                       '',            0,
                       '',            0,
                       '',            0
                       ])



def Return_list_runs_sensi_D17O_FA(L, folder_name, list_Phi, list_O3col):
    "Returns a list of sensitivity tests"

    for val_Phi in list_Phi:

        for val_O3col in list_O3col:
            if val_O3col != 'flat_DC':
                
                # Normal run
                L.append([folder_name,
                               'Phi',         val_Phi,
                               'f_cage',         0.15,
                               '',                  0,
                               '',                  0,
                               '',                  0,
                               'O3_col',      val_O3col
                               ])

                # Normal run, no cage effect
                L.append([folder_name,
                               'Phi',         val_Phi,
                               'f_cage',            0,
                               '',                  0,
                               '',                  0,
                               '',                  0,
                               'O3_col',      val_O3col
                               ])

                # All 0 but D17O(FT) = 30 permil
                L.append([folder_name,
                               'Phi',         val_Phi,
                               'NO2_cycling',    9999,
                               'FT_D17O',          30,
                               'FS_D17O',           0,
                               '',                  0,
                               'O3_col',      val_O3col
                               ])

                # All 0 but D17O(FS) = 42 permil
                L.append([folder_name,
                               'Phi',         val_Phi,
                               'NO2_cycling',    9999,
                               'FT_D17O',           0,
                               'FS_D17O',          42,
                               '',                  0,
                               'O3_col',      val_O3col
                               ])

                # All 0 but D17O(NO2, PSS) = 35 permil
                L.append([folder_name,
                               'Phi',         val_Phi,
                               'NO2_cycling',     -35,
                               'FT_D17O',           0,
                               'FS_D17O',           0,
                               '',                  0,
                               'O3_col',      val_O3col
                               ])

            else:
                
                # Same for realistic DC O3col
                # Normal run
                L.append([folder_name,
                               'Phi',         val_Phi,
                               'f_cage',         0.15,
                               '',                  0,
                               '',                  0,
                               '',                  0,
                               '',                  0
                               ])

                # Normal run, no cage effect
                L.append([folder_name,
                               'Phi',         val_Phi,
                               'f_cage',            0,
                               '',                  0,
                               '',                  0,
                               '',                  0,
                               '',                  0
                               ])

                # All 0 but D17O(FT) = 30 permil
                L.append([folder_name,
                               'Phi',         val_Phi,
                               'NO2_cycling',    9999,
                               'FT_D17O',          30,
                               'FS_D17O',           0,
                               '',                  0,
                               '',                  0
                               ])

                # All 0 but D17O(FS) = 42 permil
                L.append([folder_name,
                               'Phi',         val_Phi,
                               'NO2_cycling',    9999,
                               'FT_D17O',           0,
                               'FS_D17O',          42,
                               '',                  0,
                               '',                  0
                               ])

                # All 0 but D17O(NO2, PSS) = 35 permil
                L.append([folder_name,
                               'Phi',         val_Phi,
                               'NO2_cycling',     -35,
                               'FT_D17O',           0,
                               'FS_D17O',           0,
                               '',                  0,
                               '',                  0
                               ])
    

    return L


def Return_list_Vostok_ice_cases_1_and_2(L, folder_name, list_Phi, list_O3col):
    "Returns a list of runs"

    for val_Phi in list_Phi:
        for val_O3col in list_O3col:

            L.append([folder_name,
                      'local_recycling',    0,
                      'O3_col',     val_O3col,
                      'Phi',          val_Phi,
                      '',                   0,
                      '',                   0,
                      '',                   0
                      ])

            L.append([folder_name,
                      'local_recycling',    1,
                      'O3_col',     val_O3col,
                      'Phi',          val_Phi,
                      '',                   0,
                      '',                   0,
                      '',                   0
                      ])

    return L


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
    return 1.4*10**-12 * exp(-1310/T)

def BrO_NO_220_430K(T):
    'rate constant in cm3 molecule-1 s-1'

    # BrO + NO (Atkinson et al., 2007, preferred value)
    # T range = 220-430 K
    # Atkinson et al., 2007, ACP Evaluated kinetic and photochemical data for atmospheric chemistry: Volume III – gas phase reactions of inorganic halogens
    # Reaction # 76
    # k = 8.7*E−12 * exp(260/T) cm3 molecule−1 s−1   
    return 8.7*10**-12 * exp(260/T)

def HO2_NO_200_400K(T):
    'rate constant in cm3 molecule-1 s-1'

    # HO2 + NO (Atkinson et al., 2004, preferred value)
    # T range = 200-400 K
    # Atkinson et al., 2004, ACP Evaluated kinetic and photochemical data for atmospheric chemistry: Volume I - gas phase reactions of Ox, HOx, NOx and SOx species
    # Reaction # 45
    # k = 3.6*E−12 * exp(270/T) cm3 molecule−1 s−1   
    return 3.6*10**-12 * exp(270/T)

def CH3O2_NO_200_430K(T):
    'rate constant in cm3 molecule-1 s-1'

    # CH3O2 + NO (Atkinson et al., 2006, preferred value)
    # T range = 200-430 K
    # Atkinson et al., 2006, ACP Evaluated kinetic and photochemical data for atmospheric chemistry: Volume II – gas phase reactions of organic species
    # Reaction # 120
    # k = 2.3*E−12 * exp(360/T) cm3 molecule−1 s−1   
    return 2.3*10**-12 * exp(360/T)

def OH_H2O_iso_300_420K(T):
    'rate constant in cm3 molecule-1 s-1'

    # HO + H2O (isotopic exchange)
    # T range = 300-420 K
    # Dubey et al., 1997, J. Phys. Chem.
    # k = 2.3 +/- 1.0 * 10^-13 exp[-(2100 +/- 250)/T] cm3 molecule-1 s-1
    return 2.3*10**-13 * exp(-2100/T)

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
    return 1.85*10**-12 * exp(-1690/T)

def calc_H2O_mixing_ratio_Bolton(P, T, RH):
    'calculates H2O mixing ratio, given RH, P and T'

    # From Bolton et al., 1980
    P_H2O = RH/100. * 6.112*exp(17.67*(T-273)/(T-29.5)) # in hPa
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
