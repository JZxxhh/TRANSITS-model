#################################################################################################
#################################################################################################
### This is the TRANSITS (TRansfert of Atmospheric Nitrate Stable Isotopes To the Snow) model ###
#################################################################################################
#################################################################################################

### IMPORT MODULES ###\
# We here import commonly used modules
import numpy as np                  # array module                        
from pylab import *
from time import *                  # module to measure time
# We here import home-made modules and all the functions contained in them
import function as f
import model as m
import os
import glob
import time
import warnings
warnings.filterwarnings("ignore")   # Ignore warnings like division by 0
import netCDF4


# Read the Input file
input_file_to_read = 'input_time_REALISTIC_simulation.txt'
print('The file \"' + input_file_to_read + '\" is used to generate the inputs')
# Model's version
version  = '4.2'

###############################################################################
##### BASIC MODEL PARAMETERS ##################################################
# PARAMETERS object : parameters are constant values and are not scenarios as are the inputs
### Creation of the object
param = f.Params()
param.cp_S =                    1.              # m^2 
# Below are constants
param.M_N =                     14.             # g.mol^(-1)         nitrogen molar mass
param.M_O =                     16.             # g.mol^(-1)         oxygen molar mass
param.calc_M_NO3()              # CALCULATED      g.mol^(-1)         nitrate molar mass
param.N_A =                     6.02214179E+23  # molecules.mol^(-1) Avogadro's number
param.Pi =                      3.1415926535898
###############################################################################
### OPTIONS ###################################################################
# OPTIONS object : options are integers which set the model to specific states
### Creation of the object
option = f.Options()
### The object is filled in with attributes
### The options are changed below (in loop)
###############################################################################
### INPUTS CREATION ###########################################################
### Inputs vary according to the scenarios choosen for the simulation
# Creation of the object
inp = f.Inputs()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Read and import SZA_rep file
load = open("inputs_SZA_rep\inputs_SZA_rep.txt", "r")
inp.SZA_rep = np.genfromtxt(load, skip_header=0, delimiter="")
inp.SZA_rep = inp.SZA_rep * 60  # converted from minutes to seconds

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Read and import J_O3 file
# READ for every ozone column, FLIP (light to right), TRANSPOSE and TRANSFER to arrays
file = open("inputs_TUV\inputs_JO3.txt","r")
importJO3 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
# Create object in the input object
inp.JO3 = f.Ozone_column()
inp.JO3.DU0025 = importJO3[1, :]
inp.JO3.DU0050 = importJO3[2, :]
inp.JO3.DU0075 = importJO3[3, :]
inp.JO3.DU0100 = importJO3[4, :]
inp.JO3.DU0125 = importJO3[5, :]
inp.JO3.DU0150 = importJO3[6, :]
inp.JO3.DU0175 = importJO3[7, :]
inp.JO3.DU0200 = importJO3[8, :]
inp.JO3.DU0225 = importJO3[9, :]
inp.JO3.DU0250 = importJO3[10, :]
inp.JO3.DU0275 = importJO3[11, :]
inp.JO3.DU0300 = importJO3[12, :]
inp.JO3.DU0325 = importJO3[13, :]
inp.JO3.DU0350 = importJO3[14, :]
inp.JO3.DU0375 = importJO3[15, :]
inp.JO3.DU0400 = importJO3[16, :]
inp.JO3.DU0425 = importJO3[17, :]
inp.JO3.DU0450 = importJO3[18, :]
inp.JO3.DU0475 = importJO3[19, :]
inp.JO3.DU0500 = importJO3[20, :]
inp.JO3.DU0750 = importJO3[21, :]
inp.JO3.DU1000 = importJO3[22, :]
inp.JO3.DU1500 = importJO3[23, :]
inp.JO3.DU2000 = importJO3[24, :]
inp.JO3.DU3000 = importJO3[25, :]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Read and import J files
inp.J14NO3 = f.Ozone_column()
inp.J15NO3 = f.Ozone_column()
# Read header (depths) only
load = open("inputs_TUV\inputs_J14NO3_0100DU.txt", "r")
read_for_header = np.genfromtxt(load, skip_header=0, delimiter="")
inp.J14NO3.header = read_for_header[0, :]
# READ for every ozone column, FLIP (light to right), TRANSPOSE and TRANSFER to arrays
file = open("inputs_TUV\inputs_J14NO3_0025DU.txt","r")
inp.J14NO3.DU0025 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0025DU.txt","r")
inp.J15NO3.DU0025 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0050DU.txt","r")
inp.J14NO3.DU0050 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0050DU.txt","r")
inp.J15NO3.DU0050 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0075DU.txt","r")
inp.J14NO3.DU0075 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0075DU.txt","r")
inp.J15NO3.DU0075 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0100DU.txt","r")
inp.J14NO3.DU0100 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0100DU.txt","r")
inp.J15NO3.DU0100 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0125DU.txt","r")
inp.J14NO3.DU0125 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0125DU.txt","r")
inp.J15NO3.DU0125 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0150DU.txt","r")
inp.J14NO3.DU0150 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0150DU.txt","r")
inp.J15NO3.DU0150 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0175DU.txt","r")
inp.J14NO3.DU0175 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0175DU.txt","r")
inp.J15NO3.DU0175 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0200DU.txt","r")
inp.J14NO3.DU0200 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0200DU.txt","r")
inp.J15NO3.DU0200 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0225DU.txt","r")
inp.J14NO3.DU0225 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0225DU.txt","r")
inp.J15NO3.DU0225 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0250DU.txt","r")
inp.J14NO3.DU0250 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0250DU.txt","r")
inp.J15NO3.DU0250 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0275DU.txt","r")
inp.J14NO3.DU0275 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0275DU.txt","r")
inp.J15NO3.DU0275 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0300DU.txt","r")
inp.J14NO3.DU0300 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0300DU.txt","r")
inp.J15NO3.DU0300 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0325DU.txt","r")
inp.J14NO3.DU0325 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0325DU.txt","r")
inp.J15NO3.DU0325 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0350DU.txt","r")
inp.J14NO3.DU0350 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0350DU.txt","r")
inp.J15NO3.DU0350 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0375DU.txt","r")
inp.J14NO3.DU0375 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0375DU.txt","r")
inp.J15NO3.DU0375 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0400DU.txt","r")
inp.J14NO3.DU0400 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0400DU.txt","r")
inp.J15NO3.DU0400 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0425DU.txt","r")
inp.J14NO3.DU0425 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0425DU.txt","r")
inp.J15NO3.DU0425 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0450DU.txt","r")
inp.J14NO3.DU0450 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0450DU.txt","r")
inp.J15NO3.DU0450 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0475DU.txt","r")
inp.J14NO3.DU0475 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0475DU.txt","r")
inp.J15NO3.DU0475 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0500DU.txt","r")
inp.J14NO3.DU0500 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0500DU.txt","r")
inp.J15NO3.DU0500 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_0750DU.txt","r")
inp.J14NO3.DU0750 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_0750DU.txt","r")
inp.J15NO3.DU0750 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_1000DU.txt","r")
inp.J14NO3.DU1000 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_1000DU.txt","r")
inp.J15NO3.DU1000 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_1500DU.txt","r")
inp.J14NO3.DU1500 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_1500DU.txt","r")
inp.J15NO3.DU1500 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_2000DU.txt","r")
inp.J14NO3.DU2000 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_2000DU.txt","r")
inp.J15NO3.DU2000 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J14NO3_3000DU.txt","r")
inp.J14NO3.DU3000 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV\inputs_J15NO3_3000DU.txt","r")
inp.J15NO3.DU3000 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
# In J14NO3 and J15NO3, rows are increasing depths and cols are increasing SZA (from 50 to 90)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Temporal series of inputs are found in the input_file_to_read
f.Read_inputs(inp, param, input_file_to_read)
#################################################################################################
#################################################################################################
######## LOOPS
### Sensitivity or simulation test loop
### 3 variables can vary independently
#>> Single run

#run01 = ['_no_sensi/',
#         'accu_mass',   75,
#         'nb_yr',            25,
#         '',            0]
         
run01 = ['_No_sensi/',
         '',            0,
         '',            0,
         '',            0]         

set_runs = [run01]

###>> East Antarctica runs
##set_runs = f.Return_list_of_runs('East_Antarctica')

###>> Sensitivity tests runs
#set_runs = f.Return_list_of_runs('sensitivity')

count_runs = 0
total_runs = len(set_runs)
for run in set_runs:
    count_runs = count_runs + 1
    if (run[1] == '') & (run[3] == '') & (run[5] == ''):
        print( 'Single run of the model : started on ' + time.ctime(time.time()))
    else:
        print( 'Run ' + str(count_runs) + ' out of ' + str(total_runs) + ' : started on ' + time.ctime(time.time()))
        print( '>>>>> Par1/' + run[1] + ' = ' + str(run[2]) + ', Par2/' + str(run[3]) + ' = ' + str(run[4]) + ', Par3/' + str(run[5]) + ' = ' + str(run[6]))

    # Read inputs again
    f.Read_inputs(inp, param, input_file_to_read)

    ### Set values to each parameter in the general case :
    param.nb_yr =                   25              # years                     number of computed years in the simulation
    param.nb_ts =                   52              # time steps                number of time steps in one year
    param.calc_cp_dur()             # CALCULATED        no unit                     number of time steps in simulation
    param.Dt =                      (364*24*60*60)/param.nb_ts      # seconds                   time step duration (default = 1 year / 52) 10107 is the number of minutes in such a time step
    param.cp_S =                    1.              # m^2                       computed column surface
    param.SN_set_height =           1               # meters                    height of the total active snowpack
    param.SN_set_thick =            0.001           # meters                    thickness of each layer
    param.SN_d =                    0.3            # adimensional              snow density
    param.calc_rho()                                # CALCULATED        kg.m^(-3)                   snow volumic mass
    param.SN_SSA =                  38.             # m2/kg                     Surface Specific Area of snow
    param.SN_tau =                  0.67            # no unit                   Tortuosity (here taken from Pinzer 2010, otherwise fresh = 0.7, old = 0.5)
    param.Phi =                     0.026           # molecules.photon^(-1)     photolytic quantum yield
    param.f_cage =                  0.1             # fraction                  cage effect parameter
    param.D17O_water =              0.              # permil                    D17O of water in snow
    param.JO3_threshold =           1.0E-09         # s-1                      JO3 threshold value. For this value : 100% oxidation by OH
    param.JO3_threshold =           param.JO3_threshold*param.Dt
    param.SL_thick =                0.005           # meters                    skinlayer thickness
    param.photic_zone =             1.              # adimensional              photic zone compression factor (1 is 10 cm e-folding as usual)
    param.actinic_flux =            1.              # adimensional              change in the incoming actinic flux (because of solar variability or orbital parameters. Default is 1)
    param.accu_mass_yr =            28              # kg.m^(-2).yr^(-1)         accumulated snow mass per year
    param.FPI_mass_yr =             8.2E-6               # kgN.m^(-2).yr^(-1)        TOTAL primary input mass deposited per year 
    param.FS_FPI_frac =             0.5               # adimensional              fraction of the primary input that deposits as stratospheric nitrate per year
    param.FT_FPI_frac =             0.5               # adimensional              fraction of the primary input that comes as long distance nitrate per year
    param.AT_eps15_dep =            +10.            # permil                    15N/14N isotope fractionation at deposition of nitrate from AT to snow

    option.NO2_oxidation         =           1               # adimensional              NO2 oxidation option : 1=local oxi, 0=no local oxi
    option.dep_diff_scenario     =           2               # Diffusion option : 0 = 100% deposition on top layer + no diffusion at all,
                                                             #             1 = 100% deposition on top layer + diffusion from this deposited nitrate only
                                                             #             2 = 100% deposition on top layer + diffusion in 1-m snowpack at all time step
    option.Diff                  =    1.33e-11               # Diffusion coefficient to be used in case 2 (m2/s)
    option.O3col_variations      =           0               # adimensional              variations in O3col inputs :
                                                             #             0 = no variations from year to year, read the time series in input file
                                                             #             1 = O3col is varied from year to year and the values are looked for in a specific input
    option.sastrugi_effect       =           0               # adimensional              sudden deposition of a sastrugi :
                                                             #             0 = no sastrugi
                                                             #             1 = a sastrugi deposits suddently at given time step and of a given thicknes
    option.sastrugi_thick        =           0.1             # meters
    option.sastrugi_TS           =           20              # time step number (0 is 21 June)  
                                                            

    ##### Transfer appropriately the values of the parameters to be changed
    pars_names  = [run[1], run[3], run[5]]
    pars_values = [run[2], run[4], run[6]]
    for i in range(3):
        # change parameters
        if pars_names[i] == 'nb_yr':
            param.nb_yr = pars_values[i]
            param.calc_cp_dur()
        elif pars_names[i] == 'SN_rho':
            param.SN_rho = pars_values[i]
        elif pars_names[i] == 'SN_SSA':
            param.SN_SSA = pars_values[i]
        elif pars_names[i] == 'SN_tau':
            param.SN_tau = pars_values[i]
        elif pars_names[i] == 'Phi':
            param.Phi = pars_values[i]
        elif pars_names[i] == 'f_cage':
            param.f_cage = pars_values[i]
        elif pars_names[i] == 'JO3_threshold':
            param.JO3_threshold = pars_values[i] * param.Dt
        elif pars_names[i] == 'photic_zone':
            param.photic_zone = pars_values[i]
        elif pars_names[i] == 'actinic_flux':
            param.actinic_flux = pars_values[i]
        elif pars_names[i] == 'accu_mass':
            param.accu_mass_yr = pars_values[i]
        elif pars_names[i] == 'FPI_mass_yr':
            param.FPI_mass_yr = pars_values[i] * param.FPI_mass_yr
        elif pars_names[i] == 'FS_FPI_frac':
            param.FS_FPI_frac = pars_values[i]
            param.FT_FPI_frac = 1 - pars_values[i]
        elif pars_names[i] == 'FT_FPI_frac':
            param.FT_FPI_frac = pars_values[i]
            param.FS_FPI_frac = 1 - pars_values[i]
        elif pars_names[i] == 'AT_eps15_dep':
            param.AT_eps15_dep = pars_values[i]
        # change inputs
        elif pars_names[i] == 'AT_conc':
            for j in range(len(inp.AT_conc)):
                inp.AT_conc[j] = pars_values[i]
        elif pars_names[i] == 'AT_height':
            for j in range(len(inp.AT_height)):
                inp.AT_height[j] = pars_values[i]
        elif pars_names[i] == 'FS_d15N':
            for j in range(len(inp.FS_d15N)):
                inp.FS_d15N[j] = pars_values[i]
        elif pars_names[i] == 'FS_D17O':
            for j in range(len(inp.FS_D17O)):
                inp.FS_D17O[j] = pars_values[i]
        elif pars_names[i] == 'FT_d15N':
            for j in range(len(inp.FT_d15N)):
                inp.FT_d15N[j] = pars_values[i]
        elif pars_names[i] == 'FT_D17O':
            for j in range(len(inp.FT_D17O)):
                inp.FT_D17O[j] = pars_values[i]
        elif pars_names[i] == 'f_exp':
            for j in range(len(inp.f_exp)):
                inp.f_exp[j] = pars_values[i]
        elif pars_names[i] == 'O3col_DU':
            for j in range(len(inp.O3colDU)):
                inp.O3colDU[j] = pars_values[i]
        elif pars_names[i] == 'alpha':
            for j in range(len(inp.alpha)):
                inp.alpha[j] = pars_values[i]
        elif pars_names[i] == 'D17O_O3':
            for j in range(len(inp.D17O_O3)):
                inp.D17O_O3[j] = pars_values[i]
        elif pars_names[i] == 'x_fact':
            for j in range(len(inp.x_fact)):
                inp.x_fact[j] = pars_values[i]
        # change options
        elif pars_names[i] == 'NO2_oxidation':
            option.NO2_oxidation = pars_values[i]
        elif pars_names[i] == 'dep_diff_scenario':
            option.dep_diff_scenario = pars_values[i]
        elif pars_names[i] == 'Diff':
            option.Diff = pars_values[i]
        elif pars_names[i] == 'sastrugi_effect':
            option.sastrugi_effect = pars_values[i]
        elif pars_names[i] == 'sastrugi_thick':
            option.sastrugi_thick = pars_values[i]
        elif pars_names[i] == 'sastrugi_TS':
            option.sastrugi_TS = pars_values[i]

        else:
            pass

    ### EXTRA PARAMETERS CALCULATION
    # We here calculate the number of snow layers in the active snowpack zone
    param.calc_layers_nb()

    ### Warn user if active snowpack is too small
    if param.photic_zone*1 > param.SN_set_height:
        print( "!!! Be careful! Active zone in the snowpack (%s meters) is too small for the defined photic zone (factor %s) !!!" % (param.SN_set_height, param.photic_zone))
        print( "!!! Computation stopped")
        break
    
    ### INPUTS PREPARATION
    # Some input time series are further calculated
    inp.calc_accu_mass(param)   # snow accumulated mass in each time step              (kg.m^(-2))
    inp.calc_FS_mass(param)     # stratospheric denitrification flux at each time step (kgN.m^(-2).s^(-1))
    inp.calc_FT_mass(param)     # long distance transportflux at each time step        (kgN.m^(-2).s^(-1))

    ### RESET STATE VARIABLES + INITIALIZE
    # Build a snowpack with a set depth and thickness for each layer
      # Creation of the object
    state = f.State(param)

    ### RESET OUTPUTS VARIABLES + INITIALIZE
      # Creation of the object
    out = f.Outputs(param)

    #########################################################################################
    ### MODEL'S RUN
    # Below, we run the model
    m.model(param, option, state, inp, out)
    #########################################################################################

    ### ADDITIONAL OUTPUTS CALCULATIONS
    # Calculates the mid depth of each layer
    out.calc_mid_depth()
    # Calculates nitrate concentrations in snow
    out.calc_conc(param)
    # Converts the photolytically produced NOx flux (originally in kgN.m^(-2).s^(-1)) in molecules.m^(-2).s^(-1)
    out.calc_NOx_flux_molecule(param)
    # Converts the photolytically produced NOx flux (originally in kgN.m^(-2).s^(-1)) in nmol.m^(-2).h^(-1)
    out.calc_NOx_flux_nmol(param)
    # Convert snow layer dates
    out.FA_date[0, :] = out.FA_date[0, :]

    ###### ADD TO EXPORTED DATA in Netcdf
    # No more resampling
    # Extracts skinlayer
    out.calc_skinlayer(param)

    #########################################################################################
    ### Write in NetCDF file ################################################################
    #########################################################################################
    ### Give the file an automatic name
    if count_runs < 10:
        h_nr = '0' + str(count_runs)
    else:
        h_nr = str(count_runs)

    name123 = ['', '', '']

    for i in [1, 3, 5]:
        if run[i] == '':
            name123[int((i-1)/2)] = ''
        else:
            name123[int((i-1)/2)] = run[i]+'='+str(run[i+1])

    if not os.path.exists(os.getcwd()+'\\_simulated_data\\' + run[0]):
        os.makedirs(os.getcwd()+'\\_simulated_data\\' + run[0])

    name = './_simulated_data/' + run[0] + '/TRANSITS_simul_' + name123[0] + '_' + name123[1] + '_' + name123[2] + '.nc'
    nc = netCDF4.Dataset(name, 'w')

    nc.description = 'TRANSITS model v' + version
    nc.history     = 'Created ' + time.ctime(time.time())

    ### CREATE GROUPS and describe them
    grpout = nc.createGroup('outputs')
    grpout.description = 'Group with all the outputs of the model'
    grpinp = nc.createGroup('inputs')
    grpinp.description = 'Group with all the inputs to the model'
    grppar = nc.createGroup('parameters')
    grppar.description = 'Group with all the parameters to the model'
    grpopt = nc.createGroup('options')
    grpopt.description = 'Group with all the options of the model'

    ### CREATE DIMENSION in each group
    # createDimension(name, size... 0 or None is unlimited... means in can be appended to)
    grpout.createDimension('LAnumber', param.SN_final_nb)   # Total number of layers in the snowpack
    grpout.createDimension('TSyearnr', param.nb_ts)         # Number of time steps per year
    grpout.createDimension('TStotal_number', param.cp_dur)  # Total number of time step
    
    grpinp.createDimension('TSyearnr', param.nb_ts)         # Number of time steps per year
    
    grppar.createDimension('DIMzero', 1)          # Dimension 0 for single value
    grpopt.createDimension('DIMzero', 1)          # Dimension 0 for single value  

    ###################### Write outputs ####################################################
    # createVariable (variable's name, type, dimensions, zlib=True means compression is ON)
    # If 1 dimension : do not forget the ","
    ncSN_depth = grpout.createVariable('SN_depth', 'f', ('LAnumber', 'TSyearnr'), zlib=True)
    ncSN_depth.long_name = 'Snow : lower depth of the layer'
    ncSN_depth.units     = 'meters from snowpack surface'
    ncSN_depth_mid = grpout.createVariable('SN_depth_mid', 'f', ('LAnumber', 'TSyearnr'), zlib=True)
    ncSN_depth_mid.long_name = 'Snow : mid depth of the layer'
    ncSN_depth_mid.units     = 'meters from snowpack surface'
    ncSN_thick = grpout.createVariable('SN_thick', 'f', ('LAnumber', 'TSyearnr'), zlib=True)
    ncSN_thick.long_name = 'Snow : thickness of the layer'
    ncSN_thick.units     = 'meters'
    ncSN_mass  = grpout.createVariable('SN_mass', 'f', ('LAnumber', 'TSyearnr'), zlib=True)
    ncSN_mass.long_name  = 'Snow : nitrate mass in the layer'
    ncSN_mass.units      = 'kgN'
    ncSN_date  = grpout.createVariable('SN_date', 'f', ('LAnumber', 'TSyearnr'), zlib=True)
    ncSN_date.long_name  = 'Snow : date of each layer in the layer'
    ncSN_date.units      = 'years from beginning of simulation'
    ncSN_conc  = grpout.createVariable('SN_conc', 'f', ('LAnumber', 'TSyearnr'), zlib=True)
    ncSN_conc.long_name  = 'Snow : nitrate conc in the layer'
    ncSN_conc.units      = 'ngNO3-/g'
    ncSN_d15N  = grpout.createVariable('SN_d15N', 'f', ('LAnumber', 'TSyearnr'), zlib=True)
    ncSN_d15N.long_name  = 'Snow : d15N in nitrate in the layer'
    ncSN_d15N.units      = 'permil'
    ncSN_D17O  = grpout.createVariable('SN_D17O', 'f', ('LAnumber', 'TSyearnr'), zlib=True)
    ncSN_D17O.long_name  = 'Snow : D17O in nitrate in the layer'
    ncSN_D17O.units      = 'permil'
    ncSN_eps15  = grpout.createVariable('SN_eps15', 'f', ('LAnumber', 'TSyearnr'), zlib=True)
    ncSN_eps15.long_name = 'Snow : eps15 photolytic fractionation constant in the layer'
    ncSN_eps15.units     = 'permil'
    ncSN_J14    = grpout.createVariable('SN_J14', 'f', ('LAnumber', 'TSyearnr'), zlib=True)
    ncSN_J14.long_name = 'Snow : J14 photolytic rate constant in the layer'
    ncSN_J14.units     = '/s'
    ncSN_JO3    = grpout.createVariable('SN_JO3', 'f', ('LAnumber', 'TSyearnr'), zlib=True)
    ncSN_JO3.long_name = 'Snow : JO3 photolytic rate constant at time of the deposition of the layer'
    ncSN_JO3.units     = '/s'
    ncSN_FP_mass_contri  = grpout.createVariable('SN_FP_mass_contri', 'f', ('LAnumber', 'TSyearnr'), zlib=True)
    ncSN_FP_mass_contri.long_name  = 'Snow : NOx flux emitted from each layer'
    ncSN_FP_mass_contri.units      = 'kgN/m2/s'
    ncSN_FP_d15N_contri  = grpout.createVariable('SN_FP_d15N_contri', 'f', ('LAnumber', 'TSyearnr'), zlib=True)
    ncSN_FP_d15N_contri.long_name  = 'Snow : d15N in NOx flux emitted from each layer'
    ncSN_FP_d15N_contri.units      = 'permil'
    ncSN_FP_D17O_contri  = grpout.createVariable('SN_FP_D17O_contri', 'f', ('LAnumber', 'TSyearnr'), zlib=True)
    ncSN_FP_D17O_contri.long_name  = 'Snow : D17O in NOx flux emitted from each layer'
    ncSN_FP_D17O_contri.units      = 'permil'
    ncSL_depth  = grpout.createVariable('SL_depth', 'f', ('TSyearnr', ))
    ncSL_depth.long_name  = 'Output : lower depth of the skinlayer throughout the year (in each time step)'
    ncSL_depth.units      = 'meters'
    ncSL_thick  = grpout.createVariable('SL_thick', 'f', ('TSyearnr', ))
    ncSL_thick.long_name  = 'Output : thickness of the skinlayer throughout the year (in each time step)'
    ncSL_thick.units      = 'meters'
    ncSL_mass  = grpout.createVariable('SL_mass', 'f', ('TSyearnr', ))
    ncSL_mass.long_name  = 'Output : nitrate mass in the skinlayer throughout the year (in each time step)'
    ncSL_mass.units      = 'kgN'
    ncSL_conc  = grpout.createVariable('SL_conc', 'f', ('TSyearnr', ))
    ncSL_conc.long_name  = 'Output : nitrate conc in the skinlayer throughout the year (in each time step)'
    ncSL_conc.units      = 'ngNO3-/g'
    ncSL_d15N  = grpout.createVariable('SL_d15N', 'f', ('TSyearnr', ))
    ncSL_d15N.long_name  = 'Output : d15N in nitrate in the skinlayer throughout the year (in each time step)'
    ncSL_d15N.units      = 'permil'
    ncSL_D17O  = grpout.createVariable('SL_D17O', 'f', ('TSyearnr', ))
    ncSL_D17O.long_name  = 'Output : D17O in nitrate in the skinlayer throughout the year (in each time step)'
    ncSL_D17O.units      = 'permil'
    ncFD_mass  = grpout.createVariable('FD_mass', 'f', ('TSyearnr', ))
    ncFD_mass.long_name  = 'Deposited flux : nitrate mass from the atmosphere'
    ncFD_mass.units      = 'kgN/m2/s'
    ncFD_d15N  = grpout.createVariable('FD_d15N', 'f', ('TSyearnr', ))
    ncFD_d15N.long_name  = 'Deposited flux : d15N in nitrate'
    ncFD_d15N.units      = 'permil'
    ncFD_D17O  = grpout.createVariable('FD_D17O', 'f', ('TSyearnr', ))
    ncFD_D17O.long_name  = 'Deposited flux : D17O in nitrate'
    ncFD_D17O.units      = 'permil'
    ncFE_mass  = grpout.createVariable('FE_mass', 'f', ('TSyearnr', ))
    ncFE_mass.long_name  = 'Exported flux : nitrate mass leaving the atmosphere'
    ncFE_mass.units      = 'kgN/m2/s'
    ncFE_d15N  = grpout.createVariable('FE_d15N', 'f', ('TSyearnr', ))
    ncFE_d15N.long_name  = 'Exported flux : d15N in nitrate'
    ncFE_d15N.units      = 'permil'
    ncFE_D17O  = grpout.createVariable('FE_D17O', 'f', ('TSyearnr', ))
    ncFE_D17O.long_name  = 'Exported flux : D17O in nitrate'
    ncFE_D17O.units      = 'permil'
    ncFP_mass  = grpout.createVariable('FP_mass', 'f', ('TSyearnr', ))
    ncFP_mass.long_name  = 'Photolyzed flux : nitrate mass transported off the snowpack'
    ncFP_mass.units      = 'kgN/m2/s'
    ncFP_d15N  = grpout.createVariable('FP_d15N', 'f', ('TSyearnr', ))
    ncFP_d15N.long_name  = 'Photolyzed flux : d15N in nitrate'
    ncFP_d15N.units      = 'permil'
    ncFP_D17O  = grpout.createVariable('FP_D17O', 'f', ('TSyearnr', ))
    ncFP_D17O.long_name  = 'Photolyzed flux : D17O in nitrate after local oxidation'
    ncFP_D17O.units      = 'permil'
    ncFP_d15N_NOx  = grpout.createVariable('FP_d15N_NOx', 'f', ('TSyearnr', ))
    ncFP_d15N_NOx.long_name  = 'Photolyzed flux : d15N in photolyzed flux prior to local oxidation'
    ncFP_d15N_NOx.units      = 'permil'
    ncFP_D17O_NOx  = grpout.createVariable('FP_D17O_NOx', 'f', ('TSyearnr', ))
    ncFP_D17O_NOx.long_name  = 'Photolyzed flux : D17O in photolyzed flux prior to local oxidation'
    ncFP_D17O_NOx.units      = 'permil'
    ncAT_mass  = grpout.createVariable('AT_mass', 'f', ('TSyearnr', ))
    ncAT_mass.long_name  = 'nitrate mass in atmospheric box'
    ncAT_mass.units      = 'kgN/m2'
    ncAT_d15N  = grpout.createVariable('AT_d15N', 'f', ('TSyearnr', ))
    ncAT_d15N.long_name  = 'd15N in nitrate'
    ncAT_d15N.units      = 'permil'
    ncAT_D17O  = grpout.createVariable('AT_D17O', 'f', ('TSyearnr', ))
    ncAT_D17O.long_name  = 'Atmosphere : D17O in nitrate'
    ncAT_D17O.units      = 'permil'
    ncAT_D17O_addO  = grpout.createVariable('AT_D17O_addO', 'f', ('TSyearnr', ))
    ncAT_D17O_addO.long_name  = 'Atmosphere : D17O in the additional oxygen atom transferred to NO2'
    ncAT_D17O_addO.units      = 'permil'
    ncAT_D17O_NO2_PSS  = grpout.createVariable('AT_D17O_NO2_PSS', 'f', ('TSyearnr', ))
    ncAT_D17O_NO2_PSS.long_name  = 'Atmosphere : D17O of NO2 after having applied the PSS'
    ncAT_D17O_NO2_PSS.units      = 'permil'
    ncFA_mass  = grpout.createVariable('FA_mass', 'f', ('TStotal_number', ))
    ncFA_mass.long_name  = 'Archived flux : nitrate mass in the archived flux'
    ncFA_mass.units      = 'kgN/m2/s'
    ncFA_d15N  = grpout.createVariable('FA_d15N', 'f', ('TStotal_number', ))
    ncFA_d15N.long_name  = 'Archived flux : d15N in nitrate'
    ncFA_d15N.units      = 'permil'
    ncFA_D17O  = grpout.createVariable('FA_D17O', 'f', ('TStotal_number', ))
    ncFA_D17O.long_name  = 'Archived flux : D17O in nitrate'
    ncFA_D17O.units      = 'permil'
    ncFA_thick  = grpout.createVariable('FA_thick', 'f', ('TStotal_number', ))
    ncFA_thick.long_name  = 'Archived flux : thickness of each layer'
    ncFA_thick.units      = 'meters'
    ncFA_date  = grpout.createVariable('FA_date', 'f', ('TStotal_number', ))
    ncFA_date.long_name  = 'Archived flux : date of the snow layer'
    ncFA_date.units      = 'date in years'
    ncFA_JO3  = grpout.createVariable('FA_JO3', 'f', ('TStotal_number', ))
    ncFA_JO3.long_name  = 'Archived flux : JO3 photolytic rate constant at time of the deposition of the layer'
    ncFA_JO3.units      = '/s'
    ncFA_conc  = grpout.createVariable('FA_conc', 'f', ('TStotal_number', ))
    ncFA_conc.long_name  = 'Archived flux : nitrate conc in the snow layer'
    ncFA_conc.units      = 'ngNO3-/g'
    ncFDall_mass  = grpout.createVariable('FDall_mass', 'f', ('TStotal_number', ))
    ncFDall_mass.long_name  = 'Deposited flux : nitrate mass deposited to the snow for ALL computation (not only the last year)'
    ncFDall_mass.units      = 'kgN/m2/s'
    ncFDall_d15N  = grpout.createVariable('FDall_d15N', 'f', ('TStotal_number', ))
    ncFDall_d15N.long_name  = 'Deposited flux : nitrate d15N deposited to the snow for ALL computation (not only the last year)'
    ncFDall_d15N.units      = 'permil'
    ncFDall_D17O  = grpout.createVariable('FDall_D17O', 'f', ('TStotal_number', ))
    ncFDall_D17O.long_name  = 'Deposited flux : nitrate D17O deposited to the snow for ALL computation (not only the last year)'
    ncFDall_D17O.units      = 'permil'
    ncFPall_mass  = grpout.createVariable('FPall_mass', 'f', ('TStotal_number', ))
    ncFPall_mass.long_name  = 'Photolyzed flux : nitrate mass emitted from snowpack for ALL computation (not only the last year)'
    ncFPall_mass.units      = 'kgN/m2/s'
    ncFPall_d15N  = grpout.createVariable('FPall_d15N', 'f', ('TStotal_number', ))
    ncFPall_d15N.long_name  = 'Photolyzed flux : nitrate d15N emitted from snowpack for ALL computation (not only the last year)'
    ncFPall_d15N.units      = 'permil'
    ncFPall_D17O  = grpout.createVariable('FPall_D17O', 'f', ('TStotal_number', ))
    ncFPall_D17O.long_name  = 'Photolyzed flux : nitrate D17O emitted from snowpack for ALL computation (not only the last year)'
    ncFPall_D17O.units      = 'permil'
    ncFEall_mass  = grpout.createVariable('FEall_mass', 'f', ('TStotal_number', ))
    ncFEall_mass.long_name  = 'Exported flux : nitrate mass transported off the snowpack for ALL computation (not only the last year)'
    ncFEall_mass.units      = 'kgN/m2/s'
    ncFEall_d15N  = grpout.createVariable('FEall_d15N', 'f', ('TStotal_number', ))
    ncFEall_d15N.long_name  = 'Exported flux : nitrate d15N transported off the snowpack for ALL computation (not only the last year)'
    ncFEall_d15N.units      = 'permil'
    ncFEall_D17O  = grpout.createVariable('FEall_D17O', 'f', ('TStotal_number', ))
    ncFEall_D17O.long_name  = 'Exported flux : nitrate D17O transported off the snowpack for ALL computation (not only the last year)'
    ncFEall_D17O.units      = 'permil'

    # Write values
    ncSL_depth[:]        = out.SL_depth
    ncSL_thick[:]        = out.SL_thick
    ncSL_mass[:]         = out.SL_mass
    ncSL_conc[:]         = out.SL_conc
    ncSL_d15N[:]         = out.SL_d15N
    ncSL_D17O[:]         = out.SL_D17O
    ncSN_depth[:]        = out.SN_depth
    ncSN_depth_mid[:]    = out.SN_depth_mid
    ncSN_thick[:]        = out.SN_thick
    ncSN_mass[:]         = out.SN_mass
    ncSN_conc[:]         = out.SN_conc
    ncSN_date[:]         = out.SN_date
    ncSN_d15N[:]         = out.SN_d15N
    ncSN_D17O[:]         = out.SN_D17O
    ncSN_eps15[:]        = out.SN_eps15
    ncSN_J14[:]          = out.SN_J14
    ncSN_JO3[:]          = out.SN_JO3
    ncSN_FP_mass_contri[:]      = out.SN_FP_mass_contri
    ncSN_FP_d15N_contri[:]      = out.SN_FP_d15N_contri
    ncSN_FP_D17O_contri[:]      = out.SN_FP_D17O_contri
    ncFD_mass[:]         = out.FD_mass
    ncFD_d15N[:]         = out.FD_d15N
    ncFD_D17O[:]         = out.FD_D17O
    ncFE_mass[:]         = out.FE_mass
    ncFE_d15N[:]         = out.FE_d15N
    ncFE_D17O[:]         = out.FE_D17O
    ncFP_mass[:]         = out.FP_mass
    ncFP_d15N[:]         = out.FP_d15N
    ncFP_D17O[:]         = out.FP_D17O
    ncFP_d15N_NOx[:]         = out.FP_d15N_NOx
    ncFP_D17O_NOx[:]         = out.FP_D17O_NOx
    ncAT_mass[:]         = out.AT_mass
    ncAT_d15N[:]         = out.AT_d15N
    ncAT_D17O[:]         = out.AT_D17O
    ncAT_D17O_addO[:]    = out.AT_D17O_addO
    ncAT_D17O_NO2_PSS[:] = out.AT_D17O_NO2_PSS
    ncFA_mass[:]         = out.FA_mass
    ncFA_d15N[:]         = out.FA_d15N
    ncFA_D17O[:]         = out.FA_D17O
    ncFA_thick[:]        = out.FA_thick
    ncFA_date[:]         = out.FA_date
    ncFA_JO3[:]          = out.FA_JO3
    ncFA_conc[:]         = out.FA_conc
    ncFDall_mass[:]      = out.FDall_mass
    ncFDall_d15N[:]      = out.FDall_d15N
    ncFDall_D17O[:]      = out.FDall_D17O
    ncFEall_mass[:]      = out.FEall_mass
    ncFEall_d15N[:]      = out.FEall_d15N
    ncFEall_D17O[:]      = out.FEall_D17O
    ncFPall_mass[:]      = out.FPall_mass
    ncFPall_d15N[:]      = out.FPall_d15N
    ncFPall_D17O[:]      = out.FPall_D17O

    ###################### Write inputs ####################################################
    # createVariable (variable's name, type, dimensions)
    # If 1 dimension : do not forget the ","
    ncaccu_mass  = grpinp.createVariable('accu_mass', 'f', ('TSyearnr', ))
    ncaccu_mass.long_name  = 'Input : mass of snow accumulated throughout the year (in each time step)'
    ncaccu_mass.units      = 'kg/m2/time_step'
    ncFS_mass  = grpinp.createVariable('FS_mass', 'f', ('TSyearnr', ))
    ncFS_mass.long_name  = 'Input : nitrate mass of stratospheric origin throughout the year (in each time step)'
    ncFS_mass.units      = 'kgN/m2/time_step'
    ncFS_d15N  = grpinp.createVariable('FS_d15N', 'f', ('TSyearnr', ))
    ncFS_d15N.long_name  = 'Input : d15N in nitrate of stratospheric origin throughout the year (in each time step)'
    ncFS_d15N.units      = 'permil'
    ncFS_D17O  = grpinp.createVariable('FS_D17O', 'f', ('TSyearnr', ))
    ncFS_D17O.long_name  = 'Input : D17O in nitrate of stratospheric origin throughout the year (in each time step)'
    ncFS_D17O.units      = 'permil'
    ncFT_mass  = grpinp.createVariable('FT_mass', 'f', ('TSyearnr', ))
    ncFT_mass.long_name  = 'Input : nitrate mass of long distance origin throughout the year (in each time step)'
    ncFT_mass.units      = 'kgN/m2/time_step'
    ncFT_d15N  = grpinp.createVariable('FT_d15N', 'f', ('TSyearnr', ))
    ncFT_d15N.long_name  = 'Input : d15N in nitrate of long distance origin throughout the year (in each time step)'
    ncFT_d15N.units      = 'permil'
    ncFT_D17O  = grpinp.createVariable('FT_D17O', 'f', ('TSyearnr', ))
    ncFT_D17O.long_name  = 'Input : D17O in nitrate of long distance origin throughout the year (in each time step)'
    ncFT_D17O.units      = 'permil'
    ncAT_height  = grpinp.createVariable('AT_height', 'f', ('TSyearnr', ))
    ncAT_height.long_name  = 'Input : height of the atmospheric boundary layer throughout the year (in each time step)'
    ncAT_height.units      = 'meters from ground'
    ncAT_mass  = grpinp.createVariable('AT_mass', 'f', ('TSyearnr', ))
    ncAT_mass.long_name  = 'Input : nitrate mass in atmosphere throughout the year (in each time step)'
    ncAT_mass.units      = 'kgN'
    ncAT_conc  = grpinp.createVariable('AT_conc', 'f', ('TSyearnr', ))
    ncAT_conc.long_name  = 'Input : nitrate concentration in the atmosphere'
    ncAT_conc.units      = 'ngNO3-.m-3'
    ncf_exp  = grpinp.createVariable('f_exp', 'f', ('TSyearnr', ))
    ncf_exp.long_name  = 'Input : horizontal exported fraction throughout the year (in each time step)'
    ncf_exp.units      = 'fraction of the positive fluxes to the atmosphere'
    ncO3colDU  = grpinp.createVariable('O3colDU', 'f', ('TSyearnr', ))
    ncO3colDU.long_name  = 'Input : stratospheric ozone column throughout the year (in each time step)'
    ncO3colDU.units      = 'Dobson Units (DU)'
    ncalpha  = grpinp.createVariable('alpha', 'f', ('TSyearnr', ))
    ncalpha.long_name  = 'Input : alpha coefficient in the PSS cycling of NO3 in the atmosphere throughout the year (in each time step)'
    ncalpha.units      = 'no unit'
    ncD17O_O3  = grpinp.createVariable('D17O_O3', 'f', ('TSyearnr', ))
    ncD17O_O3.long_name  = 'Input : D17O of ozone in atmosphere throughout the year (in each time step)'
    ncD17O_O3.units      = 'permil'
    ncx_fact  = grpinp.createVariable('x_fact', 'f', ('TSyearnr', ))
    ncx_fact.long_name  = 'Input : x factor used in the calculation of D17O in OH in atmosphere throughout the year (in each time step)'
    ncx_fact.units      = 'no unit'
    ncSZA_mean  = grpinp.createVariable('SZA_mean', 'f', ('TSyearnr', ))
    ncSZA_mean.long_name  = 'Input : mean Solar Zenith Angle throughout the year (in each time step)'
    ncSZA_mean.units      = 'Degrees. 0 is zenith. 90 is at horizon'
    ncair_temp  = grpinp.createVariable('air_temp', 'f', ('TSyearnr', ))
    ncair_temp.long_name  = 'Input : Atmospheric air temperature throughout the year (in each time step)'
    ncair_temp.units      = 'Kelvin'
    ncair_pres  = grpinp.createVariable('air_pres', 'f', ('TSyearnr', ))
    ncair_pres.long_name  = 'Input : Atmospheric air temperature throughout the year (in each time step)'
    ncair_pres.units      = 'Pressure'

    # Write values
    ncaccu_mass[:]   = inp.accu_mass
    ncFS_mass[:]     = inp.FS_mass
    ncFS_d15N[:]     = inp.FS_d15N
    ncFS_D17O[:]     = inp.FS_D17O
    ncFT_mass[:]     = inp.FT_mass
    ncFT_d15N[:]     = inp.FT_d15N
    ncFT_D17O[:]     = inp.FT_D17O
    ncAT_height[:]   = inp.AT_height
    ncAT_mass[:]     = inp.AT_mass
    ncAT_conc[:]     = inp.AT_conc
    ncf_exp[:]       = inp.f_exp
    ncO3colDU[:]     = inp.O3colDU
    ncalpha[:]       = inp.alpha
    ncD17O_O3[:]     = inp.D17O_O3
    ncx_fact[:]      = inp.x_fact
    ncSZA_mean[:]    = inp.SZA_mean
    ncair_temp[:]    = inp.air_temp
    ncair_pres[:]    = inp.air_pres

    ###################### Write options #################################################
    # createVariable (variable's name, type, dimensions)
    # If 1 dimension : do not forget the ","
    ncNO2_oxidation  = grpopt.createVariable('NO2_oxidation', 'i', ('DIMzero', ))
    ncNO2_oxidation.long_name   = 'Option : NO2 oxidation option : 1=local oxi, 0=no local oxi'
    ncNO2_oxidation.units       = 'no unit'
    ncdep_diff_scenario  = grpopt.createVariable('dep_diff_scenario', 'i', ('DIMzero', ))
    ncdep_diff_scenario.long_name   = 'Option : Deposition and diffusion option : 0=100% top layer no diffusion, 1=100% on top layer + diffusion at top only, 2=100% on top layer + diffusion in the whole snowpack'
    ncdep_diff_scenario.units       = 'no unit'
    ncDiff  = grpopt.createVariable('Diff', 'f', ('DIMzero', ))
    ncDiff.long_name   = 'Option : Diffusion coefficient to be used in case 2'
    ncDiff.units       = 'm2/s'
    ncO3col_variations  = grpopt.createVariable('O3col_variations', 'i', ('DIMzero', ))
    ncO3col_variations.long_name   = 'Option : variations in O3col inputs : 0=no variations from year to year, read the time series in input file, 1=O3col is varied from year to year and the values are looked for in a specific input'
    ncO3col_variations.units       = 'no unit'
    ncsastrugi_effect  = grpopt.createVariable('sastrugi_effect', 'i', ('DIMzero', ))
    ncsastrugi_effect.long_name   = 'Option : sastrugy effect : 0=no sastrugi deposition, 1=a sastrugi deposits at given time step and given thickness'
    ncsastrugi_effect.units       = 'no unit'
    ncsastrugi_thick  = grpopt.createVariable('sastrugi_thick', 'f', ('DIMzero', ))
    ncsastrugi_thick.long_name   = 'Option : sastrugy thickness'
    ncsastrugi_thick.units       = 'meter'
    ncsastrugi_TS  = grpopt.createVariable('sastrugi_TS', 'i', ('DIMzero', ))
    ncsastrugi_TS.long_name   = 'Option : time step of sastrugi deposition'
    ncsastrugi_TS.units       = 'time step'

    # Write values
    ncNO2_oxidation[:]      = option.NO2_oxidation
    ncdep_diff_scenario[:]  = option.dep_diff_scenario
    ncDiff[:]               = option.Diff
    ncO3col_variations[:]   = option.O3col_variations
    ncsastrugi_effect[:]    = option.sastrugi_effect
    ncsastrugi_thick[:]     = option.sastrugi_thick
    ncsastrugi_TS[:]        = option.sastrugi_TS

    ###################### Write parameters #################################################
    # createVariable (variable's name, type, dimensions)
    # If 1 dimension : do not forget the ","
    ncSN_set_thick  = grppar.createVariable('SN_set_thick', 'f', ('DIMzero', ))
    ncSN_set_thick.long_name   = 'Param : set thickness of snowpack'
    ncSN_set_thick.units       = 'meters'
    ncSN_set_height  = grppar.createVariable('SN_set_height', 'f', ('DIMzero', ))
    ncSN_set_height.long_name   = 'Param : set height of snowpack'
    ncSN_set_height.units       = 'meters'
    ncnb_yr  = grppar.createVariable('nb_yr', 'i', ('DIMzero', ))
    ncnb_yr.long_name   = 'Param : number of computed years in the simulation'
    ncnb_yr.units       = 'years'
    ncnb_ts  = grppar.createVariable('nb_ts', 'i', ('DIMzero', ))
    ncnb_ts.long_name   = 'Param : number of time steps in one year'
    ncnb_ts.units       = 'no unit'
    nccp_dur  = grppar.createVariable('cp_dur', 'i', ('DIMzero', ))
    nccp_dur.long_name   = 'Param : number of time steps in simulation'
    nccp_dur.units       = 'no unit'
    ncDt  = grppar.createVariable('Dt', 'i', ('DIMzero', ))
    ncDt.long_name   = 'Param : time step duration'
    ncDt.units       = 'seconds'
    nccp_S  = grppar.createVariable('cp_S', 'f', ('DIMzero', ))
    nccp_S.long_name   = 'Param : computed column surface'
    nccp_S.units       = 'm2'
    ncSN_d  = grppar.createVariable('SN_d', 'f', ('DIMzero', ))
    ncSN_d.long_name   = 'Param : snow density'
    ncSN_d.units       = 'no unit'
    ncSN_rho  = grppar.createVariable('SN_rho', 'f', ('DIMzero', ))
    ncSN_rho.long_name   = 'Param : snow volumic mass'
    ncSN_rho.units       = 'kg/m3'
    ncSN_SSA  = grppar.createVariable('SN_SSA', 'f', ('DIMzero', ))
    ncSN_SSA.long_name   = 'Param : Surface Specific Area of snow'
    ncSN_SSA.units       = 'm2/kg'
    ncSN_tau  = grppar.createVariable('SN_tau', 'f', ('DIMzero', ))
    ncSN_tau.long_name   = 'Param : tortuosity of snow'
    ncSN_tau.units       = 'no unit'
    ncPhi  = grppar.createVariable('Phi', 'f', ('DIMzero', ))
    ncPhi.long_name   = 'Param : photolytic quantum yield'
    ncPhi.units       = 'molecules/photon'
    ncf_cage  = grppar.createVariable('f_cage', 'f', ('DIMzero', ))
    ncf_cage.long_name   = 'Param : cage effect parameter'
    ncf_cage.units       = 'no unit'
    ncD17O_water  = grppar.createVariable('D17O_water', 'f', ('DIMzero', ))
    ncD17O_water.long_name   = 'Param : D17O of water in snow'
    ncD17O_water.units       = 'permil'
    ncJO3_threshold  = grppar.createVariable('JO3_threshold', 'f', ('DIMzero', ))
    ncJO3_threshold.long_name   = 'Param : JO3 threshold value. For this value : 100% oxidation by OH'
    ncJO3_threshold.units       = '/s'
    ncSL_thick  = grppar.createVariable('SL_thick', 'f', ('DIMzero', ))
    ncSL_thick.long_name   = 'Param : skinlayer thickness'
    ncSL_thick.units       = 'meters'
    ncaccu_mass_yr  = grppar.createVariable('accu_mass_yr', 'f', ('DIMzero', ))
    ncaccu_mass_yr.long_name   = 'Param : snow annual accumulated mass'
    ncaccu_mass_yr.units       = 'kg/m2/yr'
    ncFPI_mass_yr  = grppar.createVariable('FPI_mass_yr', 'f', ('DIMzero', ))
    ncFPI_mass_yr.long_name   = 'Param : nitrate mass in the primary inputs (strato + long distance transport)'
    ncFPI_mass_yr.units       = 'kgN/m2/yr'
    ncFS_FPI_frac  = grppar.createVariable('FS_FPI_frac', 'f', ('DIMzero', ))
    ncFS_FPI_frac.long_name   = 'Param : fraction of the primary inputs that occures as stratospheric denitrification'
    ncFS_FPI_frac.units       = 'no unit'
    ncFT_FPI_frac  = grppar.createVariable('FT_FPI_frac', 'f', ('DIMzero', ))
    ncFT_FPI_frac.long_name   = 'Param : fraction of the primary inputs that occures as long distance transport'
    ncFT_FPI_frac.units       = 'no unit'
    ncphotic_zone  = grppar.createVariable('photic_zone', 'f', ('DIMzero', ))
    ncphotic_zone.long_name   = 'Param : photic zone compression factor. E-folding becomes efold * photic_zone'
    ncphotic_zone.units       = 'no unit'
    ncactinic_flux  = grppar.createVariable('actinic_flux', 'f', ('DIMzero', ))
    ncactinic_flux.long_name   = 'Param : actinic flux change factor. 1 is normal. <1 is less incoming actinic flux'
    ncactinic_flux.units       = 'no unit'
    ncAT_eps15_dep  = grppar.createVariable('AT_eps15_dep', 'f', ('DIMzero', ))
    ncAT_eps15_dep.long_name   = 'Param : 15N/14N isotope fractionation at deposition'
    ncAT_eps15_dep.units       = 'permil'

    # Write values
    ncnb_yr[:]           = param.nb_yr
    ncnb_ts[:]           = param.nb_ts
    nccp_dur[:]          = param.cp_dur
    nccp_S[:]            = param.cp_S
    ncDt[:]              = param.Dt
    ncSN_d[:]            = param.SN_d
    ncSN_rho[:]          = param.SN_rho
    ncSN_SSA[:]          = param.SN_SSA
    ncSN_tau[:]          = param.SN_tau
    ncPhi[:]             = param.Phi
    ncf_cage[:]          = param.f_cage
    ncD17O_water[:]      = param.D17O_water
    ncJO3_threshold[:]   = param.JO3_threshold
    ncSL_thick[:]        = param.SL_thick
    ncphotic_zone[:]     = param.photic_zone
    ncactinic_flux[:]    = param.actinic_flux
    ncaccu_mass_yr[:]    = param.accu_mass_yr
    ncFPI_mass_yr[:]     = param.FPI_mass_yr
    ncFS_FPI_frac[:]     = param.FS_FPI_frac
    ncFT_FPI_frac[:]     = param.FT_FPI_frac
    ncAT_eps15_dep[:]    = param.AT_eps15_dep
    ncSN_set_thick[:]    = param.SN_set_thick
    ncSN_set_height[:]   = param.SN_set_height

    ### Close file
    nc.close()
    
#################################################################################################
FS_mass_mean = param.FS_FPI_frac * param.FPI_mass_yr
FS_d15N_mean = (inp.FS_mass_rep*inp.FS_d15N).sum()
FS_D17O_mean = (inp.FS_mass_rep*inp.FS_D17O).sum()
FT_mass_mean = param.FT_FPI_frac * param.FPI_mass_yr
FT_d15N_mean = (inp.FT_mass_rep*inp.FT_d15N).sum()
FT_D17O_mean = (inp.FT_mass_rep*inp.FT_D17O).sum()
FE_mass_mean = out.FE_mass.sum()*param.Dt
FE_d15N_mean = ((out.FE_d15N[0, :]*out.FE_mass[0, :]).sum()*param.Dt)/FE_mass_mean
FE_D17O_mean = ((out.FE_D17O[0, :]*out.FE_mass[0, :]).sum()*param.Dt)/FE_mass_mean
FA_mass_mean = out.FA_mass[0, param.cp_dur-param.nb_ts:param.cp_dur].sum()
FA_d15N_mean = (out.FA_mass[0, param.cp_dur-param.nb_ts:param.cp_dur]*out.FA_d15N[0, param.cp_dur-param.nb_ts:param.cp_dur]).sum()/FA_mass_mean
FA_D17O_mean = (out.FA_mass[0, param.cp_dur-param.nb_ts:param.cp_dur]*out.FA_D17O[0, param.cp_dur-param.nb_ts:param.cp_dur]).sum()/FA_mass_mean
FA_conc_mean = FA_mass_mean*(param.M_NO3/param.M_N)*param.cp_S*1e9/(out.FA_thick[0, param.cp_dur-param.nb_ts:param.cp_dur].sum()*param.cp_S*param.SN_rho)
FP_mass_mean = out.FP_mass.sum()*param.Dt
if FP_mass_mean == 0:
    FP_d15N_mean = FP_D17O_mean = 0
else:
    FP_d15N_mean = ((out.FP_d15N[0, :]*out.FP_mass[0, :]).sum()*param.Dt)/FP_mass_mean
    FP_D17O_mean = ((out.FP_D17O[0, :]*out.FP_mass[0, :]).sum()*param.Dt)/FP_mass_mean
FD_mass_mean = out.FD_mass.sum()*param.Dt
FD_d15N_mean = ((out.FD_d15N[0, :]*out.FD_mass[0, :]).sum()*param.Dt)/FD_mass_mean
FD_D17O_mean = ((out.FD_D17O[0, :]*out.FD_mass[0, :]).sum()*param.Dt)/FD_mass_mean

FP_amount_max = out.FP_mass_amount[0, :].max()
FP_mol_max    = out.FP_mass_mol[0, :].max()

print( '----------------------------------------')
print( 'Mean annual fluxes :')
print( 'FS mass (10-6 kgN.m-2.a-1), d15N (permil), D17O(permil) : ' + str(round(FS_mass_mean*1e6, 3)) + ', ' + str(FS_d15N_mean) + ', ' + str(FS_D17O_mean))
print( 'FT mass (10-6 kgN.m-2.a-1), d15N (permil), D17O(permil) : ' + str(round(FT_mass_mean*1e6, 3)) + ', ' + str(FT_d15N_mean) + ', ' + str(FT_D17O_mean))
print( 'FE mass (10-6 kgN.m-2.a-1), d15N (permil), D17O(permil) : ' + str(round(FE_mass_mean*1e6, 3)) + ', ' + str(FE_d15N_mean) + ', ' + str(FE_D17O_mean))
print( 'FP mass (10-6 kgN.m-2.a-1), d15N (permil), D17O(permil) : ' + str(round(FP_mass_mean*1e6, 3)) + ', ' + str(FP_d15N_mean) + ', ' + str(FP_D17O_mean))
print( 'FD mass (10-6 kgN.m-2.a-1), d15N (permil), D17O(permil) : ' + str(round(FD_mass_mean*1e6, 3)) + ', ' + str(FD_d15N_mean) + ', ' + str(FD_D17O_mean))
print( 'FA mass (10-6 kgN.m-2.a-1), d15N (permil), D17O(permil) : ' + str(round(FA_mass_mean*1e6, 3)) + ', ' + str(FA_d15N_mean) + ', ' + str(FA_D17O_mean))

print( '----------------------------------------')
print( 'Annual mass balance in atmospheric box :')
print( 'Mass IN  : ' + str(round((FS_mass_mean+FT_mass_mean+FP_mass_mean)*10**6, 3)) + ' 10-6 kgN.m-2.a-1')
print( 'Mass OUT : ' + str(round((FE_mass_mean+FD_mass_mean)*10**6, 3)) + ' 10-6 kgN.m-2.a-1')
print( 'Mass variation between two last years : ' + str(round((out.AT_mass_int_YR[0, param.nb_yr-1]-out.AT_mass_int_YR[0, param.nb_yr-2])/out.AT_mass_int_YR[0, param.nb_yr-1]*100, 3)) + ' %')
print( 'Mass * d15N IN  : ' + str(round((FS_mass_mean*FS_d15N_mean+FT_mass_mean*FT_d15N_mean+FP_mass_mean*FP_d15N_mean)*10**6, 3)) + ' 10-6 kgN.permil.m-2.a-1')
print( 'Mass * d15N OUT : ' + str(round((FE_mass_mean*FE_d15N_mean+FD_mass_mean*FD_d15N_mean)*10**6, 3)) + ' 10-6 kgN.permil.m-2.a-1')
print( 'Mass * D17O IN  : ' + str(round((FS_mass_mean*FS_D17O_mean+FT_mass_mean*FT_D17O_mean+FP_mass_mean*FP_D17O_mean)*10**6, 3)) + ' 10-6 kgN.permil.m-2.a-1')
print( 'Mass * D17O OUT : ' + str(round((FE_mass_mean*FE_D17O_mean+FD_mass_mean*FD_D17O_mean)*10**6, 3)) + ' 10-6 kgN.permil.m-2.a-1')

print( '---------------------------------')
print( 'Annual mass balance in snowpack :')
print( 'Mass IN  : ' + str(round((FD_mass_mean)*10**6, 3)) + ' 10-6 kgN.m-2.a-1')
print( 'Mass OUT : ' + str(round((FP_mass_mean+FA_mass_mean)*10**6, 3)) + ' 10-6 kgN.m-2.a-1')
print( 'Mass variation between two last years : ' + str(round((out.SN_mass_int_YR[0, param.nb_yr-1]-out.SN_mass_int_YR[0, param.nb_yr-2])/out.SN_mass_int_YR[0, param.nb_yr-1]*100, 3)) + ' %')
print( 'Mass * d15N IN  : ' + str(round((FD_mass_mean*FD_d15N_mean)*10**6, 3)) + ' 10-6 kgN.permil.m-2.a-1')
print( 'Mass * d15N OUT : ' + str(round((FP_mass_mean*FP_d15N_mean+FA_mass_mean*FA_d15N_mean)*10**6, 3)) + ' 10-6 kgN.permil.m-2.a-1')
print( 'Mass * D17O IN  : ' + str(round((FD_mass_mean*FD_D17O_mean)*10**6, 3)) + ' 10-6 kgN.permil.m-2.a-1')
print( 'Mass * D17O OUT : ' + str(round((FP_mass_mean*FP_D17O_mean+FA_mass_mean*FA_D17O_mean)*10**6, 3)) + ' 10-6 kgN.permil.m-2.a-1')
print( '---------------------------------')
print( "Mean archived concentration in snow : " + str(round(FA_conc_mean, 2)) + " ngNO3/g")
print( "Median NO2 flux for the mean SZA    : " + str(round(FP_amount_max, 2)) + " nmol.m-2.h^-1")
print( "Median NO2 flux for the mean SZA    : " + str(round(FP_mol_max*10**-12, 2)) + " 10^12 molecules.m-2.s-1")
 

