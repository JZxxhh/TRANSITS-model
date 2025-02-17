# -*- coding: cp1252 -*-
#################################################################################################
#################################################################################################
#################################################################################################
### This is the TRANSITS (TRansfert of Atmospheric Nitrate Stable Isotopes To the Snow) model ###
#################################################################################################
#################################################################################################
#################################################################################################

#################################################################################################
### IMPORT MODULES ##############################################################################
# We here import commonly used modules
import numpy as np                  # array module                        
# We here import home-made modules and all the functions contained in them
import functions as f
import model as m
import time
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")   # Ignore warnings like division by 0

#################################################################################################
# Model's version
version        = 'WD 1.0'
print('### This is TRANSITS model ' + version)
print('-------------------------------')
# Input time file to read
input_file_to_read = 'input_time_REALISTIC_simulation.xlsx'
print('The file \"' + input_file_to_read + '\" is used to generate the inputs')

###################################################################################################
##### PARAMETERS ##################################################################################
# PARAMETERS object : parameters are constant values and are not scenarios as are the inputs
### Creation of the object
param = f.Params()
### Fill the object with the minimum parameters
### The other parameters are given BELOW (around line 420)
param.nb_ts =                   52             # time steps                number of time steps in one year, Spin-up of 27 days, run1=28*24=672, run2=30*24=720, run3=33*24=792, run4=35*24=840, run5=37*24=888, run6=39*24=936, run7=42*24=1008, run8=43*24=1032, run9=44*24=1056, run10=48*24=1152, run11=50*24=1200, run12=55*24=1320   
param.cp_S =                    1.              # m^2                       computed column surface
# Below are constants
param.M_N =                     14.             # g.mol^(-1)                nitrogen molar mass
param.M_O =                     16.             # g.mol^(-1)                oxygen molar mass
param.calc_M_NO3()           #CALCULATED        g.mol^(-1)                nitrate molar mass
param.N_A =                     6.02214179E+23  # molecules.mol^(-1)        Avogadro's number
param.R   =                     8.3144621       # J·K-1·mol-1               Gas constant
param.Pi =                      3.1415926535898
param.out_fric =                1 # export fraction

#################################################################################################
### OPTIONS    ##################################################################################
# OPTIONS object : options are integers which set the model to specific states
### Creation of the object
option = f.Options()

#################################################################################################
### INPUTS CREATION #############################################################################
### Inputs vary according to the scenarios choosen for the simulation
# Creation of the object
inp = f.Inputs()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Read and import SZA_rep file
load = open("inputs_SZA_rep/inputs_SZA_rep.txt", "r")
inp.SZA_rep = np.genfromtxt(load, skip_header=0, delimiter="")
inp.SZA_rep = inp.SZA_rep * 60  # converted from minutes to seconds    

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Read and import J_O3 file
# not a very clever method here as I have to change the TCO object everytime for a new sites
# READ for every ozone column, FLIP (light to right), TRANSPOSE and TRANSFER to arrays
#file = open("inputs_TUV\inputs_JO3.txt","r")
file = open("inputs_TUV/inputs_JO3.txt","r")
importJO3 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
# Create object in the input object 41*1
inp.JO3 = f.Ozone_column()
inp.JO3.DU150 = importJO3[1, :]
inp.JO3.DU175 = importJO3[2, :]
inp.JO3.DU200 = importJO3[3, :]
inp.JO3.DU225 = importJO3[4, :]
inp.JO3.DU250 = importJO3[5, :]
inp.JO3.DU275 = importJO3[6, :]
inp.JO3.DU300 = importJO3[7, :]
inp.JO3.DU325 = importJO3[8, :]
inp.JO3.DU350 = importJO3[9, :]
inp.JO3.DU375 = importJO3[10, :]
inp.JO3.DU400 = importJO3[11, :]
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Read and import J_NO2 file
# READ for every ozone column, FLIP (left to right), TRANSPOSE and TRANSFER to arrays
#file = open("inputs_TUV\inputs_JNO2.txt","r")
file = open("inputs_TUV/inputs_JNO2.txt","r") 
importJNO2 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
# Create object in the input object
inp.JNO2 = f.Ozone_column()
inp.JNO2.DU150 = importJNO2[1, :]
inp.JNO2.DU175 = importJNO2[2, :]
inp.JNO2.DU200 = importJNO2[3, :]
inp.JNO2.DU225 = importJNO2[4, :]
inp.JNO2.DU250 = importJNO2[5, :]
inp.JNO2.DU275 = importJNO2[6, :]
inp.JNO2.DU300 = importJNO2[7, :]
inp.JNO2.DU325 = importJNO2[8, :]
inp.JNO2.DU350 = importJNO2[9, :]
inp.JNO2.DU375 = importJNO2[10, :]
inp.JNO2.DU400 = importJNO2[11, :]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Read and import J's files
inp.J14NO3 = f.Ozone_column()
inp.J15NO3 = f.Ozone_column()
# Read header (depths) only
load = open("inputs_TUV/inputs_J14NO3_200DU.txt", "r")
read_for_header = np.genfromtxt(load, skip_header=0, delimiter="")
inp.J14NO3.header = read_for_header[0, :]
# READ for every ozone column, FLIP (light to right), TRANSPOSE and TRANSFER to arrays .120*41,depth*SZA
file = open("inputs_TUV/inputs_J14NO3_150DU.txt","r")
inp.J14NO3.DU150 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J15NO3_150DU.txt","r")
inp.J15NO3.DU150 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J14NO3_175DU.txt","r")
inp.J14NO3.DU175 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J15NO3_175DU.txt","r")
inp.J15NO3.DU175 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J14NO3_200DU.txt","r")
inp.J14NO3.DU200 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J15NO3_200DU.txt","r")
inp.J15NO3.DU200 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J14NO3_225DU.txt","r")
inp.J14NO3.DU225 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J15NO3_225DU.txt","r")
inp.J15NO3.DU225 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J14NO3_250DU.txt","r")
inp.J14NO3.DU250 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J15NO3_250DU.txt","r")
inp.J15NO3.DU250 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J14NO3_275DU.txt","r")
inp.J14NO3.DU275 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J15NO3_275DU.txt","r")
inp.J15NO3.DU275 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J14NO3_300DU.txt","r")
inp.J14NO3.DU300 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J15NO3_300DU.txt","r")
inp.J15NO3.DU300 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J14NO3_325DU.txt","r")
inp.J14NO3.DU325 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J15NO3_325DU.txt","r")
inp.J15NO3.DU325 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J14NO3_350DU.txt","r")
inp.J14NO3.DU350 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J15NO3_350DU.txt","r")
inp.J15NO3.DU350 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J14NO3_375DU.txt","r")
inp.J14NO3.DU375 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J15NO3_375DU.txt","r")
inp.J15NO3.DU375 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J14NO3_400DU.txt","r")
inp.J14NO3.DU400 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))
file = open("inputs_TUV/inputs_J15NO3_400DU.txt","r")
inp.J15NO3.DU400 = np.fliplr(np.transpose(np.genfromtxt(file, skip_header=1, delimiter="")))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Temporal series of inputs are found in the input_file_to_read
f.Read_inputs(inp, param, input_file_to_read)

snow_accu_frac = 0.019*np.ones(52)   #definded snow accumulation fraction 
#snow_accu_frac = np.zeros((52))+1/52
#################################################################################################

print('Single run started on ' + time.ctime(time.time()))


# Read inputs again
f.Read_inputs(inp, param, input_file_to_read)

### Set values to each parameter in the general case :      
param.nb_yr =                   15          # years                     number of computed years in the simulation (25 orig)
param.nb_end =                  26          # week                       number of the last computed week
param.calc_cp_dur()            #CALCULATED        no unit                    number of time steps in simulation
param.Dt =      (7*52*24*60*60)/param.nb_ts      # seconds                 time step duration (default = 1 year / 52) 10107 is the number of minutes in such a time step
#param.Dt =                      3600             # seconds                   time step duration (default = 1 year / 52) 10107 is the number of minutes in such a time step, run1=28*24=672, run2=30*24=720, run3=33*24=792, run4=35*24=840, run5=37*24=888, run6=39*24=936, run7=42*24=1008, run8=43*24=1032, run9=44*24=1056, run10=48, run11=50*24=1200, run12=55*24=1320
param.cp_S =                    1              # m^2                      computed column surface
param.SN_set_height =           1              # meters                height of the total active snowpack (pit1=18 cm, pit2=9cm, pit3=14 cm, pit4=18cm, pit5=11 cm, pit6=18 cm, pit7=20 cm, pit8=15 cm, pit9=14 cm (try 9cm), pit10=24 cm (try 19 cm), pit11=21 cm (try 16 cm), pit12=20 cm (try 15 cm))
param.SN_set_thick =            0.001            # meters                     thickness of each layer
param.SN_d =                    0.35            # adimensional             snow density
param.calc_rho()               #CALCULATED        kg.m^(-3)                  snow volumic mass
param.cal_SN_SSA()             # CALCULATED                m2/kg           Surface Specific Area of snow
param.SN_tau =                  1.08            # no unit                   Tortuosity (here taken from Pinzer 2010, otherwise fresh = 1.08, old = 1.15) 
param.Phi =                     0.0054           #0.002           # molecules.photon^(-1)     (0.0028) photolytic quantum yield
param.f_cage =                  0.15              # fraction                  cage effect parameter
param.D17O_water =              0.              # permil                    D17O of water in snow
param.SL_thick =                0.001           # meters                    skinlayer thickness
param.photic_zone =             1.              # adimensional              photic zone compression factor (1 is 10 cm e-folding as usual)
param.actinic_flux =            1.              # adimensional              change in the incoming actinic flux (because of solar variability or orbital parameters. Default is 1)
param.accu_mass_yr =            112             # kg.m^(-2).yr^(-1)         accumulated snow mass per year (change to 18 when snow, poss. 15),  0.00000001 otherwise
param.FPI_mass_yr =             2.2e-6  #8.06*1E-6          # kgN.m^(-2).yr^(-1)        TOTAL primary input mass deposited per year (pit1=4.6E-6, pit2=4.9E-6, pit3=5.4E-6, pit4=5.7E-6, pit5=6.0E-6 , pit6=6.4E-6, pit7=6.9E-6, pit8=7.0E-6, pit9=7.2E-6, pit10=7.8E-6, pit11=8.2E-6, pit12=9.5E-6) 
param.FS_FPI_frac =             0.0             # adimensional              fraction of the primary input that deposits as stratospheric nitrate per year
param.FT_FPI_frac =             1.0             # adimensional              fraction of the primary input that comes as long distance nitrate per year
param.AT_eps15_dep =            10.              # permil                    15N/14N isotope fractionation at deposition of nitrate from AT to snow
param.AT_BrO_conc =             2E-12         # mol/mol                   BrO mixing ratio in atmosphere
param.AT_CO_conc =              40E-9           # mol/mol                   [Not used] CO mixing ratio in atmosphere
param.AT_CH4_conc =             1.8E-6          # mol/mol                   [Not used] CH4 mixing ratio in atmosphere
param.AT_RH =                   70              # %                           [Not used] relative humidity in atmosphere
param.AT_HO2_conc_scale   =     1               # adimensional              Scaling factor for HO2 concentrations (modified in sensitivity tests)
param.AT_CH3O2_conc_scale =     1               # adimensional                Scaling factor for HO2 concentrations (modified in sensitivity tests)

option.local_recycling       =           1               # adimensional              Nitrate recycling option :
                                                                                     # 1    = recycling is 100 % ON. Chemistry is much faster than transport
                                                                                     #        Therefore, HNO3 is reformed locally
                                                                                     # 0    = no recycling. Chemistry is much slower than transport
                                                                                     #        and NO2 is completely lost from the atmospheric box
option.NO2_cycling           =           0               # adimensional              NO2 cycling (Leighton) option :
                                                                                     # 0    = cycling is OFF, NO2 keeps the snowpack D17O signature
                                                                                     # 1    = cycling is ON, NO2 inherits D17O signature from local cycling (assuming PSS)
                                                                                     # 2    = cycling is ON, set alpha value to option.NO2_alpha (constant value throughout the year)
                                                                                     # 3    = cycling is ON, set the minimum alpha value to option.NO2_alpha
                                                                                     # 9999 = cycling is ON but D17O(NO2) is set to 0 (for sensitivity tests)
                                                                                     # x<0  = cycling is ON but D17O(NO2) is set to abs(x) (for sensitivity tests)
option.NO2_set_alpha         =         0.7               # adimensional                Value of constant alpha (when option.NO2_cycling = 2)
                                                                                     # Value of minimum alpha (when option.NO2_cycling = 3)    
option.NO2_oxidation         =           1               # adimensional              NO2 local oxidation option :
                                                                                     # 0    = oxidation is OFF, HNO3 keeps the NO2 signature
                                                                                     # 1    = oxidation is ON, HNO3 gets 2/3 of the D17O of NO2 and 1/3 of D17O of the additionnal O atom
                                                                                     #        and add. O can come from OH or O3 depending on reaction rates
                                                                                     # 2    = force oxidation by OH only. D17O(add. O) = D17O(OH)
                                                                                     # 3    = force oxidation by O3 only. D17O(add. O) = D17O(O3)
                                                                                     # 9999 = oxidation is ON but D17O(add. O) is set to 0 (for sensitivity tests)
                                                                                     # x<0  = oxidation is ON but D17O(add. O) is set to abs(x) (for sensitivity tests)
option.dep_diff_scenario     =           0              # Diffusion option : 0 = 100% deposition on top layer + no diffusion at all,
                                                         #             1 = 100% deposition on top layer + diffusion from this deposited nitrate only (HNO3 diffusion)
                                                         #             2 = 100% deposition on top layer + diffusion in 1-m snowpack at all time step with fixed Diff coefficient (NO3- diffusion)
                                                         #             3 = 100% deposition on top layer + diffusion in 1-m snowpack at all time step with Temp-dependent Diff coefficient  (NO3- diffusion)
                                                         #             4 = 100% deposition on top layer + diffusion in 1-m snowpack at all time step with T/P/rho/SSA/tau-dependent Diff coefficient (NO3- diffusion)
option.Diff                  =   1e-11              # Diffusion coefficient to be used in case 2 (m2 s-1)
option.sastrugi_TS           =         1032              # time step number (0 is 21 June)  
option.SN_accu               =          0                 # option = 0: use snow accumulation rates in input file; else used defined input
option.phi                   =         1                 # option = 1: use quantum yield value as given value; else use phi(T) relationship
                                                       

if option.SN_accu == 1:
    inp.accu_mass_rep = snow_accu_frac 

### EXTRA PARAMETERS CALCULATION
# We here calculate the number of snow layers in the active snowpack zone
param.calc_layers_nb()

### INPUTS PREPARATION
# Some input time series are further calculated
inp.calc_accu_mass(param)   # snow accumulated mass in each time step              (kg.m^(-2))
inp.calc_FS_mass(param)     # stratospheric denitrification flux at each time step (kgN.m^(-2).s^(-1))
inp.calc_FT_mass(param)     # long distance transportflux at each time step        (kgN.m^(-2).s^(-1))

### RESET STATE VARIABLES + INITIALIZE
# Build a snowpack with a set depth and thickness for each layer
  # Creation of the object
state = f.State(param,inp)

### RESET OUTPUTS VARIABLES + INITIALIZE
  # Creation of the object
out   = f.Outputs(param)                                                                                                                                                                                                                                                    
#########################################################################################
### MODEL'S RUN
# Below, we run the model
print('ready for simulation!')
m.model(param, option, state, inp, out)
#########################################################################################

### ADDITIONAL OUTPUTS CALCULATIONS
# Calculates the mid depth of each layer
out.mid_depth = out.SN_depth-0.5*out.SN_thick
# Calculates nitrate concentrations in snow
out.SN_conc = (param.M_NO3 / param.M_N) * out.SN_mass / (out.SN_thick * param.cp_S * param.SN_rho) * pow(10, 9)
# Converts the photolytically produced NO2 flux (originally in kgN.m^(-2).s^(-1)) in molecules.m^(-2).s^(-1)
out.calc_NO2_flux_molecule(param)
# Converts the photolytically produced NO2 flux (originally in kgN.m^(-2).s^(-1)) in nmol.m^(-2).h^(-1)
out.calc_NO2_flux_nmol(param)
    #f.plot_weekly(out.SN_conc,out.SN_thick,out.SN_d15N,out.SN_D17O,out.SN_date,inp.accu_mass_rep)
#################################################################################################

# write output
new_ma = np.zeros((len(out.mid_depth),6))
new_ma[:,0] = out.mid_depth
new_ma[:,1] = out.SN_conc
new_ma[:,2] = out.SN_d15N
new_ma[:,3] = out.SN_D17O
new_ma[:,4] = out.SN_date
new_ma[:,5] = out.SN_thick
d15N_total = (out.SN_thick*out.SN_conc*out.SN_d15N).sum()/(out.SN_thick*out.SN_conc).sum()

##################################################################################################
# calculate the remaining fraction from FA and FD
nb_tt = out.FD_mass.shape[1]
FA_mass = np.zeros((nb_tt))
FA_d15N = np.zeros((nb_tt))
FA_D17O = np.zeros((nb_tt))

for i in range(nb_tt):
    idx = np.where(out.SN_date==i) # find the index
    FA = out.SN_thick*param.SN_d*1e3*out.SN_conc*1e-9*14/62/param.Dt # becasue we assume uniform snow accumulation, in kgN m-2 s-1 for each snow layer
    FA_mass[i] = FA[idx].sum()
    FA_d15N[i] = (FA[idx]*out.SN_d15N[idx]).sum()/FA[idx].sum()
    FA_D17O[i] = (FA[idx]*out.SN_D17O[idx]).sum()/FA[idx].sum()

frem_flux = FA_mass/out.FD_mass[0,:]

# resample into 5 cm interval and 2 m 
# res_SN = np.zeros(())


# weekly variables
# include FP, FD, MASS, d15N, D17O, surface eps, frac. total 3*2+2 = 8
frem_f = f.frem_frac(out.frem)
yr = len(frem_f)
res_weekly = np.zeros((11,yr))
res_weekly[0,:] = out.FP_mass[0,:]
res_weekly[1,:] = out.FP_d15N[0,:]
res_weekly[2,:] = out.FP_D17O[0,:]
res_weekly[3,:] = out.FD_mass[0,:]
res_weekly[4,:] = out.FD_d15N[0,:]
res_weekly[5,:] = out.FD_D17O[0,:]
res_weekly[6,:] = out.eps_surface[0,:]
res_weekly[7,:] = frem_f
res_weekly[8,:] = out.surface_conc[0,:]
res_weekly[9,:] = out.surface_d15N[0,:]
res_weekly[10,:] = out.surface_D17O[0,:]
###############################################################################################
# calculate d15N_total
row = 0
for i in range(len(new_ma[:,0])):
    if new_ma[i,4]==-1 and new_ma[i-1,4]!=-1:
        row = i
trans = new_ma[:row,:] # exclude the initial layers

d15N_total = (trans[:,5]*trans[:,1]*trans[:,2]).sum()/(trans[:,5]*trans[:,1]).sum()
print('total d15N= %s'%d15N_total)
D17O_total = (trans[:,5]*trans[:,1]*trans[:,3]).sum()/(trans[:,5]*trans[:,1]).sum()
print( 'total D17O= %s'%D17O_total)
d15N_conc = (trans[:,5]*trans[:,1]).sum()/trans[:,5].sum()
print('total conc= %s'%d15N_conc)
##############################################################################################
# np.savetxt('observe/result_all.txt',new_ma)  #save date

plt.plot(out.SN_conc, -out.mid_depth)
plt.show()
plt.plot(out.SN_d15N, -out.mid_depth)
plt.show()
plt.plot(out.SN_D17O, -out.mid_depth)
plt.show()

yrindex = []
for i in range(len(out.SN_date)):
    if state.SN_date[i]%52 == 51 and state.SN_date[i-1]%52 != 51:
        yrindex.append(i)
def weekly(new_ma):
    res = np.zeros((52,3))
    for week in range(52):
        conc = []
        d15N = []
        D17O = []
        thickness = []
        for layer in range(len(new_ma)):
            if new_ma[layer,4]%52 == week:
                conc.append(new_ma[layer,1])
                d15N.append(new_ma[layer,2])
                D17O.append(new_ma[layer,3])
                thickness.append(new_ma[layer,5])
        res[week,0] = (np.array(conc)*np.array(thickness)).sum()/np.array(thickness).sum()
        res[week,1] = (np.array(conc)*np.array(thickness)*np.array(d15N)).sum()/(np.array(conc)*np.array(thickness)).sum()
        res[week,2] = (np.array(conc)*np.array(thickness)*np.array(D17O)).sum()/(np.array(conc)*np.array(thickness)).sum()
    return res

fig, (ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)

for yr in range(len(yrindex)-1):        
    ax1.plot(weekly(new_ma[yrindex[yr]:yrindex[yr+1],:])[:,1],label = 'year %s'%(len(yrindex)-yr-1))
    ax2.plot(weekly(new_ma[yrindex[yr]:yrindex[yr+1],:])[:,2],label = 'year %s'%(len(yrindex)-yr-1))
    ax3.plot(weekly(new_ma[yrindex[yr]:yrindex[yr+1],:])[:,0],label = 'year %s'%(len(yrindex)-yr-1))
ax1.legend()
plt.show() 

# write the output data into .nc file
parameters = [param.out_fric,param.f_cage,param.Phi]
filename = 'ncoutput\%s.nc'%('Glacial_base_case_1980_oz_with_phiT')
f.create_ncfile(filename,new_ma,res_weekly,parameters)  














