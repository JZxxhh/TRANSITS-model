#################################################################################################
### IMPORT MODULES ##############################################################################
# We here import commonly used modules
import numpy as np                                          # array module                       
from pylab import *
from time import *                                          # module to measure time
import function as f
import operator


#################################################################################################
# List all O3col values used to run the TUV-model
list_O3 = list(range(25, 500, 25))
for i in [500, 750, 1000, 1500, 2000, 3000]:
    list_O3.append(i)

#################################################################################################
### Below is the model !

def model(p, opt, s, i, o):
    ### INITIALIZE the model
    # Option : if the model is run with a varying ozone column per year :
    if opt.O3col_variations == 1:
        # Read the ozone column data
        file_op = open("inputs_O3col\inputs_O3col_SP_1978-2005.txt","r")
        read_O3col_variations = np.genfromtxt(file_op, skip_header=0, delimiter="")
        # If number of years in the O3 column variations table is smaller than the number of computed years, warn the user :
        if p.nb_yr < shape(read_O3col_variations)[1]:
            print('!!! WARNING (option O3 col) !!! The number of computed years is smaller than the available number of years in the O3 col variation file')
        elif p.nb_yr > shape(read_O3col_variations)[1]:
            print('!!! WARNING (option O3 col) !!! The number of computed years is higher than the available number of years in the O3 col variation file')
            print('                                The last year in the file will be used for the extra computed years')

    ### For each time step
    for t in range(p.cp_dur):
        # we here first calculate the time step to look at in the original input file (that only represents one year)
        # example : when t=54, the model will look at tt=54%52=2 in the input file, i.e. the 3rd week
        tt = t%p.nb_ts

        # Measure the duration of one computed year
        if tt == 0:
            tps_deb = time()

        ### STEP 1 : SNOW BURIAL ################################################################
        # If sastrugi option is set, apply the defined thickness
        if operator.and_(opt.sastrugi_effect == 1, t == opt.sastrugi_TS + (p.nb_yr-1)*p.nb_ts):
            accu_thick = opt.sastrugi_thick
        else:
            # we first convert the accumulation into a snow layer thickness
            accu_thick = i.accu_mass[tt] / p.SN_rho
        if accu_thick > 0:
            # All layers are pushed lower
            dim = len(s.SN_thick)
            # Set thick at top to being equal to accumulation rate
            s.SN_thick[1:] = s.SN_thick[0:dim-1]
            s.SN_thick[0]  = accu_thick
            # Set the date to the time step number + 1 (since this time step is achieved)
            s.SN_date[1:, 0]  = s.SN_date[0:dim-1, 0]
            s.SN_date[0,  0]  = (t+1)/p.nb_ts*52*7
            # Set nitrate mass and isotopic composition at top to be null
            s.SN_mass[1:, 0]  = s.SN_mass[0:dim-1, 0]
            s.SN_mass[0,  0]  = 0
            s.SN_d15N[1:, 0]  = s.SN_d15N[0:dim-1, 0]
            s.SN_d15N[0,  0]  = 0
            s.SN_D17O[1:, 0]  = s.SN_D17O[0:dim-1, 0]
            s.SN_D17O[0,  0]  = 0
            # Set depth at top to be equal to accumulation rate
            s.SN_depth[1:] = s.SN_depth[0:dim-1]
            s.SN_depth[:]  = s.SN_depth[:] + accu_thick
            # Set JO3 at top to be null
            s.SN_JO3[1:] = s.SN_JO3[0:dim-1]
            s.SN_depth[0]  = 0
            # >>> At this stage, SN_J14 and SN_J15 do not need to be pushed down since they are calculated at the begining of STEP 2
        else:
            # if the accumulation thickness is null, we warn the user
            print("WARNING at time step %s : snow accumulation thickness is NULL!" % (tt))

        ### STEP 2 : DEPOSITION FLUX DISTRIBUTION ###############################################
        # Reset deposited mass arrays
        s.SN_massD       = np.zeros((len(s.SN_depth), 1))
        s.SN_fracD       = np.zeros((len(s.SN_depth), 1))
        # The deposited mass (kgN on computed surface) used here was calculated at the end of the previous time step
        # Now, this flux is deposited given a deposition/diffusion scenario
        if opt.dep_diff_scenario in [0, 2]:
            # 100% of the nitrate deposition occurs on the top layer
            s.SN_massD[0, 0]    = s.FD_mass * p.Dt
            # Dry deposited d15N :
            # The dry deposition does not induce any isotopic fractionation on N
            s.SN_d15ND       = s.FD_d15N
            # Dry deposited D17O :
            # The dry deposition does not induce any isotopic fractionation on O
            s.SN_D17OD       = s.FD_D17O
        elif opt.dep_diff_scenario == 1:
            # 100% of the nitrate deposition occurs on the top layer
            # BUT, we distribute the dry deposited nitrate flux at depth
            # For this, we calculate the effective diffusivity given atmospheric (air temp and pressure) and snow (SSA, tau, rho) properties
            Deff = f.DgeffHNO3(p.SN_SSA, p.SN_rho, p.SN_tau, i.air_temp[tt], i.air_pres[tt])
            # We prepare the calculation of the fraction deposited at each depth
            tau                   = 2 * (Deff * p.Dt)**0.5              # the parameter used in the erfc function
            int_erfc_tau_infinite = f.int_erfc_tau(1000000, tau, p.Pi)  # the infinite value of the integral of this erfc function
            # Deposited fractions calculation
            for j in range(len(s.SN_depth)):
                s.SN_fracD[j, 0] = (f.int_erfc_tau(s.SN_depth[j, 0], tau, p.Pi) - f.int_erfc_tau(s.SN_depth[j, 0] - s.SN_thick[j, 0], tau, p.Pi)) / int_erfc_tau_infinite        
            # Dry deposited mass in each snow layer (kgN.m^(-2)) :
            s.SN_massD[:]        = s.FD_mass * p.Dt * s.SN_fracD[:]
            # Dry deposited d15N :
            # The dry deposition does not induce any isotopic fractionation on N
            s.SN_d15ND       = s.FD_d15N
            # Dry deposited D17O :
            # The dry deposition does not induce any isotopic fractionation on O
            s.SN_D17OD       = s.FD_D17O

        # A mass and isotope balance IN SNOW is achieved after this deposition
        # The isotope mass balance gives
        # For N :
        s.SN_d15N[:] = (s.SN_d15N[:] * s.SN_mass[:] + s.SN_d15ND * s.SN_massD[:]) / (s.SN_mass[:] + s.SN_massD[:])
        # For O :
        s.SN_D17O[:] = (s.SN_D17O[:] * s.SN_mass[:] + s.SN_D17OD * s.SN_massD[:]) / (s.SN_mass[:] + s.SN_massD[:])
        # The mass balance in snow gives :
        s.SN_mass[:] = s.SN_mass[:] + s.SN_massD[:]
        # Replace NAN values (occuring in the case of divisions by 0) by 0.
        s.SN_d15N[np.isnan(s.SN_d15N)]=0.
        s.SN_D17O[np.isnan(s.SN_D17O)]=0.
        s.SN_mass[np.isnan(s.SN_mass)]=0.
        
        ### STEP 3 : DIFFUSION IN THE SNOWPACK    ###############################################
        # Resample snow profile
        # Calculate the resampling array
        RESAMP = f.resample_after_deposition(s.SN_depth[:], accu_thick, p.SN_set_thick)
        # Apply the exchange array (which has been calculated above)
        s.SN_thick[:]          = np.dot(RESAMP, s.SN_thick[:])
        s.SN_depth[:]          = np.dot(np.tril(np.ones((len(s.SN_depth[:]), len(s.SN_depth[:]))), 0), s.SN_thick[:])
        s.SN_d15N_mass[:]      = s.SN_d15N[:]*s.SN_mass[:]
        s.SN_D17O_mass[:]      = s.SN_D17O[:]*s.SN_mass[:]
        s.SN_date_mass[:]      = s.SN_date[:]*s.SN_mass[:]
        s.SN_JO3_mass[:]       = s.SN_JO3[:]*s.SN_mass[:]
        s.SN_mass[:]           = np.dot(RESAMP[:], s.SN_mass[:])
        s.SN_d15N_mass[:]      = np.dot(RESAMP[:], s.SN_d15N_mass[:])
        s.SN_D17O_mass[:]      = np.dot(RESAMP[:], s.SN_D17O_mass[:])
        s.SN_date_mass[:]      = np.dot(RESAMP[:], s.SN_date_mass[:])
        s.SN_JO3_mass[:]       = np.dot(RESAMP[:], s.SN_JO3_mass[:])
        s.SN_d15N[:]           = s.SN_d15N_mass[:] / s.SN_mass[:]
        s.SN_D17O[:]           = s.SN_D17O_mass[:] / s.SN_mass[:]
        s.SN_date[:]           = s.SN_date_mass[:] / s.SN_mass[:]
        s.SN_JO3[:]            = s.SN_JO3_mass[:] / s.SN_mass[:]
        # Replace NAN values (occuring in the case of divisions by 0) by 0. 
        s.SN_d15N[np.isnan(s.SN_d15N)]=0.
        s.SN_D17O[np.isnan(s.SN_D17O)]=0.
        s.SN_date[np.isnan(s.SN_date)]=0.
        s.SN_JO3[np.isnan(s.SN_JO3)]  =0.
        
        # record the depth and thick state
        o.depth_re[:,t]        = s.SN_depth[:,0]
        o.thick_re[:,t]        = s.SN_thick[:,0]

        # The thickness of the archived layer is equal to the snow accumulated mass
        o.FA_thick[0, t]                  = accu_thick
        # The archived mass, date, JO3 and isotopic composition are those of the lower layer
        o.FA_mass[0, t] = s.SN_mass[len(s.SN_mass)-1, 0] # kgN.m-2.TS-1 (per time step)
        o.FA_d15N[0, t] = s.SN_d15N[len(s.SN_mass)-1, 0]
        o.FA_D17O[0, t] = s.SN_D17O[len(s.SN_mass)-1, 0]
        o.FA_JO3[0, t]  = s.SN_JO3[len(s.SN_mass)-1, 0]
        o.FA_date[0, t] = s.SN_date[len(s.SN_mass)-1, 0]
        # Calculate concentration in the archived flux
        o.FA_conc[0, t] = o.FA_mass[0, t] * (p.M_NO3/p.M_N) * 1E12 / (1000 * o.FA_thick[0, t] * p.cp_S * p.SN_rho)
        
        # Diffusion occures only in the case of a specific deposition/diffusion scenario
        if opt.dep_diff_scenario == 2:
            # Calculate diffusion matrix
            # Apply to mass and isotopes in snow
            # Assume no fractionation
            # Assume no loss from snowpack to the atmosphere
            Dt_short = 4.*60*60 # seconds (4 hours)
            f.calc_diffusion(s, opt.Diff, p.Dt, Dt_short, p.SN_set_thick)

        ### STEP 4 : NITRATE PHOTOLYSIS IN THE SNOWPACK       ###################################
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PHOTOLYSIS (1)
        ### Option : if the model is run with a varying ozone column per year :
        if opt.O3col_variations == 1:
        # For each year, we have attributed the ozone column time series
            # At each beginning of year, read the O3 column data
            if tt == 0:
                # If the computed year is greater than the number of years in O3 col variation file
                if t/p.nb_ts > shape(read_O3col_variations)[1]-1:
                    # Use data from the last year
                    i.O3colDU[:] = read_O3col_variations[1:, shape(read_O3col_variations)[1]-1]
                    print('... use O3 column data from year ' + str(int(read_O3col_variations[0, shape(read_O3col_variations)[1]-1])))
                else:
                    i.O3colDU[:] = read_O3col_variations[1:, t/p.nb_ts]
                    print('... use O3 column data from year ' + str(int(read_O3col_variations[0, t/p.nb_ts])))
        # For each time step, associate the OZONE COLUMN value to the closest in O3col values computed with TUV
        # Check position in list_O3 and return closest value
        lc = list(sqrt(pow((array(list_O3) - i.O3colDU[tt]), 2)))
        valO3 = list_O3[f.indexMin(lc)]
        # Therefore [safety check], if the OZONE COLUMN value is null (means it's winter time), then the O3col value is set to smallest value in list-O3 (i.e. 25 DU) 
        # Now, determine which J14 and J15 arrays to use
        # read the appropriate TUV prepared file depending on the OZONE COLUMN value (valO3)
        if valO3 == 25:
            J14 = i.J14NO3.DU0025
            J15 = i.J15NO3.DU0025
            JO3 = i.JO3.DU0025
        elif valO3 == 50:
            J14 = i.J14NO3.DU0050
            J15 = i.J15NO3.DU0050
            JO3 = i.JO3.DU0050
        elif valO3 == 75:
            J14 = i.J14NO3.DU0075
            J15 = i.J15NO3.DU0075
            JO3 = i.JO3.DU0075
        elif valO3 == 100:
            J14 = i.J14NO3.DU0100
            J15 = i.J15NO3.DU0100
            JO3 = i.JO3.DU0100
            
        elif valO3 == 125:
            J14 = i.J14NO3.DU0125
            J15 = i.J15NO3.DU0125
            JO3 = i.JO3.DU0125
        elif valO3 == 150:
            J14 = i.J14NO3.DU0150
            J15 = i.J15NO3.DU0150
            JO3 = i.JO3.DU0150
        elif valO3 == 175:
            J14 = i.J14NO3.DU0175
            J15 = i.J15NO3.DU0175
            JO3 = i.JO3.DU0175
        elif valO3 == 200:
            J14 = i.J14NO3.DU0200
            J15 = i.J15NO3.DU0200
            JO3 = i.JO3.DU0200
        elif valO3 == 225:
            J14 = i.J14NO3.DU0225
            J15 = i.J15NO3.DU0225
            JO3 = i.JO3.DU0225
        elif valO3 == 250:
            J14 = i.J14NO3.DU0250
            J15 = i.J15NO3.DU0250
            JO3 = i.JO3.DU0250
        elif valO3 == 275:
            J14 = i.J14NO3.DU0275
            J15 = i.J15NO3.DU0275
            JO3 = i.JO3.DU0275
        elif valO3 == 300:
            J14 = i.J14NO3.DU0300
            J15 = i.J15NO3.DU0300
            JO3 = i.JO3.DU0300
        elif valO3 == 325:
            J14 = i.J14NO3.DU0325
            J15 = i.J15NO3.DU0325
            JO3 = i.JO3.DU0325
        elif valO3 == 350:
            J14 = i.J14NO3.DU0350
            J15 = i.J15NO3.DU0350
            JO3 = i.JO3.DU0350
        elif valO3 == 375:
            J14 = i.J14NO3.DU0375
            J15 = i.J15NO3.DU0375
            JO3 = i.JO3.DU0375
        elif valO3 == 400:
            J14 = i.J14NO3.DU0400
            J15 = i.J15NO3.DU0400
            JO3 = i.JO3.DU0400
        elif valO3 == 425:
            J14 = i.J14NO3.DU0425
            J15 = i.J15NO3.DU0425
            JO3 = i.JO3.DU0425
        elif valO3 == 450:
            J14 = i.J14NO3.DU0450
            J15 = i.J15NO3.DU0450
            JO3 = i.JO3.DU0450
        elif valO3 == 475:
            J14 = i.J14NO3.DU0475
            J15 = i.J15NO3.DU0475
            JO3 = i.JO3.DU0475
        elif valO3 == 500:
            J14 = i.J14NO3.DU0500
            J15 = i.J15NO3.DU0500
            JO3 = i.JO3.DU0500
        elif valO3 == 750:
            J14 = i.J14NO3.DU0750
            J15 = i.J15NO3.DU0750
            JO3 = i.JO3.DU0750
        elif valO3 == 1000:
            J14 = i.J14NO3.DU1000
            J15 = i.J15NO3.DU1000
            JO3 = i.JO3.DU1000
        elif valO3 == 1500:
            J14 = i.J14NO3.DU1500
            J15 = i.J15NO3.DU1500
            JO3 = i.JO3.DU1500
        elif valO3 == 2000:
            J14 = i.J14NO3.DU2000
            J15 = i.J15NO3.DU2000
            JO3 = i.JO3.DU2000
        elif valO3 == 3000:
            J14 = i.J14NO3.DU3000
            J15 = i.J15NO3.DU3000
            JO3 = i.JO3.DU3000
        # List depths where J data are available (1D array)
        # In general, depths range from 0 (surface) to 1m depth
        J_depth = i.J14NO3.header[1:]*p.photic_zone
         
        # The SZA repartition file has been previously read and converted from minutes to seconds.
        # Depending on the time step, we extract the row with the repartition in the SZA's (resolution is 1 degree).
        # we eventually calculate the J depth profile at this specific time step by multiplying the J14NO3 (and J15NO3, in seconds) vs depth vs SZA files by the SZA_rep (in seconds)
        # the obtained J14_TS and J14_TS which are expressed in photons hitting the considered nitrate ions during the time step.	
        # indeed : J14NO3_TS = integral of Actinic_flux x Xsection on wavelengths 280-350 x Dt versus depth
        J14_TS = np.zeros((len(J_depth), 1))
        J15_TS = np.zeros((len(J_depth), 1))
        J14_TS[:, 0] = (J14[1:,:] * i.SZA_rep[tt, :]).sum(axis=1)   # in seconds, Phi=1 ,add the rows into an array
        J15_TS[:, 0] = (J15[1:,:] * i.SZA_rep[tt, :]).sum(axis=1)   # in seconds, Phi=1
        
        # J14 and J15 are linearily interpolated at depth depending on the snowpack layer resolution
        # For depths greater than 1 meter : J = 0.
        J14_TS_inter = np.zeros((len(s.SN_depth), 1))
        J15_TS_inter = np.zeros((len(s.SN_depth), 1))
        J14_TS_inter[:, 0] = interp(s.SN_depth[0:len(s.SN_depth), 0]-0.5*s.SN_thick[0:len(s.SN_depth), 0], J_depth[:], p.actinic_flux*J14_TS[:, 0], left=0, right=0)
        J15_TS_inter[:, 0] = interp(s.SN_depth[0:len(s.SN_depth), 0]-0.5*s.SN_thick[0:len(s.SN_depth), 0], J_depth[:], p.actinic_flux*J15_TS[:, 0], left=0, right=0)

        ### J14 calculation (expressed in number of hit nitrate ions) with the desired Phi
        s.SN_J14[:, 0] = p.Phi * J14_TS_inter[:, 0]
        if (exp(-s.SN_J14)).max() > 1:
            a = str(int(t))
            print("WARNING at time step %s : the J*Dt term gets greater than 1... The model photolyses more nitrate than available!" % (a))

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PHOTOLYSIS (2)
        # we now calculate the remaining fraction in each snow layer
        frem = exp(-s.SN_J14[:])
        o.frem[:,t]=frem[:,0]
        # remaining fraction as a function of depth
        # we now calculate the remaining and emitted nitrate mass fractions
        s.SN_massR[:] = frem[:] * s.SN_mass[:]                  # the remaining fraction in the snow
        s.SN_massE[:] = (1 - frem[:]) * s.SN_mass[:]            # the NOx produced at the vicinity of the snow grain (not emitted yet)
        # PHOTOLYSIS : we calculate the d15N in both fractions
        # first, we calculate the 15eps at this time step and at depth :
        s.SN_eps15[:] = ((J15_TS_inter[:]/J14_TS_inter[:]) - 1) * 1000
        # Replace NAN values (occuring in the case of divisions by 0) by 0. 
        s.SN_eps15[np.isnan(s.SN_eps15)]=0.
        # This corresponds to the Rayleigh distillation of an opened system (Rayleigh, 1902 in Morin, S., 2008)
        s.SN_d15NR[:] = ((1 + s.SN_d15N[:]/1000) * pow(frem[:], (s.SN_eps15[:]/1000)) - 1) * 1000
        s.SN_d15NE[:] = ((1 + s.SN_d15N[:]/1000) * (1 - pow(frem[:], s.SN_eps15[:]/1000 + 1))/(1 - frem[:]) - 1) * 1000
        # Replace NAN and INFinite values (occuring in the case of divisions by 0) by 0.
        s.SN_d15NR[np.isnan(s.SN_d15NR)]=0.
        s.SN_d15NE[np.isnan(s.SN_d15NE)]=0.
        s.SN_d15NR[np.isinf(s.SN_d15NR)]=0.
        s.SN_d15NE[np.isinf(s.SN_d15NE)]=0.
        # PHOTOLYSIS : we calculate the D17O in both fractions
        # Photolysis does not induce any isotopic fractionation in O otams so that :
        s.SN_D17OR[:] = s.SN_D17O[:]
        s.SN_D17OE[:] = s.SN_D17O[:]

        ### STEP 5 : CAGE EFFECT #########################
        # a small fraction of the emitted fraction (which is still in the vicinity of the snow grain) undergoes isotopic exchange with water in the surrounding
        # We first calculate which mass fraction undergoes cage effect :
        s.SN_massC[:] = p.f_cage * s.SN_massE[:]
        # The emitted mass fraction is therefore smaller :
        s.SN_massE[:] = (1 - p.f_cage) * s.SN_massE[:]
        # We calculate the d15N in this new fraction
        # Nitrogen isotopes are preserved by this process so that :
        s.SN_d15NC[:] = s.SN_d15NE[:]
        # We calculate the D17O in this new fraction
        # Oxygen isotopes are not preserved by this process. Each molecule undergoing cage effect isotopically
        # exchanges one oxygen atom with water (which features D17O_water), so that :
        s.SN_D17OC[:] = (2 * s.SN_D17OE[:] + 1 * p.D17O_water) / 3

        ### STEP 6 : MASS and ISOTOPE BALANCE IN SNOW ##################################"
        # The mass balance in snow gives :        
        s.SN_mass[:] = s.SN_massR[:] + s.SN_massC[:]
        # ISOTOPIC MASS BALANCE IN SNOW
        # For nitrogen isotopes :
        s.SN_d15N[:] = (s.SN_d15NR[:]*s.SN_massR[:] + s.SN_d15NC[:]*s.SN_massC[:]) / s.SN_mass[:]
        # For oxygen isotopes :
        s.SN_D17O[:] = (s.SN_D17OR[:]*s.SN_massR[:] + s.SN_D17OC[:]*s.SN_massC[:]) / s.SN_mass[:]
        # Replace NAN and INFinite values (occuring in the case of divisions by 0) by 0. 
        s.SN_d15N[np.isnan(s.SN_d15N)]=0.
        s.SN_D17O[np.isnan(s.SN_D17O)]=0.
        s.SN_d15N[np.isinf(s.SN_d15N)]=0.
        s.SN_D17O[np.isinf(s.SN_D17O)]=0.

        # Store values of NOx emitted in each layer
        if t > (p.cp_dur - p.nb_ts):
            o.SN_FP_mass_contri[:, t%p.nb_ts] = s.SN_massE[0:len(s.SN_depth)-1, 0] / p.Dt     # NOx mass flux in kgN.m-2.s-1
            o.SN_FP_d15N_contri[:, t%p.nb_ts] = s.SN_d15NE[0:len(s.SN_depth)-1, 0]     # d15N in flux, permil
            o.SN_FP_D17O_contri[:, t%p.nb_ts] = s.SN_D17OE[0:len(s.SN_depth)-1, 0]     # D17Oin flux, permil
        
        ### STEP 7 : NITRATE MASS AND ISOTOPE BALANCE IN THE ATMOSPHERE #########################
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PHOTOLYTIC FLUX FP
        # we calculate the mass of this flux
        s.FP_mass = s.SN_massE.sum(dtype=float) / p.Dt        # the flux is expressed in kgN.m^(-2).s^(-1)
        # PHOTOLYTIC FLUX FP : we calculate the d15N and D17O of this flux
        if s.FP_mass == 0.:
            s.FP_d15N = 0.
            s.FP_D17O = 0.
        else:
            s.FP_d15N = (s.SN_d15NE[:] * s.SN_massE[:] / p.Dt).sum(dtype=float) / s.FP_mass
            s.FP_D17O = (s.SN_D17OE[:] * s.SN_massE[:] / p.Dt).sum(dtype=float) / s.FP_mass
        # PHOTOLYTIC FLUX FP : we here store the data prior to the local oxidation calculation to track the isotopic composition of the emitted NOx
        if t > (p.cp_dur - p.nb_ts):
            o.FP_d15N_NOx[0, t%p.nb_ts] = s.FP_d15N
            o.FP_D17O_NOx[0, t%p.nb_ts] = s.FP_D17O

        # MASS BALANCE IN THE ATMOSPHERE
        # Calculate the deposited flux FD (will be deposited at next time step)
        s.FD_mass = (i.AT_mass[tt-1] - i.AT_mass[tt]) / p.Dt + (1 - i.f_exp[tt]) * (s.FP_mass + i.FS_mass[tt] + i.FT_mass[tt])  # kgN.m-2.s-1
        s.FE_mass = i.f_exp[tt] * (s.FP_mass + i.FS_mass[tt] + i.FT_mass[tt])                                                   # kgN.m-2.s-1

        # Mass balance in atmosphere
        if s.FP_mass == 0:
            s.AT_d15N = (i.AT_mass[tt-1] * s.AT_d15N + p.Dt * (1 - i.f_exp[tt]) * (
                i.FS_d15N[tt] * i.FS_mass[tt] + i.FT_d15N[tt] * i.FT_mass[tt]) - p.Dt * s.FD_mass * p.AT_eps15_dep
                         ) / (i.AT_mass[tt] + p.Dt * s.FD_mass)
        else:
            s.AT_d15N = (i.AT_mass[tt-1] * s.AT_d15N + p.Dt * (1 - i.f_exp[tt]) * (
                s.FP_d15N * s.FP_mass + i.FS_d15N[tt] * i.FS_mass[tt] + i.FT_d15N[tt] * i.FT_mass[tt]) - p.Dt * s.FD_mass *  p.AT_eps15_dep
                         ) / (i.AT_mass[tt] + p.Dt * s.FD_mass)

        # Isotopic composition of FD 
        s.FD_d15N = s.AT_d15N + p.AT_eps15_dep    # same d15N value as in the atmosphere but fractionation at deposition
        s.FD_D17O = s.AT_D17O                     # same D17O value as in the atmosphere but NO fractionation at deposition

        # Isotopic composition of FE
        s.FE_d15N = (s.FP_d15N * s.FP_mass + i.FS_d15N[tt] * i.FS_mass[tt] + i.FT_d15N[tt] * i.FT_mass[tt])/(
                s.FP_mass + i.FS_mass[tt] + i.FT_mass[tt])
        s.FE_D17O = (s.FP_D17O * s.FP_mass + i.FS_D17O[tt] * i.FS_mass[tt] + i.FT_D17O[tt] * i.FT_mass[tt])/(
                s.FP_mass +i.FS_mass[tt] + i.FT_mass[tt])
        
        # Local oxidation of the photolytic NOx flux
        # We calculate the photochemical steady-state (PSS) to obtain D17O of NOx
        if opt.NO2_oxidation == 1:
            s.AT_D17O_NO2_PSS = i.alpha[tt] * (1.18 * i.D17O_O3[tt] + 6.6)
            # We calculate below the mean based on the SZA repartition and the JO3 (SZA dependent) measured with the current ozone column in time step
            # J_O3 is taken at the surface of the snow pack
            JO3_mean = (JO3 * i.SZA_rep[tt, :]).sum() / p.Dt
            # J_O3 is stored for the new layer (archived parameter. Allows to track the seasonnality at depth)
            s.SN_JO3[0, 0] = JO3_mean
            # We here calculate the D17O in the additionnal O atom
            if JO3_mean > p.JO3_threshold:
                s.AT_D17O_addO = i.x_fact[tt] * 0.5 * 1.5 *i.D17O_O3[tt]
            else:
                s.AT_D17O_addO = 1.5 * i.D17O_O3[tt] + JO3_mean/p.JO3_threshold * (i.x_fact[tt] * 0.5 * 1.5 *i.D17O_O3[tt] - 1.5 * i.D17O_O3[tt])
            # Only if the FP flux is non-zero
            if s.FP_mass == 0.:
                s.FP_D17O = 0.
            else:
                s.FP_D17O = (2 * s.AT_D17O_NO2_PSS + 1 * s.AT_D17O_addO) / 3
        # Only if the FP flux is non-zero
        else:
            s.AT_D17O_NO2_PSS = 0.
            if s.FP_mass == 0.:
                s.FP_D17O = 0.                
        # The mass balance gives the expression of D17O in atmospheric nitrate
        # If the FP flux is zero, we remove the photolytic flux in the budget
        if s.FP_mass == 0.:
            s.AT_D17O = (s.AT_mass * s.AT_D17O + p.Dt * (1 - i.f_exp[tt]) * (i.FS_D17O[tt] * i.FS_mass[tt] + i.FT_D17O[tt] * i.FT_mass[tt])) / (s.AT_mass + p.Dt * s.FD_mass)
        else:
            s.AT_D17O = (s.AT_mass * s.AT_D17O + p.Dt * (1 - i.f_exp[tt]) * (s.FP_D17O * s.FP_mass + i.FS_D17O[tt] * i.FS_mass[tt] + i.FT_D17O[tt] * i.FT_mass[tt]))/ (s.AT_mass + p.Dt * s.FD_mass)

        ### STEP 8 : OUTPUTS EXPORT AND CALCULATION #############################################
        # Exportation of each column array to the (n x t) arrays
        if t >= (p.cp_dur - p.nb_ts):
            o.SN_depth[:, t%p.nb_ts]          = s.SN_depth[0:len(s.SN_depth)-1, 0]
            o.SN_thick[:, t%p.nb_ts]          = s.SN_thick[0:len(s.SN_depth)-1, 0]
            o.SN_mass[:, t%p.nb_ts]           = s.SN_mass[0:len(s.SN_depth)-1, 0]
            o.SN_date[:, t%p.nb_ts]           = s.SN_date[0:len(s.SN_depth)-1, 0]
            o.SN_JO3[:, t%p.nb_ts]            = s.SN_JO3[0:len(s.SN_depth)-1, 0] / p.Dt
            o.SN_J14[:, t%p.nb_ts]            = s.SN_J14[0:len(s.SN_depth)-1, 0] / p.Dt
            o.SN_eps15[:, t%p.nb_ts]          = s.SN_eps15[0:len(s.SN_depth)-1, 0]
            o.SN_d15N[:, t%p.nb_ts]           = s.SN_d15N[0:len(s.SN_depth)-1, 0]
            o.SN_D17O[:, t%p.nb_ts]           = s.SN_D17O[0:len(s.SN_depth)-1, 0]
            o.FD_mass[0, t%p.nb_ts-1]         = s.FD_mass
            o.FD_d15N[0, t%p.nb_ts-1]         = s.FD_d15N
            o.FD_D17O[0, t%p.nb_ts-1]         = s.FD_D17O
            o.FE_mass[0, t%p.nb_ts-1]         = s.FE_mass
            o.FE_d15N[0, t%p.nb_ts-1]         = s.FE_d15N
            o.FE_D17O[0, t%p.nb_ts-1]         = s.FE_D17O 
            o.FP_mass[0, t%p.nb_ts]           = s.FP_mass
            o.FP_d15N[0, t%p.nb_ts]           = s.FP_d15N
            o.FP_D17O[0, t%p.nb_ts]           = s.FP_D17O
            o.AT_mass[0, t%p.nb_ts]           = s.AT_mass
            o.AT_d15N[0, t%p.nb_ts]           = s.AT_d15N
            o.AT_D17O[0, t%p.nb_ts]           = s.AT_D17O
            o.AT_D17O_addO[0, t%p.nb_ts]      = s.AT_D17O_addO
            o.AT_D17O_NO2_PSS[0, t%p.nb_ts]   = s.AT_D17O_NO2_PSS

        o.FDall_mass[0, t]         = s.FD_mass
        o.FDall_d15N[0, t]         = s.FD_d15N
        o.FDall_D17O[0, t]         = s.FD_D17O
        o.FEall_mass[0, t]         = s.FE_mass
        o.FEall_d15N[0, t]         = s.FE_d15N
        o.FEall_D17O[0, t]         = s.FE_D17O
        o.FPall_mass[0, t]         = s.FP_mass
        o.FPall_d15N[0, t]         = s.FP_d15N
        o.FPall_D17O[0, t]         = s.FP_D17O

        # [Diagnostic tool] Calculate and store atmospheric nitrate mass along the whole computation
        ### at each time step
        o.AT_mass_int_TS[0, t]              = s.AT_mass
        ### for each computed year
        if tt == 0:
            count_AT_sum = s.AT_mass
        elif tt == p.nb_ts - 1:
            o.AT_mass_int_YR[0, int(t/p.nb_ts)] = count_AT_sum + s.AT_mass
        else:
            count_AT_sum = count_AT_sum + s.AT_mass

        # [Diagnostic tool] Calculate and store nitrate mass in snow pack
        ### at each time step
        o.SN_mass_int_TS[0, t] = sum(s.SN_mass[0:len(s.SN_depth)-1, 0])
        ### for each computed year
        if tt == 0:
            count_SN_sum = sum(s.SN_mass[0:len(s.SN_depth)-1, 0])
        elif tt == p.nb_ts - 1:
            o.SN_mass_int_YR[0, int(t/p.nb_ts)] = count_SN_sum + sum(s.SN_mass[0:len(s.SN_depth)-1, 0])
        else:
            count_SN_sum = count_SN_sum + sum(s.SN_mass[0:len(s.SN_depth)-1, 0])


        ### Returns computed duration of the computed year
        if tt == p.nb_ts-1:
            tps_fin = time()
            print('Year ' + str(int(t/p.nb_ts+1)) + ' computed in ' + str(str(int((tps_fin - tps_deb)*10)/10)) + ' s')      
### This is the end of the computational loop
#################################################################################################

