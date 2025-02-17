#################################################################################################
### IMPORT MODULES ##############################################################################
# We here import commonly used modules
import numpy as np                                          # array module                     
import time
import functions as f
from pylab import interp,exp


#################################################################################################
# List all O3col values used to run the TUV-model
list_O3 = range(200, 425, 25)

#################################################################################################
### Below is the model !
def model(p, opt, s, i, o):
    ### For each time step
    for t in range(p.cp_dur):
        tt = t%p.nb_ts
        if tt == 0:
            tps_deb = time.time() #start time
        ### STEP 1 : SNOW BURIAL ################################################################
        accu_thick = i.accu_mass[tt] / p.SN_rho   # in m       
        if accu_thick > 0 :
            s.SN_thick  = np.append(accu_thick,s.SN_thick)
            # Set the date to the time step number + 1 (since this time step is achieved)
            s.SN_date  = np.append(t,s.SN_date)
            # Set nitrate mass and isotopic composition at top to be null
            s.SN_mass  = np.append(0,s.SN_mass)
            #s.SN_mass[0,  0]  = 200.0E-9 #uncomment after fresh snowfall
            s.SN_d15N  = np.append(0,s.SN_d15N)
            #s.SN_d15N[0,  0]  = 7.6 #uncomment after fresh snowfall
            s.SN_D17O   = np.append(0,s.SN_D17O)
            #s.SN_D17O[0,  0]  = 28.
            # Set depth at top to be equal to accumulation rate 
            s.SN_depth =  np.append(0,s.SN_depth)
            s.SN_depth[:]  = s.SN_depth[:] + accu_thick
        else:
            # if the accumulation thickness is null, we warn the user
            print( "WARNING during year %s, at time step %s : snow accumulation thickness is negitive!" % (int(t/p.nb_ts)+1, tt))

        ### STEP 2 : DEPOSITION FLUX DISTRIBUTION ###############################################
        # Reset deposited mass arrays
        s.SN_massD       = np.zeros((len(s.SN_depth)))
        s.SN_fracD       = np.zeros((len(s.SN_depth)))
        # The deposited mass (kgN on computed surface) used here was calculated at the end of the previous time step
        # Now, this flux is deposited given a deposition/diffusion scenar
        if opt.dep_diff_scenario in [0, 2, 3, 4]:
            # 100% of the nitrate deposition occurs on the top layer
            s.SN_massD[0]    = s.FD_mass * p.Dt
            # Dry deposited d15N :
            # The dry deposition does not induce any isotopic fractionation on N
            s.SN_d15ND       = s.FD_d15N
            # Dry deposited D17O :
            # The dry deposition does not induce any isotopic fractionation on 0
            s.SN_D17OD       = s.FD_D17O
        elif opt.dep_diff_scenario == 1:
            # 100% of the nitrate deposition occurs on the top layer
            # BUT, we distribute the dry deposited nitrate flux at depth
            # For this, we calculate the effective diffusivity given atmospheric (air temp and pressure) and snow (SSA, tau, rho) properties
#            Deff = f.DgeffHNO3(p.SN_SSA, p.SN_rho, p.SN_tau, i.AT_Temp[tt], i.AT_Pres[tt])
            Deff = opt.Diff
            # We prepare the calculation of the fraction deposited at each depth
            tau                   = 2 * (Deff * p.Dt)**0.5              # the parameter used in the erfc function
            int_erfc_tau_infinite = f.int_erfc_tau(1000000, tau, p.Pi)  # the infinite value of the integral of this erfc function
            # Deposited fractions calculation
            for j in range(len(s.SN_depth)):
                s.SN_fracD[j] = (f.int_erfc_tau(s.SN_depth[j], tau, p.Pi) - f.int_erfc_tau(s.SN_depth[j] - s.SN_thick[j], tau, p.Pi)) / int_erfc_tau_infinite        
            # Dry deposited mass in each snow layer (kgN.m^(-2)) 
            s.SN_massD[:]        = s.FD_mass * p.Dt * s.SN_fracD[:]
            # Dry deposited d15N :
            # The dry deposition does not induce any isotopic fractionation on N
            s.SN_d15ND       = s.FD_d15N
            # Dry deposited D17O :
            # The dry deposition does not induce any isotopic fractionation on 0
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
        s.SN_thick = f.my_resample(s.SN_thick,accu_thick, p.SN_set_thick,'thickness')
        s.SN_depth = f.my_resample(s.SN_depth,accu_thick, p.SN_set_thick,'depth')
        s.SN_d15N = f.my_resample(s.SN_d15N,accu_thick, p.SN_set_thick,'d15N')
        s.SN_D17O = f.my_resample(s.SN_D17O,accu_thick, p.SN_set_thick,'D17O')
        s.SN_date = f.my_resample(s.SN_date,accu_thick, p.SN_set_thick,'date')
        s.SN_mass = f.my_resample(s.SN_mass,accu_thick, p.SN_set_thick,'mass')
        s.SN_year = f.my_resample(s.SN_year,accu_thick, p.SN_set_thick,'year')
        
        o.surface_conc[0,t] = (p.M_NO3 / p.M_N) * s.FD_mass * p.Dt / (accu_thick * p.cp_S * p.SN_rho) * pow(10, 9)
        o.surface_d15N[0,t] = s.FD_d15N
        o.surface_D17O[0,t] = s.FD_D17O
        ### STEP 4 : NITRATE PHOTOLYSIS IN THE SNOWPACK       ###################################
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PHOTOLYSIS (1)
        idx = np.abs(list_O3 - i.O3colDU[tt]).argmin()
        valO3 = list_O3[idx]
            # Therefore [safety check], if the OZONE COLUMN value is null (means it's winter time), then the O3col value is set to smallest value in list-O3 (i.e. 25 DU) 
        # Now, determine which J14 and J15 arrays to use
        # read the appropriate TUV prepared file depending on the OZONE COLUMN value (valO3)
        if valO3 == 150:
            J14  = i.J14NO3.DU150
            J15  = i.J15NO3.DU150
            JO3  = i.JO3.DU150
            JNO2 = i.JNO2.DU150
        elif valO3 == 175:
            J14  = i.J14NO3.DU175
            J15  = i.J15NO3.DU175
            JO3  = i.JO3.DU175
            JNO2 = i.JNO2.DU175
        if valO3 == 200:
            J14  = i.J14NO3.DU200
            J15  = i.J15NO3.DU200
            JO3  = i.JO3.DU200
            JNO2 = i.JNO2.DU200
        elif valO3 == 225:
            J14  = i.J14NO3.DU225
            J15  = i.J15NO3.DU225
            JO3  = i.JO3.DU225
            JNO2 = i.JNO2.DU225
        elif valO3 == 250:
            J14  = i.J14NO3.DU250
            J15  = i.J15NO3.DU250
            JO3 = i.JO3.DU250
            JNO2 = i.JNO2.DU250
        elif valO3 == 275:
            J14  = i.J14NO3.DU275
            J15  = i.J15NO3.DU275
            JO3  = i.JO3.DU275
            JNO2 = i.JNO2.DU275
        elif valO3 == 300:
            J14  = i.J14NO3.DU300
            J15  = i.J15NO3.DU300
            JO3  = i.JO3.DU300
            JNO2 = i.JNO2.DU300
        elif valO3 == 325:
            J14  = i.J14NO3.DU325
            J15  = i.J15NO3.DU325
            JO3  = i.JO3.DU325
            JNO2 = i.JNO2.DU325
        elif valO3 == 350:
            J14  = i.J14NO3.DU350
            J15  = i.J15NO3.DU350
            JO3  = i.JO3.DU350
            JNO2 = i.JNO2.DU350
        elif valO3 == 375:
            J14  = i.J14NO3.DU375
            J15  = i.J15NO3.DU375
            JO3  = i.JO3.DU375
            JNO2 = i.JNO2.DU375
        elif valO3 == 400:
            J14  = i.J14NO3.DU400
            J15  = i.J15NO3.DU400
            JO3  = i.JO3.DU400
            JNO2 = i.JNO2.DU400
        # List depths where J data are available (1D array)
        # In general, depths range from 0 (surface) to 1m depth
        J_depth = i.J14NO3.header[1:]/p.photic_zone
        # Give JNO2 value
        # The SZA repartition file has been previously read and converted from minutes to seconds.
        # Depending on the time step, we extract the row with the repartition in the SZA's (resolution is 1 degree).
        # we eventually calculate the J depth profile at this specific time step by multiplying the J14NO3 (and J15NO3, in seconds) vs depth vs SZA files by the SZA_rep (in seconds)
        # the obtained J14_TS and J14_TS which are expressed in photons hitting the considered nitrate ions during the time step.	
        # indeed : J14NO3_TS = integral of Actinic_flux x Xsection on wavelengths 280-350 x Dt versus depth
        J14_TS = np.zeros((len(J_depth)))
        J15_TS = np.zeros((len(J_depth)))
        J14_TS[:] = (J14[1:, :] * i.SZA_rep[tt, :]).sum(axis=1)   # in seconds, Phi=1
        J15_TS[:] = (J15[1:, :] * i.SZA_rep[tt, :]).sum(axis=1)   # in seconds, Phi=1
        
        # J14 and J15 are linearily interpolated at depth depending on the snowpack layer resolution
        # For depths greater than 1 meter : J = 0.
        J14_TS_inter = np.zeros((len(s.SN_depth)))
        J15_TS_inter = np.zeros((len(s.SN_depth)))
        J14_TS_inter[:] = interp(s.SN_depth[0:len(s.SN_depth)]-0.5*s.SN_thick[0:len(s.SN_depth)], J_depth[:], p.actinic_flux*J14_TS[:])
        J15_TS_inter[:] = interp(s.SN_depth[0:len(s.SN_depth)]-0.5*s.SN_thick[0:len(s.SN_depth)], J_depth[:], p.actinic_flux*J15_TS[:])
        J14_TS_inter[np.where(J14_TS_inter/J14_TS_inter[0] < 0.001)] = 0
        J15_TS_inter[np.where(J15_TS_inter/J15_TS_inter[0] < 0.001)] = 0
        

        ### J14 calculation (expressed in number of hit nitrate ions) with the desired Phi
        if opt.phi == 0:        
            p.Phi = np.exp(-2400/i.AT_Temp[tt]+3.6)    # USE equations from chu et al, 2003 for quantum yield
        Jdt = np.zeros_like(s.SN_depth)
        #p.Phi = np.exp(-2400/i.AT_Temp[tt]+3.6)
        #print p.Phi
        Jdt = p.Phi * J14_TS_inter[:] #sigema J*dt
#        o.J_NO3[0,t] = Jdt/i.SZA_rep[tt, :].sum()
        if (exp(-Jdt)).max() > 1:
            a = str(int(t))
            print( "WARNING at time step %s : the J*Dt term gets greater than 1... The model photolyses more nitrate than available!" % (a))

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PHOTOLYSIS (2)
        # we now calculate the remaining fraction in each snow layer
        frem = exp(-Jdt)
        # print('%s    %.3f'%(J14_TS_inter[0]*p.Phi/(7*24*3600), frem[0]))
        # remaining fraction as a function of depth
        # we now calculate the remaining and emitted nitrate mass fractions
        s.SN_massR = np.zeros_like(s.SN_depth) #created the array seam dimension as depth
        s.SN_massE = np.zeros_like(s.SN_depth)
        s.SN_eps15 = np.zeros_like(s.SN_depth)
        s.SN_d15NR = np.zeros_like(s.SN_depth)
        s.SN_d15NE = np.zeros_like(s.SN_depth)
        s.SN_D17OR = np.zeros_like(s.SN_depth)
        s.SN_D17OE = np.zeros_like(s.SN_depth)
        
        s.SN_massR[:] = frem[:] * s.SN_mass[:]                  # the remaining fraction in the snow
        s.SN_massE[:] = (1 - frem[:]) * s.SN_mass[:]            # the NO2 produced at the vicinity of the snow grain (not emitted yet)
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
        # Photolysis does not induce any isotopic fractionation so that :
        s.SN_D17OR[:] = s.SN_D17O[:]
        s.SN_D17OE[:] = s.SN_D17O[:]

        ### STEP 5 : CAGE EFFECT #########################
        # a small fraction of the emitted fraction (which is still in the vicinity of the snow grain) undergoes isotopic exchange with water in the surrounding
        # We first calculate which mass fraction undergoes cage effect :
        s.SN_massC=np.zeros_like(s.SN_depth)
        s.SN_d15NC=np.zeros_like(s.SN_depth)
        s.SN_D17OC=np.zeros_like(s.SN_depth)
            
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
        s.SN_d15N[np.isnan(s.SN_d15N)] = 0.
        s.SN_D17O[np.isnan(s.SN_D17O)] = 0.
        s.SN_d15N[np.isinf(s.SN_d15N)] = 0.
        s.SN_D17O[np.isinf(s.SN_D17O)] = 0.

        ### STEP 7 : PHOTOLYTIC FLUX BUDGET #####################################################
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
        # PHOTOLYTIC FLUX FP : we here store the data prior to the local oxidation calculation to track the isotopic composition of the emitted NO2
        if t > (p.cp_dur - p.nb_ts):
            o.FP_d15N_NO2[0, t%p.nb_ts] = s.FP_d15N
            o.FP_D17O_NO2[0, t%p.nb_ts] = s.FP_D17O
        ### STEP 8 : NO2 CYCLING AND LOCAL OXIDATION ############################################
        # We calculate below the mean J_NO2 value based on the SZA repartition and the JNO2 (SZA dependent) measured with the current ozone column in time step
        # J_NO2 is calculated at the surface of the snow pack
        #assuming PSS
        JNO2_mean = (JNO2 * i.SZA_rep[tt, :]).sum() / p.Dt
        # J_NO2 is stored for the new layer (archived parameter. Allows to track the seasonnality at depth)
        # J_NO2 is also stored in the atmosphere data
        s.AT_JNO2 = JNO2_mean
        # Estimate [OH] = OH_conc based on J_NO2 and relationship given in Kukui et al. (2014), figure 3a
        s.AT_OH_conc  = 2.46E9 * s.AT_JNO2     # molecule cm-3
#        print 'OH conc : ',s.AT_OH_conc
        # Estimate [RO2] = RO2_conc based on J_NO2 and relationship given in Kukui et al. (2014), figure 3b
        s.AT_RO2_conc = 93.8 * 1E9 * s.AT_JNO2  # molecule cm-3
        # Estimate [HO2] = HO2_conc based on [HO2]/[RO2] ratio given in Kukui et al. (2014)
        ratio_HO2_RO2 = 2/3
        s.AT_HO2_conc = ratio_HO2_RO2 * s.AT_RO2_conc  # molecule cm-3
        if p.AT_HO2_conc_scale != 1:
            s.AT_HO2_conc = p.AT_HO2_conc_scale * s.AT_HO2_conc
        # Estimate [CH3O2] = CH3O2_conc based on [HO2]/[RO2] ratio given in Kukui et al. (2014) and [RO2] = [HO2] + [CH3O2]
        s.AT_CH3O2_conc = (1 - ratio_HO2_RO2) * s.AT_RO2_conc  # molecule cm-3
        if p.AT_CH3O2_conc_scale != 1:
            s.AT_CH3O2_conc = p.AT_CH3O2_conc_scale * s.AT_CH3O2_conc
        # Calculate O3 mixing ratios given data in input file
        air_conc = (i.AT_Pres[tt] * 1E2 / (p.R * i.AT_Temp[tt])) * 1E-6 * p.N_A
        s.AT_O3_conc  = i.AT_O3_conc[tt] * air_conc * 1E-9      # i.AT_O3_conc[tt] is in ppbv
        s.AT_BrO_conc = p.AT_BrO_conc * air_conc                # p.AT_BrO_conc is in mol/mol
        # Calculate alpha value given mixing ratios and kinetic rates
        s.AT_alpha = (f.O3_NO_195_308K(i.AT_Temp[tt]) * s.AT_O3_conc +
                      f.BrO_NO_220_430K(i.AT_Temp[tt]) * s.AT_BrO_conc
                      ) / (f.O3_NO_195_308K(i.AT_Temp[tt]) * s.AT_O3_conc +
                           f.BrO_NO_220_430K(i.AT_Temp[tt]) * s.AT_BrO_conc +
                           f.HO2_NO_200_400K(i.AT_Temp[tt]) * s.AT_HO2_conc +
                           f.CH3O2_NO_200_430K(i.AT_Temp[tt]) * s.AT_CH3O2_conc
                           )
        # Calculate k_OH, the kinetic rate constant of the reaction NO2 + OH + M -> HNO3 + M
        N2_conc = 0.78 * i.AT_Pres[tt] * 0.001 * 100000 * p.N_A / (p.R * i.AT_Temp[tt] * 1000000)
        k_OH = 1.8E-30 * (i.AT_Temp[tt]/300)**-3.0 * N2_conc
        # Calculate k_O3, the kinetic rate constant of the reaction NO2 + O3 -> NO3 + O2
        k_O3 = 1.2E-13 * exp(-2450/i.AT_Temp[tt])
        s.AT_frac_OH_oxidation = (k_OH * s.AT_OH_conc) / (k_OH * s.AT_OH_conc + k_O3 * s.AT_O3_conc)
        ## Local cycling (Leigthon) of the photolytic NO2 flux
        # Store D17O value in FP before cycling (PSS)
        s.FP_D17O_before_PSS = s.FP_D17O

        # We calculate the photochemical steady-state (PSS) to obtain D17O of NO2
        if opt.NO2_cycling != 0:
            if opt.NO2_cycling == 1:
                # Cycling is ON, we calculate the photochemical steady-state (PSS) to obtain D17O of NO2
                s.AT_D17O_NO2_PSS = s.AT_alpha * (1.18 * i.AT_D17O_O3b[tt] + 6.6)
                #print 'at PSS ,D17O: '+str(s.AT_D17O_NO2_PSS)
            elif opt.NO2_cycling == 9999:
                # Cycling is turned ON but with the option D17O(NO2) = 0 (for sensitivity tests)
                s.FP_D17O = 0
                s.AT_D17O_NO2_PSS = 0
            elif opt.NO2_cycling < 0:
                # Cycling is turned ON but with the option D17O(NO2) = abs(opt.NO2_cycling) (for sensitivity tests)
                s.FP_D17O         = float(abs(opt.NO2_cycling))
                s.AT_D17O_NO2_PSS = float(abs(opt.NO2_cycling))
            elif opt.NO2_cycling == 2:
                # Cycling is turned ON but we set alpha to a given constant value (for sensitivity tests)
                s.AT_alpha = opt.NO2_set_alpha
                s.AT_D17O_NO2_PSS = s.AT_alpha * (1.18 * i.AT_D17O_O3b[tt] + 6.6)
            elif opt.NO2_cycling == 3:
                # Cycling is turned ON but we set alpha not to get below a given minimum value (for sensitivity tests)
                if s.AT_alpha < opt.NO2_set_alpha:
                    s.AT_alpha = opt.NO2_set_alpha
                s.AT_D17O_NO2_PSS = s.AT_alpha * (1.18 * i.AT_D17O_O3b[tt] + 6.6)
                # Set FP D17O value
            # Only if the FP flux is non-zero
            if s.FP_mass == 0.:
                s.FP_D17O = 0.
            else:
                s.FP_D17O = s.AT_D17O_NO2_PSS
        elif opt.NO2_cycling == 0:
            # Cycling is OFF, NO2 keeps the snowpack D17O signature
            # s.FP_D17O is unchanged except if FP flux is zero then set values to 0
            if s.FP_mass == 0.:
                s.FP_D17O = 0.
            s.AT_D17O_NO2_PSS = 0

        ## Local oxidation of the photolytic NO2 flux
        if opt.NO2_oxidation != 0:
            if opt.NO2_oxidation == 1:
                # Oxidation is turned ON, add. O can come from OH or O3 depending on relative reaction rates
                # Assume 60% of OH comes from photolysis of HONO and HONO comes from NO2- + H+ and NO2- comes from photolysis of NO3-
                # Assume stochastic repartition of 17O-excess in NO3-, NO2- and HONO
                # Assume 40% of OH comes from sources with D17O=0. O3 pathway neglected
                s.AT_D17O_OH   = 0 #0.6 * s.FP_D17O_before_PSS
                s.AT_D17O_O3   = 1.23 * i.AT_D17O_O3b[tt] + 9.02 # Berhanu et al. (2012)
                s.AT_D17O_addO = s.AT_frac_OH_oxidation * s.AT_D17O_OH + (1 - s.AT_frac_OH_oxidation) * s.AT_D17O_O3
#                print 'the add O: '+str(s.AT_D17O_addO)
            elif opt.NO2_oxidation == 2:
                # Oxidation is turned ON, add. O comes from OH only
                # Assume 60% of OH comes from photolysis of HONO and HONO comes from NO2- + H+ and NO2- comes from photolysis of NO3-
                # Assume stochastic repartition of 17O-excess in NO3-, NO2- and HONO
                # Assume 40% of OH comes from sources with D17O=0. O3 pathway neglected
                s.AT_D17O_OH   = 0.6 * s.FP_D17O_before_PSS
                s.AT_D17O_O3   = 0
                s.AT_D17O_addO = s.AT_D17O_OH
            elif opt.NO2_oxidation == 3:
                # Oxidation is turned ON, add. O comes from O3 only
                s.AT_D17O_OH   = 0
                s.AT_D17O_O3   = 1.23 * i.AT_D17O_O3b[tt] + 9.02 # Berhanu et al. (2012)
                s.AT_D17O_addO = s.AT_D17O_O3
            elif opt.NO2_oxidation == 9999:
                # Oxidation is turned ON but with the option D17O(add. O) = 0 (for sensitivity tests)
                s.AT_D17O_OH   = 0
                s.AT_D17O_O3   = 0
                s.AT_D17O_addO = 0
            elif opt.NO2_oxidation < 0:
                # Oxidation is turned ON but with the option D17O(add. O) = abs(opt.NO2_oxidation) (for sensitivity tests)
                s.AT_D17O_OH   = 0
                s.AT_D17O_O3   = 0
                s.AT_D17O_addO = abs(opt.NO2_oxidation)
            # Set value to D17O(NO2)
            if s.FP_mass == 0.:
                s.FP_D17O = 0.
            else:
                s.FP_D17O = (2 * s.FP_D17O +
                             1 * s.AT_D17O_addO
                             ) / 3
        elif opt.NO2_oxidation == 0:
            # No local oxidation, D17O of the FP flux remains the same (D17O(HNO3) = D17O(NO2))
            s.AT_D17O_addO = 0
            s.AT_D17O_OH   = 0
            s.AT_D17O_O3   = 0

        ### STEP 9 : MASS BALANCE IN THE ATMOSPHERE #############################################    
        # MASS BALANCE IN THE ATMOSPHERE
        # Calculate the deposited flux FD (will be deposited at next time step)
        # Depends on nitrate recycling option.
        # If ON, means that FD is locally oxidized to reform HNO3
        if opt.local_recycling == 1:
            s.FD_mass = (i.AT_mass[tt-1] - i.AT_mass[tt]) / p.Dt + (1 - i.AT_f_exp[tt]) * s.FP_mass + (i.FS_mass[tt] + i.FT_mass[tt])  # kgN.m-2.s-1
            s.FE_mass = i.AT_f_exp[tt] * s.FP_mass                                                   # kgN.m-2.s-1

            # Mass balance in atmosphere
            if s.FP_mass == 0:
                s.AT_d15N = (i.AT_mass[tt-1] * s.AT_d15N + p.Dt * (i.FS_d15N[tt] * i.FS_mass[tt] + i.FT_d15N[tt] * i.FT_mass[tt]) - p.Dt * s.FD_mass * p.AT_eps15_dep ) / (i.AT_mass[tt] + p.Dt * s.FD_mass)
            else:
                s.AT_d15N = (i.AT_mass[tt-1] * s.AT_d15N + p.Dt * ((1 - i.AT_f_exp[tt]) * 
                        s.FP_d15N * s.FP_mass + i.FS_d15N[tt] * i.FS_mass[tt] + i.FT_d15N[tt] * i.FT_mass[tt]) - p.Dt * s.FD_mass *  p.AT_eps15_dep
                             ) / (i.AT_mass[tt] + p.Dt * s.FD_mass)

            # Isotopic composition of FE
            s.FE_d15N = s.FP_d15N
            s.FE_D17O = s.FP_D17O
            # The mass balance gives the expression of D17O in atmospheric nitrate
            # If the FP flux is zero, we remove the photolytic flux in the budget
            if s.FP_mass == 0.:
                s.AT_D17O = (s.AT_mass * s.AT_D17O + p.Dt * (i.FS_D17O[tt] * i.FS_mass[tt] + i.FT_D17O[tt] * i.FT_mass[tt])) / (s.AT_mass + p.Dt * s.FD_mass)
            else:
                s.AT_D17O = (s.AT_mass * s.AT_D17O + p.Dt * ((1 - i.AT_f_exp[tt]) * s.FP_D17O * s.FP_mass + i.FS_D17O[tt] * i.FS_mass[tt] + i.FT_D17O[tt] * i.FT_mass[tt]))/ (s.AT_mass + p.Dt * s.FD_mass)
            s.FD_d15N = s.AT_d15N + p.AT_eps15_dep    # same d15N value as in the atmosphere but fractionation at deposition
            s.FD_D17O = s.AT_D17O                     # same D17O value as in the atmosphere but NO fractionation at deposition
        # If OFF, means that FP is removed from the atmospheric box
        elif opt.local_recycling == 0:
            s.FD_mass = (i.AT_mass[tt-1] - i.AT_mass[tt]) / p.Dt +  (i.FS_mass[tt] + i.FT_mass[tt])  # kgN.m-2.s-1
            s.FE_mass = 0                                                  # kgN.m-2.s-1
            # Mass balance in atmosphere

            s.AT_d15N = (i.AT_mass[tt-1] * s.AT_d15N + p.Dt * (
                    i.FS_d15N[tt] * i.FS_mass[tt] + i.FT_d15N[tt] * i.FT_mass[tt]) - p.Dt * s.FD_mass * p.AT_eps15_dep
                             ) / (i.AT_mass[tt] + p.Dt * s.FD_mass)
            s.AT_D17O = (s.AT_mass * s.AT_D17O + p.Dt * (1 - i.AT_f_exp[tt]) * (i.FS_D17O[tt] * i.FS_mass[tt] + i.FT_D17O[tt] * i.FT_mass[tt])) / (s.AT_mass + p.Dt * s.FD_mass)   
            # Isotopic composition of FD 
            s.FD_d15N = s.AT_d15N + p.AT_eps15_dep    # same d15N value as in the atmosphere but fractionation at deposition
            s.FD_D17O = s.AT_D17O                     # same D17O value as in the atmosphere but NO fractionation at deposition
            # Isotopic composition of FE
            s.FE_d15N = s.FP_d15N
            s.FE_D17O = s.FP_D17O
            
        ### STEP 10 : OUTPUTS EXPORT AND CALCULATION #############################################
        # Exportation of each column array to the (n x t) arrays
        #the new output 
        #nitrate flux
        o.FD_mass[0, t]         = s.FD_mass
        o.FD_d15N[0, t]         = s.FD_d15N
        o.FD_D17O[0, t]         = s.FD_D17O
        o.FP_mass[0, t]           = s.FP_mass
        o.FP_d15N[0, t]           = s.FP_d15N
        o.FP_D17O[0, t]           = s.FP_D17O
        o.AT_mass[0, t]           = s.AT_mass
        o.AT_d15N[0, t]           = s.AT_d15N
        o.AT_D17O[0, t]           = s.AT_D17O
        o.AT_D17O_addO[0, t]      = s.AT_D17O_addO
        o.AT_D17O_NO2_PSS[0, t]   = s.AT_D17O_NO2_PSS
        o.AT_alpha[0, t]          = s.AT_alpha
        o.AT_JNO2[0, t]           = s.AT_JNO2
        o.AT_OH_conc[0, t]        = s.AT_OH_conc
        o.AT_HO2_conc[0, t]       = s.AT_HO2_conc
        o.AT_CH3O2_conc[0, t]     = s.AT_CH3O2_conc
        o.AT_D17O_OH[0, t]        = s.AT_D17O_OH
        o.eps_surface[0, t]       = s.SN_eps15[0]

        ### Returns computed duration of the computed year
        if tt == p.nb_ts-1:
            tps_fin = time.time()
            print( 'Year ' + str(int(t/p.nb_ts+1)) + ' computed in ' + str(str(int((tps_fin - tps_deb)*10)/10)) + ' s')
##
##            print s.SN_mass[0:len(s.SN_mass)-1, 0].sum(), (s.SN_mass[0:len(s.SN_mass)-1, 0]*s.SN_d15N[0:len(s.SN_mass)-1, 0]).sum(), (s.SN_mass[0:len(s.SN_mass)-1, 0]*s.SN_D17O[0:len(s.SN_mass)-1, 0]).sum(),(s.SN_mass[0:len(s.SN_mass)-1, 0]*s.SN_d15N[0:len(s.SN_mass)-1, 0]).sum()/(s.SN_mass[0:len(s.SN_mass)-1, 0].sum()), (s.SN_mass[0:len(s.SN_mass)-1, 0]*s.SN_D17O[0:len(s.SN_mass)-1, 0]).sum()/(s.SN_mass[0:len(s.SN_mass)-1, 0].sum())
        if t == p.cp_dur-1: #in the last time step:
            o.SN_depth = s.SN_depth
            o.SN_thick = s.SN_thick
            o.SN_mass = s.SN_mass
            o.SN_d15N = s.SN_d15N
            o.SN_D17O = s.SN_D17O
            o.SN_date = s.SN_date
        o.frem[:t+1,t] = (f.find_eps(s.SN_date,frem,t))*(1-p.f_cage)+p.f_cage
### This is the end of the computational loop
#################################################################################################