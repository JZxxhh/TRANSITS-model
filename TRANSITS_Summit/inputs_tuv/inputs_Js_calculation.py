 ### Import modules
from __future__ import division
import numpy as np
from StringIO import StringIO
from pylab import *

import os
import glob
#########################################################################################

#########################################################################################
### PARAMETERS
#########################################################################################
N_A    = 6.02214179E+23 # Avogadro Number
factor = 1e-6

new_nm = np.arange(280, 351, 1)

# Calculate XS14 using the fit in Ayalneh et al 2014 on Chu et al 2003 278K data

A = 192.5 # M-1
C = 34052 # cm-1
W = 3573  # cm-1
S = 0.9

wl         = np.arange(280, 351, 1)
E          = 1/wl * 1e7
X          = (E - C) / W
XSoverE    = (1000*np.log(10)/N_A) * factor * A * (1 - S*X) * np.exp(-X*X * (1 - S*X + 0.5*(S*X)*(S*X)))
XS14_abs   = XSoverE * E


# Calculate XS15 using Ayalneh et al 2014 and a width reduction factor of 1.0% and a center shift of -32.5 cm-1 to model the 14N to 15N isotopic substitution effetcs

wrf = 0.01  # width reduction factor
DC  = -32.5 # cm-1

A = A/(1-wrf)  # M-1
C = C - DC # cm-1
W = W*(1-wrf)  # cm-1
S = 0.895

E          = 1/wl * 1e7
X          = (E - C) / W
XSoverE    = (1000*np.log(10)/N_A) * factor * A * (1 - S*X) * np.exp(-X*X * (1 - S*X + 0.5*(S*X)*(S*X)))
XS15_abs   = XSoverE * E


### CROSS SECTION O3 ::: Reading the input file : inputs_Xsection_O3_Malicet1995_243K.txt
inp = open("O3_cross_section/inputs_Xsection_O3_Malicet1995_243K.txt","r")
XSO3 = np.genfromtxt(inp, skip_header=1, delimiter="")
XSO3_nm  = new_nm
XSO3_abs = np.interp(new_nm, XSO3[:, 0], XSO3[:, 1])

##### CROSS SECTION NO2 ::: Reading the input file : inputs_Xsection_NO2_Bogumil2003_243K_interpol.txt
inp = open("NO2_cross_section/inputs_Xsection_NO2_Bogumil2003_260K.txt","r")
XSNO2 = np.genfromtxt(inp, skip_header=1, delimiter="")
XSNO2_nm  = new_nm
XSNO2_abs = np.interp(new_nm, XSNO2[:, 0], XSNO2[:, 1])
#########################################################################################

path= "TUV_actinic_flux"

### Create output of the J_O3 calculation
J_O3 = np.zeros((41, len(glob.glob( os.path.join(path, '*.txt') )) + 1))
JO3_header = open("JO3_header.txt", "w" )
JO3_header.write("SZA\Ozone_column_DU")         # write first column in header
JO3_header.write( "\t" )

##### Create output of the J_NO2 calculation
J_NO2 = np.zeros((41, len(glob.glob( os.path.join(path, '*.txt') )) + 1))
JNO2_header = open("JNO2_header.txt", "w" )
JNO2_header.write("SZA\Ozone_column_DU")         # write first column in header
JNO2_header.write( "\t" )

#########################################################################################
### Input reading and preparation
#########################################################################################
count = 0    # counter for J_O3 extraction
                   
for infile in glob.glob( os.path.join(path, '*.txt') ):
    ### Gives the file a name
    infile_name = filter(lambda x: x.isdigit(), infile) + 'DU'
    
    ### ACTINIC FLUX ::: Reading the input file : infile
    inp = open(infile, "r")
    read_for_header = np.genfromtxt(open('header_ini.txt','r'))
    header = read_for_header[:]      # This is only to read the header (= depths)
    inp = open(infile,"r")
    TUV_DU = np.genfromtxt(inp, skip_header=0, delimiter="")

    # Write data in J_O3 header file
    JO3_header.write(infile_name)    # write the ozone column value
    JO3_header.write( "\t" )

##    # Write data in J_NO2 header file
    JNO2_header.write(infile_name)    # write the ozone column value
    JNO2_header.write( "\t" )

    count = count + 1

    ### PREPARATION
    (rows, cols) = TUV_DU.shape  #(71*41)*119

    J_14NO3  = np.zeros((41, cols+1))
    J_15NO3  = np.zeros((41, cols+1))
    I_abs    = np.zeros((71, 1))      # transfer array for actinic flux values

    for i in range(41):
        # Write the SZA value
        J_14NO3[i, 0] = J_15NO3[i, 0] = J_O3[i, 0] = 90-i
        J_NO2[i, 0] = 90-i
        # Extract the actinic flux at surface
        for j in range(71):
            I_abs[j, 0] = TUV_DU[i*71 + j, 0]
        # Calculate J_O3 and J_NO2 at surface at a given SZA
        J_O3[i, count]      = (I_abs[:, 0] * XSO3_abs).sum()
        J_NO2[i, count]     = (I_abs[:, 0] * XSNO2_abs).sum()
        # Extract the actinic flux at different depths
        for d in range(cols):
            for j in range(71):         # 71 is the number of wavelengths
                I_abs[j, 0] = TUV_DU[i*71 + j, d]
            # Calculate J14 and J15 at different depth and at a given SZA
            J_14NO3[i, d+1] = (I_abs[:, 0] * XS14_abs).sum()  # the factor 1000 is because of the conversion of liters to cm3
            J_15NO3[i, d+1] = (I_abs[:, 0] * XS15_abs).sum()
    
    # flipud the array
    J_14NO3[:,1:] = np.flipud(J_14NO3[:,1:])
    J_15NO3[:,1:] = np.flipud(J_15NO3[:,1:])
    
    ### OUTPUT EXPORT in files
    # First create the header file
    out = open('header.txt', "w" )
    out.write("SZA\depth(m)")
    out.write("\t")
    for i in range(len(header)-2):
        out.write(str(header[i+2] * 0.01))   # depth conversion to meters
        out.write( "\t" )
    out.write( "\n" )
    out.close()

    # Then save data in another file
    np.savetxt('data14.txt', J_14NO3, delimiter="\t")
    np.savetxt('data15.txt', J_15NO3, delimiter="\t")

    # Eventually merge the two files and delete the intermediate ones
    open('inputs_J14NO3_' + infile_name +'.txt','w').write(open('header.txt','r').read()+open('data14.txt','r').read())
    open('inputs_J15NO3_' + infile_name +'.txt','w').write(open('header.txt','r').read()+open('data15.txt','r').read())
    os.remove('header.txt')
    os.remove('data14.txt')
    os.remove('data15.txt')

# Finish writing in J_O3 header file
JO3_header.write( "\n" )
JO3_header.close()

J_O3[:,1:] = np.flipud(J_O3[:,1:])
J_NO2[:,1:] = np.flipud(J_NO2[:,1:])

### Finish writing in J_NO2 header file
JNO2_header.write( "\n" )
JNO2_header.close()

##for i in range(41):
##    I_310[i, 0] = 90 - i

### OUTPUT EXPORT in fil" for the 310 nm extraction
# save data in temporary file
np.savetxt('JO3_TEMP.txt', J_O3, delimiter="\t")
np.savetxt('JNO2_TEMP.txt', J_NO2, delimiter="\t")

# Eventually merge the two files and delete the intermediate ones
open('inputs_JO3.txt','w').write(open('JO3_header.txt','r').read()+open('JO3_TEMP.txt','r').read())
os.remove('JO3_header.txt')
os.remove('JO3_TEMP.txt')

open('inputs_JNO2.txt','w').write(open('JNO2_header.txt','r').read()+open('JNO2_TEMP.txt','r').read())
os.remove('JNO2_header.txt')
os.remove('JNO2_TEMP.txt')


