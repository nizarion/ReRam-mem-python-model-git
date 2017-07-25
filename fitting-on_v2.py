"""
Python Code for Fitting ON switching data to find device parameters:

This program finds the saturation conductivity, activation temperature, electrode thermal conductivity and filament radius from the ON switching data of an I-V hysteresis loop
Written by Andrew J. Lohn and Patrick R. Mickel.
Created on Tue May 20 15:26:10 2014.
Related theory is described in "Isothermal Switching and Detailed Filament Evolution in Memristive Systems" published in Advanced Materials.
"""
import numpy as np
from scipy.optimize import curve_fit
import pylab as pl
import csv

########  Set the Device Parameters  ###########
de = 30 * 10**-9            #Electrode thickness in meters #done
do = 6 * 10**-9             #Oxide thickness in meters #done
Lwf = 2.44 * 10**-8         #Wiedemann-Franz constant
Trt = 296                   #Ambient temperature in Kelvin #done
######## End of Parameters Section ###########                 #Ambient temperature in Kelvin #done

def onSwitch(R,ke,sigSat,Tc):               #This function calculates Powers for a list of Resistances given the device properties
    Ar = 2*ke*do/(sigSat*de)                #This is a prefactor in the power-resistance relation
    Rmin = ke / (4*3.14*(sigSat**2)*Lwf*Tc*de)      #This is the minimum resistance
    
    P = Ar * (Tc - Trt) / (R - Rmin)        #Calculate power as a function of resistance
    if ke > 10 and ke < 500 and sigSat > 1*10**3 and sigSat < 1*10**6 and Tc > 500 and Tc < 3000:       #if the device parameters are reasonable then return the power
        return P
    else:                                   #If the parameters are unreasonable then return a large number so the error will be large in the fitting routine.
        return 1*10**6
    
def PRconvert(IV):                          #This function converts current-voltage data to power-resistance data
    PR = []
    for pair in IV:
        power = pair[0]*pair[1]             #Power is current times voltage
        resistance = pair[1]/pair[0]        #Resistance is voltage divided by current
        PR.append([power,resistance])
    return PR


###########  Load and Fit Data  ##############
IVdata = csv.reader(open('onSwitch-sn7cuc13-2ndpos-309-390v2.csv','rb'))      #A csv file is used where the columns are voltage then current. 
IV = []                                             #This data is only the ON switching while the resistance is changing.
for row in IVdata:
    IV.append([float(row[1]),float(row[0])])        #This reverses the order to match with the order in the variable name - IV is I then V
IV = np.array(IV)

PR = PRconvert(IV)                                  #Convert the current-voltage data to power-resistance data
PR = np.array(PR)

pinit = [105.5,5*10**5,1750]                        #Initialize the values for [electrode thermal conductivity (ke), saturation conductivity (sigSat), activation temperature (Tc)]
popt,pconv = curve_fit(onSwitch,PR[:,1],PR[:,0],p0=pinit)       #Fit the data to the power-resistance relation to determine the ke, sigSat, and Tc values
print "ke, sigSat, and Tc values:"
print popt

PRfit = onSwitch(PR[:,1],popt[0],popt[1],popt[2])   #Create a list of power values from the fitting parameters for comparison to data

resistance = IV[-1][1]/IV[-1][0]
print "radius at the end of the ON switching:"
print (do/(popt[1]*3.14*resistance))**0.5           #print out the radius at the end of the ON switching

pl.plot(PR[:,1],PR[:,0],'bo')                       #plot the data and the fit
pl.plot(PR[:,1],PRfit,'b')
pl.show()
