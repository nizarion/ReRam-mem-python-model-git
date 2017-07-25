"""
Python Code for Fitting OFF state data to find Poole-Frenkel a and b parameters:

This program finds the Poole-Frenkel coefficients to describe nonlinearity in the simulator
Written by Andrew J. Lohn and Patrick R. Mickel.
Created on Tue May 20 15:26:10 2014.
"""
import numpy as np
from scipy.optimize import curve_fit
import pylab as pl
import csv

f = open('Nonlinearity-sn7cuc13-2ndneg-257-308v3.csv','rb')           #This is the experimental current-voltage data from only the OFF state. 
expData = csv.reader(f)                     #The resistance should be changing due to nonlinearity only, not state change.
expVs = []
expIs = []
for row in expData:
    expVs.append(float(row[0]))
    expIs.append(float(row[1]))
expVs = np.array(expVs)
expIs = np.array(expIs)
f.close()

#These should all be the same as will be used in the simulation code
########  Set the Device Parameters  ###########
sigSat = 2.47 * 10**5                       #Saturation conductivity was previously found from fitting the ON switching data @fitting-ON
numberOfShells = 200                            #Number of Grids in the Radial Direction
thickness = 6 * 10**-9                         #Vertical Thickness in the z direction
simRadiusRange = 20 * 10**-9                    #The full spacial range of the simulation
shellSize = simRadiusRange/numberOfShells       #Calculate GridSpacing

rmax = 4.39 * 10**-10                     #@ON #This is the radius that was previously found from fitting the ON switching data @fitting-ON
######## End of Parameters Section ###########
numAtrmax = int(rmax/shellSize)

shellRadii = [(i+1)*shellSize for i in xrange(numAtrmax)]
shellGeometries = []                                                #Build an Array of Geometries for Each Shell
for radius in shellRadii:
    shellGeometries.append(3.14*(radius**2 - (radius-shellSize)**2) / thickness)

shellConcentrations = 50*np.ones(numAtrmax)             #This sets the filament to the minimum concentration for all shells up to the largest radius during ON switching
minConc = 50

a = 1                       #These are initialized values with no real significance was0.001
b = 1 * 10**-5                           #was 7

def findShellConds(shellConcentrations,voltageVal):         #This function finds the shell conductivities given their concnentrations and the applied voltage
    shellConductivities = []
    for shellConc in shellConcentrations:
        if shellConc == 0:
            shellConductivities.append(0)
        else:
            shellConductivities.append(sigSat * abs(minConc-shellConc)/minConc + sigSat * (1 - abs(minConc-shellConc)/minConc) * abs(voltageVal) * a*np.exp(b*(abs(voltageVal)**0.5)))
    shellConductivities = np.array(shellConductivities)
    return shellConductivities

def nonLinFunc(Voltages,a,b):                   #This function calcuates currents the way the simulator would for a range of voltages and known a and b values
    Currents = []
    for voltageVal in Voltages:
        shellConductivities = []
        for shellConc in shellConcentrations:
            if shellConc == 0:
                shellConductivities.append(0)
            else:
                shellConductivities.append(sigSat * abs(minConc-shellConc)/minConc + sigSat * (1 - abs(minConc-shellConc)/minConc) * abs(voltageVal) * a*np.exp(b*(abs(voltageVal)**0.5)))
        shellConductivities = np.array(shellConductivities)
        resistance = 1 / np.sum(shellConductivities * shellGeometries)
        Currents.append(voltageVal/resistance)
    if a > 1*10**-21 and a < 1*10**3 and b > 1*10**-1 and b < 1*10**3:      #Return the currents if the a and b values are in a reasonable range.
        return Currents
    else:                           #If the a and be values are not reasonable then return a large value to create a large error during the fitting routine.
        return 1*10**6
    

##########  Fit the Data  ###########
pinit = [1*10**-9,1]                    #Initialize the a and b values for fitting
popt,pconv = curve_fit(nonLinFunc,expVs,expIs,p0=pinit)     #Fit to find the a and b values
print "the a and b values are:"
print popt
fitIs = nonLinFunc(expVs,popt[0],popt[1])       #Create a list of fitted currents for comparison to experimental data
pl.plot(expVs,expIs,'bo')           #Plot the experimental and fitted curves.
pl.plot(expVs,fitIs,'b-')
pl.show()
