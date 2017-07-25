"""
This program simulates the electrical response of memristors or RRAM based on the physical parameters of the device.

Written by Andrew J. Lohn and Patrick R. Mickel
Created on Tue May 20 15:26:10 2014

Related theory is described in "Isothermal Switching and Detailed Filament Evolution in Memristive Systems" published in Advanced Materials.
"""
import numpy as np
import pylab as pl
import math

########  Set the Device Parameters  ###########
delT = 1513                    #Activation temperature in degrees C from the ambient temperature
Trt = 296                       #done #Ambient temperature in degrees C
de = 30 * 10**-9                #done #Electrode thickness
do = 6 * 10**-9                	#done #Oxide thickness
ke = 195                     	#done1 #@fittin-ON #Effective electrode thermal conductivity - Usually higher than bulk values due to spreading effects 
Lwf = 2.44 * 10**-8             #Wiedemann-Franz constant #stathera

sigSat = 2.47 * 10**5           #@fittin-ON #Saturation Conductivity in the Filament 
maxConc = 500                   #The maximum concentration of oxygen vacancies or dopants in arbitrary units 
minConc = 100                    #The minimum concentration of oxygen vacancies or dopants in arbitrary units 
a = 1.06946442e-13                     #@OFF #Amplitude coefficient for Poole-Frenkel contribution
b = 7.08532181e-01                     #@OFF #Nonlinearity coefficient for Poole-Frenkel contribution 
######## End of Parameters Section ###########

Ar = 2*ke*do/(sigSat * de)                                  #Coefficient for radius change equation is composed of previously defined parameters
Rmin = ke / (4*3.14*(sigSat**2)*Lwf*(Trt+delT)*de)          #Minimum resistance value is composed of previously defined parameters

numberOfShells = 200                            #Set the number of shells for the simulation
simRadiusRange = 20 * 10**-9                    #The full spacial range of the simulation
shellSize = simRadiusRange/numberOfShells       #Calculate grid spacing radially

shellRadii = [(i+1)*shellSize for i in xrange(numberOfShells)]      #Build an array of outer edge locations for each shell
shellGeometries = []                                                #Build an array of aeometries for each shell to be used with conductivity to calculate resistance
for radius in shellRadii:
    shellGeometries.append(3.14*(radius**2 - (radius-shellSize)**2) / do)

shellConcentrations = np.ones(numberOfShells)       #Initialize the shell concentrations
#Forming
formingRadius = (2 * 10**-9)/shellSize              #Set the forming radius # was formingRadius = (2 * 10**-9)/shellSize   
for i in xrange(len(shellConcentrations)):
    if i < formingRadius:
        shellConcentrations[i] = 80                 #Set the post-forming concentration - This can be a low value if the device is turned off after forming
    else:
        shellConcentrations[i] = 0                  #Set the concentration to 0 outside the formed region. The minimum concentration is higher than this region which is as-deposited.


##########  Establish Functions for Changing Resistance  #############
def findShellConds(shellConcentrations,voltageVal):             #Calculate the shell conductivities from the shell concentrations and the applied voltage
    shellConductivities = []
    for shellConc in shellConcentrations:
        if shellConc == 0:
            shellConductivities.append(0)                       #Enforce that the as-deposited region does not contribute.
        else:                                                   #Otherwise use a combination of Ohmic and Poole-Frenkel weighted by concentration
            shellConductivities.append(sigSat * abs(minConc-shellConc)/minConc + sigSat * (1 - abs(minConc-shellConc)/minConc) * abs(voltageVal) * a*np.exp(b*(abs(voltageVal)**0.5)))
    shellConductivities = np.array(shellConductivities)
    return shellConductivities

#Define a function for sourcing current
cSourceVoltages = np.linspace(0,2,100)                          #Voltage is needed to calculate conductivity. These values are stepped through to find the voltage for that calculation.
def CurrentSourcing(currentVal):
    global shellConcentrations
    if currentVal < 0:                                          #Negative current implies conductivity change  
        for voltageVal in cSourceVoltages:                                                      #Find the voltage that corresponds to the sourced current by stepping through values
            shellConductivities = findShellConds(shellConcentrations,voltageVal)
            resistance = 1/np.sum(shellConductivities*shellGeometries)
            currentTest = -voltageVal/resistance
            if currentTest <= currentVal:                       #If the calculated current is more negative than the sourced current then stop stepping through voltages
                break
        for i in xrange(len(shellConcentrations)):                                              #Find the location of the edge of the inner-most filament
            if shellConcentrations[i] < np.max(shellConcentrations):
                break
        shellConductivities = findShellConds(shellConcentrations,voltageVal)
        resistance = 1 / np.sum(shellConductivities[0:i+1] * shellGeometries[0:i+1])                #Calculate the resistance of the inner-most filament        
        rop = (i+1)*shellSize                                                                   #Approximate operating radius as current filament radius
        Asig = 8*(do**2)*Lwf*(Trt+delT)/(rop**2)                                                #Calculate new Asig with new operating radius
        Rmax = 4*(do**2)*Lwf*(Trt+delT)*de/(3.14*ke*(rop**4))                                   #Calculate new Rmax with new operating radius
        while (currentVal**2)*resistance >= Asig*delT/(Rmax - resistance) and resistance < Rmax:                       #While sufficient power is supplied for activation
            shellMaxConc = np.max(shellConcentrations)
            if shellMaxConc <= minConc:
                break
            for i in xrange(len(shellConcentrations)):                                          #Decrease the conductivity by one step for every shell in the inner-most filament
                if shellConcentrations[i] == shellMaxConc:
                    shellConcentrations[i] += -1
                else: 
                    break
            shellConductivities = findShellConds(shellConcentrations,voltageVal)
            resistance = 1 / np.sum(shellConductivities[0:i+1] * shellGeometries[0:i+1])
            rop = (i+1)*shellSize                                                               #Calculate new operating radius and associated parameters
            Asig = 8*(do**2)*Lwf*(Trt+delT)/(rop**2)
            Rmax = 4*(do**2)*Lwf*(Trt+delT)*de/(3.14*ke*(rop**4))

    if currentVal >= 0:                                         #Positive current implies radius change
        for voltageVal in cSourceVoltages:                                                      #Find the voltage that corresponds to the sourced current by stepping through values
            shellConductivities = findShellConds(shellConcentrations,voltageVal)
            resistance = 1/np.sum(shellConductivities*shellGeometries)
            currentTest = voltageVal/resistance
            if currentTest >= currentVal:                       #If the calculated current is more negative than the sourced current then stop stepping through voltages
                break
        for i in xrange(len(shellConcentrations)):                                              #Find the location of the edge of the inner-most filament
            if shellConcentrations[i] < np.max(shellConcentrations):
                break
        shellConductivities = findShellConds(shellConcentrations,voltageVal)
        resistance = 1 / np.sum(shellConductivities[0:i+1] * shellGeometries[0:i+1])                #Calculate the resistance of the inner-most filament
        while (currentVal**2)*resistance > Ar*delT / (resistance-Rmin) and resistance > Rmin:
            for i in xrange(len(shellConductivities)):
                if shellConcentrations[i] < maxConc:
                    shellConcentrations[i] = maxConc                                            #Find the first shell that is not saturated and move it to saturation
                    break                                                                       #Break so only one shell saturates at a time
            shellConductivities = findShellConds(shellConcentrations,voltageVal)
            resistance = 1 / np.sum(shellConductivities[0:i+1] * shellGeometries[0:i+1])            #Calculate the new resistance of the inner-most filament
            if i == len(shellConductivities)-1:                                                 #If all the shells are saturated then exit
                break
    shellConductivities = findShellConds(shellConcentrations,voltageVal)
    resistance = 1 / np.sum(shellConductivities * shellGeometries)                              #Calculate the resistance of the entire filament, not just the inner-most filament
    return currentVal * resistance                                                              #Return the voltage

#Define a function for sourcing voltage
def VoltageSourcing(voltageVal):
    global shellConcentrations
    if voltageVal <= 0:                                                                         #Negative voltage implies conductivity change
        for i in xrange(len(shellConcentrations)):                                              #Find the location of the edge of the inner filament
            if shellConcentrations[i] < np.max(shellConcentrations):
                break
        shellConductivities = findShellConds(shellConcentrations,voltageVal)
        resistance = 1 / np.sum(shellConductivities[0:i+1] * shellGeometries[0:i+1])                #Calculate the resistance of the inner-most filament        
        rop = (i+1)*shellSize                                                                   #Approximate operating radius as current filament radius
        Asig = 8*(do**2)*Lwf*(Trt+delT)/(rop**2)                                                #Calculate new Asig with new operating radius
        Rmax = 4*(do**2)*Lwf*(Trt+delT)*de/(3.14*ke*(rop**4))                                   #Calculate new Rmax with new operating radius
        while (voltageVal**2)/resistance >= Asig*delT/(Rmax - resistance) and resistance < Rmax:                       #While sufficient power is supplied for activation
            shellMaxConc = np.max(shellConcentrations)
            if shellMaxConc <= minConc:
                break
            for i in xrange(len(shellConcentrations)):                                          #Decrease the conductivity by one step for every shell in the inner-most filament
                if shellConcentrations[i] == shellMaxConc:
                    shellConcentrations[i] += -1
                else:
                    break
            shellConductivities = findShellConds(shellConcentrations,voltageVal)
            resistance = 1 / np.sum(shellConductivities[0:i+1] * shellGeometries[0:i+1])
            rop = (i+1)*shellSize                                                                #Calculate new operating radius and associated parameters
            Asig = 8*(do**2)*Lwf*(Trt+delT)/(rop**2)
            Rmax = 4*(do**2)*Lwf*(Trt+delT)*de/(3.14*ke*(rop**4))
    
    if voltageVal > 0:                                                                          #Positive voltage implies radius change
        for i in xrange(len(shellConcentrations)):                                              #Find the location of the edge of the inner-most filament
            if shellConcentrations[i] < np.max(shellConcentrations):
                break
        shellConductivities = findShellConds(shellConcentrations,voltageVal)
        resistance = 1 / np.sum(shellConductivities[0:i+1] * shellGeometries[0:i+1])                #Calculate the resistance of the inner-most filament
        while (voltageVal**2)/resistance > Ar*delT / (resistance-Rmin):# and resistance > Rmin:   #While sufficient power is supplied for activation
            for i in xrange(len(shellConductivities)):
                if shellConcentrations[i] < maxConc:
                    shellConcentrations[i] = maxConc                                            #Find the first shell that is not saturated and move it to saturation
                    break                                                                       #Break so only one shell saturates at a time
            shellConductivities = findShellConds(shellConcentrations,voltageVal)
            resistance = 1 / np.sum(shellConductivities[0:i+1] * shellGeometries[0:i+1])            #Calculate the new resistance of the inner-most filament
            if i == len(shellConductivities)-1:                                                 #If all the shells are saturated then exit
                break
    shellConductivities = findShellConds(shellConcentrations,voltageVal)
    resistance = 1 / np.sum(shellConductivities * shellGeometries)                              #Calculate the resistance of the entire filament, not just the inner-most filament
    return voltageVal / resistance                                                              #Return the current

#########  Create Input Signals and Generate Output Signals  #########################
Voltages = []           #Create empty lists to be populated during the voltage or current sweep
Currents = []
logCurrents = []
temp = -4
maxVoltage = 0

#Switch On
CurrSpacing = np.logspace(-12,-2,100)                #Use current sourcing to sweep from 0 to 6 mA # was 50 CurrSpacing = np.logspace(0,6*10**-3,5)
CurrentSweep = CurrSpacing
for val in reversed(CurrSpacing):
    CurrentSweep = np.hstack([CurrentSweep,val])        #Also sweep from 6 mA back to 0
for currentVal in CurrentSweep:
    voltage = CurrentSourcing(currentVal)               #Calculate the voltage by calling the current sourcing function
    Voltages.append(voltage)
    maxVoltage = max(Voltages)
    Currents.append(currentVal)
    if currentVal==0 :
        logCurrents.append(temp)
    else:
        logCurrents.append(np.log10(abs(currentVal)))
        temp = np.log10(abs(currentVal))
    
#Switch Off 
VoltSpacing = np.linspace(0,-maxVoltage,50)            #Use voltage sourcing to sweep from 0 to -1.0 V #was VoltSpacing = np.linspace(0,-1.0,50)  
VoltageSweep = VoltSpacing
for val in reversed(VoltSpacing):
    VoltageSweep = np.hstack([VoltageSweep,val])        #Also weep from -1.1 V to 0 
for voltageVal in VoltageSweep:
    current = VoltageSourcing(voltageVal)               #Call the voltage sourcing function to find the current at every voltage step
    Voltages.append(voltageVal)
    Currents.append(current) 
    if current==0 :
        logCurrents.append(temp)
    else:
        logCurrents.append(np.log10(abs(current)))
        temp = np.log10(abs(current))

pl.plot(Voltages,logCurrents) #pl.plot(Voltages,Currents OR logCurrents)
pl.show()
