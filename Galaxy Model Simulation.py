# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 20:09:29 2021

Authors: Edward Thomas
"""


#Imports
import numpy as np
import matplotlib.pyplot as plt


'''
These are all the equations used for the simulation
'''
#defining the formula for critical mass
def mcrit(age):
		return (10/age)**(1/2.5)


#defining E_0 by rearanging the total mass IMF formula
def E_0(mgas=10**11,maxm=25,minm=0.8):
	return mgas/((maxm**(-0.35)-minm**(-0.35))/-0.35)


#deining the equation for total number of stars in the galaxy by rearranging the IMF formula
def noofstars(mgas=10**11,maxm=25,minm=0.8):
	return E_0()*((maxm**(-1.35)-minm**(-1.35))/-1.35)
	
#defining how to find luminosity of a star
def luminosity(m):
	if m<1:
		a=5
	else:
		a=3.9
	return m**a

#defining the equation for luminosity as a factor of magnitude
def lum(mag):
	return lo*(10**(mag*-0.4))

#defining the equation for magnitude as a factor of luminosity (the previous equation rearranged)
def magtotal(ltotal):
	return -2.5*np.log10(ltotal/lo)



	
'''
These are all the values and arrays needed for the first billion years to find
luminosity- some are reused to find the change in V mag and B-V colour
'''
#the timestep for the first 1 billion years
initialtimestep=(10**3)
#creating an array of solar masses from 0.5-10
mass=np.linspace(0.8,25,121)
#creating an empty array to put the calculated number of star at each mass
noandmstar=np.zeros(len(mass))
#creating an array for time at each for the first 1 billion years
time=np.linspace(0,1,initialtimestep)
#finding the critical masses for each timestep for the first 1 billion years
criticalmass=mcrit(time)
#creating an empty array to put the values for luminosity for each mass in
l=np.zeros(len(noandmstar))
#creating an array to put the total luminosity for each timestep in the first 1 billion years
iltotal=np.zeros(initialtimestep)

'''
These are all the values and arrays needed for the last 11 billion years to find
luminosity- some are reused to find the change in V mag and B-V colour
'''
#a second timestep for the next 11 billion years as less activity is occuring
timestep2=10**3
#an array of the time for the next 11 Gyrs
time2=np.linspace(1,12,timestep2*11)
#finding the critical masses for this time
criticalmass2=mcrit(time2)
#creating a new array for the new luminosity values
l2=np.zeros(len(noandmstar))
ltotal=np.zeros(timestep2*11)

'''
These are all the new values and arrays needed to find how the V mag changes in
the first 1 billion years- some are reused to find the change in B-V colour
'''
#defining the base line luminosity
lo=3.0128*10**28
#creating an array for the values of mass and V magnitude given in the table
vmag=([6.5,4.93,4.2,1.0,-0.5,-1.4,-3.7,-4.3,-5.2])
newmass=([0.8,1,1.25,2,3,5,9,15,25])
#this empty array is for the luminosity values for each mass
lv=np.zeros(len(newmass))
#This empty array is for the number of stars in the cluster for each mass
nonm=np.zeros(len(newmass))
#This empty array is for the total luminosity at each timesteps
vmtotal=np.zeros(initialtimestep)

'''
These are all the new values and arrays needed to find how the V mag changes in
the remaining 11 billion years- some are reused to find the change in B-V colour
'''

lv2=np.zeros(len(newmass))

vmtotal2=np.zeros(timestep2*11)

'''
These are all the new values and arrays needed to find how the B-V colour changes in
the first 1 billion years
'''
#creating an array containing the B-V colour at each mass
bvmag=([1.02,0.69,0.59,0.1,-0.12,-0.18,-0.30,-0.32,-0.32])
#this empty array is for the luminosity values for each mass
lb=np.zeros(len(newmass))
#This empty array is for the total luminosity at each timesteps
bmtotal=np.zeros(initialtimestep)
#Adding the B-V and V mag values together to get an array of B magnitudes at each mass
bmag=np.add(bvmag,vmag)

'''
These are all the new values and arrays needed to find how the B-V colour changes in
the remaining 11 billion years
'''
#this empty array is for the luminosity values for each mass
lb2=np.zeros(len(newmass))
#This empty array is for the total luminosity at each timestep
bmtotal2=np.zeros(timestep2*11)
	
'''
								***SIMULATING THE FIRST 1 BILLION YEARS***
'''

for j in range(1,initialtimestep):
	'''
	This loop calculates the total luminosity
	
	First it finds the total number of stars formed for each mass of a star and adds them to itself
	after each timestep. This gives the total number of stars after each timestep.
	There is then a check to see if each mass is over the critical mass,if it is 
	then one timestep worth of stars formed at that mass is subtracted from the 
	total number. This method takes into account the time taken for the stars to form initially
	as it will take 1 billion years for all the stars at the mass that is 
	above the critical mass to be removed, and they are removed at the same constant rate they were formed.
	After each timestep the total luminosity is then calculated by finding the 
	luminosity at each mass and finding the sum of those values.
	
	'''
	for i in range(len(noandmstar)):
		noandmstar[i]=noandmstar[i]+((E_0(mgas=10**11/initialtimestep))*(mass[i]**(-2.35)))
		if criticalmass[j]<=mass[i]:
			noandmstar[i]=noandmstar[i]-((E_0(mgas=10**11/initialtimestep))*(mass[i]**(-2.35)))
		if noandmstar[i]<=0:
			noandmstar[i]=0
		l[i]=luminosity(mass[i])*noandmstar[i]
	iltotal[j]=sum(l)
	'''
	This loop calculates the total V and B magnitude
	
	It works similarly to the luminosity loop. The number of stars formed in each time step for a 
	mass is found and added to the previous timestep. There is then a check to see if this mass is above the 
	critical mass. If it is then a timestep worth of star formation is deducted from the total number. 
	From this the luminosity of the V and B magnitudes from each mass is calculated and multiplied by the respective
	number of stars. The values are then added together to get the total luminosity. This is then converted 
	back to magnitude to find the total magnitude of all the mass ranges at each timestep 
	'''

	for k in range(len(newmass)):
		nonm[k]=nonm[k]+((E_0(mgas=10**11/initialtimestep))*(newmass[k]**(-2.35)))
		if nonm[k]>0:
			if criticalmass[j]<=newmass[k]:
				nonm[k]=nonm[k]-((E_0(mgas=10**11/initialtimestep))*(newmass[k]**(-2.35)))
		lv[k]=lum(vmag[k])*nonm[k]
		lb[k]=lum(bmag[k])*nonm[k]
	vmtotal[j]=(magtotal((sum(lv))))	
	bmtotal[j]=(magtotal((sum(lb))))
	

#The first values in the V and B mag arrays will be 0 due to the time being zero so these values are corrected to nan
vmtotal[0]=np.nan
bmtotal[0]=np.nan
		

'''
								***SIMULATING THE REMAINING 11 BILLION YEARS***
'''
for j in range(timestep2*11):
	'''
	This loop calculates the total luminosity
	
	
	This is similar to the first luminosity loop, however, there is no star formation
	occuring so it is just checking to see if the critical mass is reached at each timestep. 
	If it is then one timestep worth of star formation at that mass range is subtracted. 
	The total luminosity is then calculated by finding the luminosity at each mass and finding
	the sum of those values.
	'''
	
	for i in range(len(noandmstar)):
		if criticalmass2[j]<=mass[i]:
			noandmstar[i]=noandmstar[i]-((E_0(mgas=10**11/timestep2))*(mass[i]**(-2.35)))
			if noandmstar[i]<=0:
				noandmstar[i]=0
		l2[i]=luminosity(mass[i])*noandmstar[i]
	ltotal[j]=sum(l2)

	""""
	This loop calculates the total V and B magnitude
	
	
	This loop finds the total magnitude for the remaining 11 billion years where there is no star formation.
	There is no increase in magnitude so it only checks to see if the critical mass is reached for each mass
	at each timestep. If it is then the number of stars it is reduced by the amount formed at each timestep.
	From this the luminosity of the V and B magnitudes from each mass is found and multiplied by the respective
	number of stars. The values are then added together to get the total luminosity. This is then converted 
	back to magnitude to find the total magnitude of all the mass ranges at each timestep 
	"""
	for k in range(len(newmass)):
		if nonm[k]>0:
			if criticalmass2[j]<=newmass[k]:
				nonm[k]=nonm[k]-((E_0(mgas=10**11/initialtimestep))*(newmass[k]**(-2.35)))
		lv2[k]=lum(vmag[k])*nonm[k]
		lb2[k]=lum(bmag[k])*nonm[k]
	vmtotal2[j]=(magtotal((sum(lv2))))	
	bmtotal2[j]=(magtotal((sum(lb2))))
	

'''
Combining the values found from the first billion years and last 11 billion years
to be able to plot 
'''
#combining the time arrays for the first 1 billion years and last 11 billion years
total_time=np.concatenate((time,time2))
#combining the luminosity arrays for the first 1 billion years and last 11 billion years
total_luminosity=np.concatenate((iltotal,ltotal))
#joining the arrays of total V magnitudes for the first billion years and last 11 billion years
total_vmagnitude=np.concatenate((vmtotal,vmtotal2))
#joining the arrays of total B magnitudes for the first billion years and last 11 billion years
total_bmagnitude=np.concatenate((bmtotal,bmtotal2))

'''
Calculating the B-V colour
'''
#subtracting the total V mag array from the total B mag array to find the total change in B-V colour 
total_bvmagnitude=np.subtract(total_bmagnitude,total_vmagnitude)


'''
The Plots
'''
#plotting change in luminosity as a function of time
plt.figure()
plt.plot(total_time,total_luminosity)
plt.xlabel('Time (Gyrs)', fontsize=15)
plt.ylabel('Luminosity (Solar Luminosity)', fontsize=15)
plt.title('The Change in Luminosity Against Time in a Simulated Galaxy', fontsize=20)
plt.grid()

#plotting change in V mag as a function of time
plt.figure()
plt.plot(total_time,total_vmagnitude)
plt.xlabel('Time (Gyrs)', fontsize=15)
plt.ylabel('Total Absolute V Magnitude', fontsize=15)
plt.title('The Change in Absolute V Magnitude Against Time in a Simulated Galaxy', fontsize=20)
plt.grid()

#plotting change in B-V Colour as a function of time
plt.figure()
plt.plot(total_time,total_bvmagnitude)
plt.xlabel('Time (Gyrs)', fontsize=15)
plt.ylabel('B-V Colour', fontsize=15)
plt.title('The Change in B-V Colour Against Time in a Simulated Galaxy', fontsize=20)
plt.grid()


