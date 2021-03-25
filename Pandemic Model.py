# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 17:42:53 2020

@author: Ed
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sc

'''For Question 3 I defined a function where the fourth order Runge-Kutta(R-K) method was used
The method I used was extrapolated from a second order R-K method 
(SOURCE:https://www.mdpi.com/2227-7390/8/7/1174/htm)
It works by first defining all the variables used inthe equations provided in the question
From there it loops it for 300 days, and as my change in time is 0.1, it therefore runs 3000 times.
All the important values are then added after each loop to the empty lists created before the loop.
I then return these lists to be able to plot for the remainder of the question
'''



def sir4ork(p=0,N=1,r=0,beta=0.14,gamma=1/14,dt=0.1,t = 0):
	#defining the equation for suseptability, including p for question 4
	i=N*0.001
	s=N-(1*i)-p

	time=[]
	suseptable=[]
	infected=[]
	recovered=[]
	#the for loop will loop this for 300 days. index is 30000 because dt
	for index in range(3000):
		# step 1
		ks1 = - ((beta*s*i)/N)
		ki1 = ((beta*s*i)/N) - (gamma*i)
		# step 2
		s2 = s + (ks1*dt)
		i2 = i + (ki1*dt)
		ks2 = - (beta*s2*i2)
		ki2 = (beta*s2*i2) - (gamma*i2)
		# step 3   
		s3 = s2 + (ks2*dt)
		i3 = i2 + (ki2*dt)
		ks3 = - (beta*s3*i3)
		ki3 = (beta*s3*i3) - (gamma*i3)
		#step 4   
		s4 = (s3 + ks1*dt)
		i4 = (i3 + ki1*dt)
		ks4 = - (beta*s4*i4)
		ki4 = (beta*s4*i4) - (gamma*i4)
		#the final values for each time period (dt)    
		s = s + (ks1+(2*ks2)+(2*ks3)+ks4)*dt/6
		i = i + (ki1+(2*ki2)+(2*ki3)+ki4)*dt/6
		r = N - s - i 
		t = t + dt
		#Adding these values to the lists created earlier   
		time.append(t)
		suseptable.append(s)
		infected.append(i)
		recovered.append(r)
	return infected,time,suseptable,recovered

infected,time,suseptable,recovered=sir4ork()
#Plotting the results
plt.figure(1)
plt.plot(time,np.array(suseptable)*100, 'r', linewidth=1, label='S')
plt.plot(time,np.array(infected)*100, 'b', linewidth=1, label='I')
plt.plot(time,np.array(recovered)*100, 'g', linewidth=1, label='R')
plt.title("Runge - Kutta SIR Model")
plt.legend()
plt.xlabel('Time (Days)')
plt.ylabel('Population as a Percentage')
plt.grid()
plt.show()







plt.figure(2)
plt.grid()
for i in range(0,6,1):
	i=i/10
	infected,timem,suseptable,recovered=sir4ork(p=i)
	plt.plot(time,np.array(infected)*100, linewidth=1, label=i*100)
	plt.title('How Different Percentages of Population Vaccinated Effects Infection Levels')
	plt.legend()
	





plt.figure(3)

for i in range(0,6,1):
	i=i/10
	infected,timem,suseptable,recovered=sir4ork(p=i)
	hospitalised=0.05*np.array(infected)
	plt.plot(time,hospitalised*100, linewidth=1, label=i*100)
	plt.title('How Different Percentages of Population Vaccinated Effects Infection Levels')
	
	
plt.axhline(y=0.25, linestyle='--', c='black', label='hospital beds available' )
plt.legend()
plt.grid()

plt.figure(4,figsize=(8,6))
for i in range(20,25,1):
	i=i/100
	infected,timem,suseptable,recovered=sir4ork(p=i)
	hospitalised=0.05*np.array(infected)
	plt.plot(time,hospitalised*100, linewidth=1, label=i*100)
	plt.title('How Different Percentages of Population Vaccinated Effects Infection Levels')

plt.axhline(y=0.25, linestyle='--', c='black', label='hospital beds available' )
plt.legend()
plt.grid()


plt.figure(5)
def trap(yfinal,yinitial):
	h=1
	return ((h/2)*(yfinal+yinitial))


def culm_trap():
    y = (1/14)

    for i in range(0,6,1):
        i=i/10
        infected,time,suseptable,recovered=sir4ork(i)
        infectedgamma = np.array(infected) * y
        vals=[]
        
        for index in range(len(infectedgamma)-1):
            if index<len(infectedgamma)-1:
                holder=trap(infectedgamma[index],infectedgamma[index+1])
                vals.append(holder)
        
        #Values found for areas, now to plot a cumulative sum against t
        plt.plot(time[:-1], np.cumsum(vals)*10)
        
    #Format plot
    plt.legend(['p=0%','p=10%','p=20%','p=30%','p=40%','p=50%'])
    plt.xlabel('Days')
    plt.ylabel('Percentage of Population')
    plt.title('Culmative Sum of The Population Infected')
    plt.grid()
    return (vals)	
    
    
    
#Call function
culm_trap()
