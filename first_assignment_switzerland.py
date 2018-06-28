# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 17:12:35 2017

@author: leoIE
"""
import math
import scipy.stats
from matplotlib import pyplot as plt
import numpy as np
import matplotlib
import matplotlib.ticker as ticker

############################################################################################
"""Open file and read from it the contents"""
with open('population_switzerland.dat','r') as fin:
    header=fin.readline()
    pop={}
    #print (content)
    for line in fin:
        A=line.strip().split('\t')
        pop[int(A[0])]=int(A[1])
        
############################################################################################
"""Growth rate"""
def GR(pop):
    temp=[]
    interval=[1,1700,1900,1950,2000,2009]
    for i in range(len(interval)-1): 
        gr1=float((pop[interval[i+1]]-pop[interval[i]])/(pop[interval[i]]*(interval[i+1]-interval[i])))   
        temp.append(gr1*100)
    return temp,interval

############################################################################################
"""Population change"""        
def pop_change(pop):
    ordered=sorted(pop)    
    temp=[]
    for i in range(len(ordered)-1): 
        A=pop[ordered[i+1]]-pop[ordered[i]]
        temp.append(A)        
    return temp     
   
############################################################################################    
"""Growth rate"""
#Pop_change=pop_change(pop)  
Growth_rate,interval=GR(pop)     
with open("growth_rate_switzerland.dat",'w') as fout:
    fout.write("Period (Years)\tGrowth rate (%)\t\n")

    i=0
    for item in Growth_rate:
        #print (item)
        fout.write("%d-%d\t%g\n"%(interval[i],interval[i+1],item))
        i+=1
        
############################################################################################    
"""Read keys and values from population dictionary"""        
population=[value for value in pop.values()]
years=[key for key in pop.keys()]

############################################################################################
"""Fitting models"""
#Coefficient of determination
def r(obs,sim):
    slope, intercept, r_value, p_value, std_err = stats.linregress(obs,sim)
    return r_value**2
    

#Percentage bias
def pbias(obs,sim):
    temp=0
    for item in obs.keys():
        #print (item)
        temp1=sim[item]-obs[item]
        temp=temp1+temp
    return (100.0*(temp/sum(obs))) 

############################################################################################
"""Hubbert"""
def Hubbert(K,r):
    Hub={}
    j=2030
    i=r
    for k in range(1,2500):
        tempo=300000+K/(1+math.exp(i*(j-k)))
        Hub[k]=tempo
    return Hub

K=16833900
r=0.01151871
H=Hubbert(K,r)

############################################################################################
"""Calculating the coefficient of determination and PBIAS"""
percentage_bias=pbias(pop,H)

#Taking only the keys in the original pop dictionary
H_r=[]
for item in years:
    H_r.append(H[item])
    
r=scipy.stats.pearsonr(population,H_r)
R=r[0]**2

############################################################################################
"""Comparing observed data with Hubbert's model"""
Delta=[]
for year in years:
    Delta.append((pop[year]-H[year]))
    
############################################################################################
"""Plotting the results"""
fig,ax=plt.subplots()
ax.plot(range(1,2500),H.values(),label='Hubbert')
ax.scatter(years,population,5,'r',label='Observed')
ax.set(xlabel='Years', ylabel='Population (Millions)',
       title='Observed population data VS Fitted model')
scale_y = 1e6
ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_y))
ax.yaxis.set_major_formatter(ticks_y)
ax.grid()
plt.legend(loc=2)
ax.set_xlim([0,2009])
plt.show()
fig.savefig("population_switzerland.png")


"""Plotting the difference between observer data and Hubbert's model"""

fig,ax=plt.subplots(1,2)
ax[1].scatter(years,Delta,5,label='Difference',color='g')
ax[1].set(xlabel='Years', ylabel='Population difference (Millions)')
scale_y = 1e6
ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_y))
ax[1].yaxis.set_major_formatter(ticks_y)
ax[1].grid()
ax[1].set_xlim([1800,2100])
ax[1].legend(loc=0)
ax[1].yaxis.tick_right()
ax[1].yaxis.set_label_position("right")

ax[0].scatter(years,population,5,label='Observed',color='r')
ax[0].plot(range(1,2500),H.values(),label='Hubbert')
ax[0].set(xlabel='Years',ylabel='Population (Millions)')
scale_y = 1e6
ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_y))
ax[0].yaxis.set_major_formatter(ticks_y)
ax[0].grid()
ax[0].set_xlim([1800,2100])
fig.suptitle('Difference between observed data and Hubbert')
ax[0].legend(loc=0)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
fig.savefig("Difference plot.png")
plt.show()























