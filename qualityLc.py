#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 09:50:48 2021

@author: silvioleccia
"""

import numpy as np
from scipy import stats as stat
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

def modelToFit(time,param):
    #time in days
    magModel=[]
    amplitudes=[]
    phases=[]
    numberOfHarmonics=int((len(param)-2)/2)
    print(numberOfHarmonics)
    zp=param[1]
    period=param[0]
    for i in range(numberOfHarmonics):
        amplitudes.append(param[2+i])
        phases.append(param[numberOfHarmonics+2+i])
#    print(numberOfHarmonics)
    for i in range(len(time)):
        y=zp
        for j in range(0,int(numberOfHarmonics)):
            y=y+amplitudes[j]*np.cos((2*np.pi/period)*(j+1)*(time[i])+phases[j])
        magModel.append(y)
    print(amplitudes)
    print(phases)
    return magModel

def qualityCheck(time,period,factor1,yfin):
    #period=param[0]
    phase= ((time-time[0])/period)%1
    indexSorted=np.argsort(phase)
   
    distances=[]
    indexStart=[]
    indexStop=[]
    leftDistance=phase[indexSorted[0]]
    rightDistance=1-phase[indexSorted[len(indexSorted)-1]]
    for i in range(len(phase)-1):
        dist=phase[indexSorted[i+1]]-phase[indexSorted[i]]
        distances.append(dist)
        
        
    #factor=sum(distances)/len(distances)*factor1
    distancesTotal=distances
    distancesTotal.append(leftDistance)
    distancesTotal.append(rightDistance)
    factor=sum(distancesTotal)/len(distancesTotal)*factor1
    maxDistance=max(distancesTotal)
    for i in range(len(phase)-1):
        dist=phase[indexSorted[i+1]]-phase[indexSorted[i]]
        distances.append(dist)
        if (dist > factor):
            indexStart.append(indexSorted[i])
            indexStop.append(indexSorted[i+1])
            
    return maxDistance,len(indexStart)#,indexStart,indexStop



"""        
time1=np.arange(0,10,2.1)
time2=np.arange(14,30,1.1)
time3=np.arange(35,40,0.2)

time=np.concatenate((np.asarray(time1),np.asarray(time2),np.asarray(time3)))
param=[3.2,3,1,1]
mag=modelToFit(time,param)

a,b=qualityCheck(time,mag,param[0],2)



phaseComputed=(((time-time[0])/param[0])%1)

fig=plt.figure()
ax1 = plt.subplot(111)  
ax1.set_ylabel('u (mag)')
ax1.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
ax1.invert_yaxis()
#ax1.set_xlim([0,1])
    #ax2.plot(phaseGAIA1,magGAIA1,'o',color='red')
ax1.plot(phaseComputed,mag,'o',color='black')      
"""