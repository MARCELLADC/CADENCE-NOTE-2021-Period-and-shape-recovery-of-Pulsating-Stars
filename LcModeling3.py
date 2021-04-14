#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 15:28:05 2021

@author: silvioleccia
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 10:14:13 2021

@author: silvioleccia
"""
import numpy as np
from scipy import stats as stat
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
from astropy.io import ascii
#from astropy.table import Table
def modelToFit(time,param):
    #time in days
    magModel=[]
    amplitudes=[]
    phases=[]
    numberOfHarmonics=int((len(param)-2)/2)
    
    zp=param[1]
    period=param[0]
    for i in range(numberOfHarmonics):
        amplitudes.append(param[2+i])
        phases.append(param[numberOfHarmonics+2+i])
    for i in range(len(time)):
        y=zp
        for j in range(0,int(numberOfHarmonics)):
            y=y+amplitudes[j]*np.cos((2*np.pi/period)*(j+1)*(time[i])+phases[j])
        magModel.append(y)
    return magModel

def modelToFit2(phase,param):
    #time in days
    magModel=[]
    amplitudes=[]
    phases=[]
    numberOfHarmonics=int((len(param)-2)/2)
    
    zp=param[1]
    period=param[0]
    for i in range(numberOfHarmonics):
        amplitudes.append(param[2+i])
        phases.append(param[numberOfHarmonics+2+i])
    for i in range(len(phase)):
        y=zp
        for j in range(0,int(numberOfHarmonics)):
            y=y+amplitudes[j]*np.cos(np.pi*(j+1)*(phase[i])+phases[j])
        magModel.append(y)
    return magModel
def chisqr(residual,Ndat,Nvariable):
    chi=sum(pow(residual,2))/(Ndat-Nvariable)
    return chi

def residuals(datax,datay,fitparameters):
    residuals=modelToFit(datax,fitparameters)-datay
    return residuals


def meanmag_antilog(mag):
    mag=np.asarray(mag)
    flux=10.**(-mag/2.5)
    if len(flux)>0:
        result=(-2.5)*np.log10(sum(flux)/len(flux))
    else:
        result=9999.
    return result

def computingLcModel(data,period,numberOfHarmonics,index,outDir):
    

    def modelToFit2_fit(coeff):
        fit = modelToFit(x,coeff)
        return (fit - y_proc)
    
    
#    time_u=data['timeu']#time must be in days
#    time_g=data['timeg']
#    time_r=data['timer']
#    time_i=data['timei']
#    time_z=data['timez']
#    time_y=data['timey']
#    mag_u=data['magu']
#    mag_g=data['magg']
#    mag_r=data['magr']
#    mag_i=data['magi']
#    mag_z=data['magz']

    time_u=data['timeu'][index['ind_notsaturated_u']]#time must be in days
    time_g=data['timeg'][index['ind_notsaturated_g']]
    time_r=data['timer'][index['ind_notsaturated_r']]
    time_i=data['timei'][index['ind_notsaturated_i']]
    time_z=data['timez'][index['ind_notsaturated_z']]
    time_y=data['timey'][index['ind_notsaturated_y']]
    mag_u=data['magu'][index['ind_notsaturated_u']]
    mag_g=data['magg'][index['ind_notsaturated_g']]
    mag_r=data['magr'][index['ind_notsaturated_r']]
    mag_i=data['magi'][index['ind_notsaturated_i']]
    mag_z=data['magz'][index['ind_notsaturated_z']]
    mag_y=data['magy'][index['ind_notsaturated_y']] 

    
    parametersForLcFit=[period,1] # period,zp
    for i in range(numberOfHarmonics):
        parametersForLcFit.append(1)#added ampl
        parametersForLcFit.append(1)#added phase
    x=time_u
    
    
    
    y_proc = np.copy(mag_u)
    if len(y_proc)>(numberOfHarmonics*2)+2:
        print('fitting u band')
        fit_u,a=leastsq( modelToFit2_fit,parametersForLcFit)
        residual=residuals(x,y_proc,fit_u)
        chi_u=chisqr(residual,len(x),len(fit_u))
    else:
        fit_u=[9999.]
        chi_u=9999.
    x=time_g
    y_proc = np.copy(mag_g)
    if len(y_proc)>(numberOfHarmonics*2)+2:
        print('fitting g band')
        fit_g,a=leastsq( modelToFit2_fit,parametersForLcFit)
        residual=residuals(x,y_proc,fit_g)
        chi_g=chisqr(residual,len(x),len(fit_g))
    else:
        fit_g=[9999.]
        chi_g=9999.
    y_proc = np.copy(mag_r)
    x=time_r
    if len(y_proc)>(numberOfHarmonics*2)+2:
        print('fitting r band')
        fit_r,a=leastsq( modelToFit2_fit,parametersForLcFit)
        residual=residuals(x,y_proc,fit_r)
        chi_r=chisqr(residual,len(x),len(fit_r))
    else:
        fit_r=[9999.]
        chi_r=9999.
    x=time_i
    y_proc = np.copy(mag_i)
    if len(y_proc)>(numberOfHarmonics*2)+2:
        print('fitting i band')
        fit_i,a=leastsq( modelToFit2_fit,parametersForLcFit)
        residual=residuals(x,y_proc,fit_i)
        chi_i=chisqr(residual,len(x),len(fit_i))
    else:
        fit_i=[9999.]
        chi_i=9999.
    x=time_z
    y_proc = np.copy(mag_z)
    if len(y_proc)>(numberOfHarmonics*2)+2:
        print('fitting z band')
        fit_z,a=leastsq( modelToFit2_fit,parametersForLcFit)
        residual=residuals(x,y_proc,fit_z)
        chi_z=chisqr(residual,len(x),len(fit_z))
    else:
        fit_z=[9999.]
        chi_z=9999.
    x=time_y
    y_proc = np.copy(mag_y)
    if len(y_proc)>(numberOfHarmonics*2)+2:
        print('fitting y band')
        fit_y,a=leastsq( modelToFit2_fit,parametersForLcFit)
        residual=residuals(x,y_proc,fit_y)
        chi_y=chisqr(residual,len(x),len(fit_y))
    else:
        fit_y=[9999.]
        chi_y=9999.
        
    results={'u':fit_u,'g':fit_g,'r':fit_r,'i':fit_i,'z':fit_z,'y':fit_y,'chi_u':chi_u,
             'chi_g':chi_g,'chi_r':chi_r,'chi_i':chi_i,'chi_z':chi_z,'chi_y':chi_y}
    
    return results


def plotting(data,fittingParameters,period,yearfin,outDir):
    phaseforModel=np.arange(0,1,0.001)
    fig=plt.figure(figsize=(10,16), dpi=80)
    plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95
                        , left = 0.1, hspace = 0.08, wspace = 0.2)
    
    ax1 = plt.subplot2grid((3,2), (0,0)) # topleft    
    ax1.set_ylabel('u (mag)')
    ax1.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax1.invert_yaxis()
    ax1.set_xlim([0,1])
    ax1.plot(((data['timeu'])/period)%1,data['magu'],'o',color='purple')
    if len(fittingParameters['u'])>1:
        timeForModel=np.arange(data['timeu'][0],data['timeu'][0]+2*period,0.01)
        magModelForPlot=modelToFit(timeForModel,fittingParameters['u'])
        ax1.plot((timeForModel/period)%1,magModelForPlot,'.',color='black')       
    ax1.set_xlabel('phase')
    ax2 = plt.subplot2grid((3,2), (1,0)) # topleft    
    ax2.set_ylabel('g (mag)')
    ax2.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax2.invert_yaxis()
    ax2.set_xlim([0,1])
    ax2.plot(((data['timeg'])/period)%1,data['magg'],'o',color='purple')
    if len(fittingParameters['g'])>1:
        timeForModel=np.arange(data['timeg'][0],data['timeg'][0]+2*period,0.01)
        magModelForPlot=modelToFit(timeForModel,fittingParameters['g'])
        ax2.plot((timeForModel/period)%1,magModelForPlot,'.',color='black')
    ax2.set_xlabel('phase')
    ax3 = plt.subplot2grid((3,2), (2,0)) # topleft    
    ax3.set_ylabel('r (mag)')
    ax3.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax3.invert_yaxis()
    ax3.set_xlim([0,1])
    ax3.plot(((data['timer'])/period)%1,data['magr'],'o',color='purple')
    if len(fittingParameters['r'])>1:
        timeForModel=np.arange(data['timer'][0],data['timer'][0]+2*period,0.01)
        magModelForPlot=modelToFit(timeForModel,fittingParameters['r'])
        ax3.plot((timeForModel/period)%1,magModelForPlot,'.',color='black')
    ax3.set_xlabel('phase')
    ax4 = plt.subplot2grid((3,2), (0,1)) # topleft    
    ax4.set_ylabel('i (mag)')
    ax4.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax4.invert_yaxis()
    ax4.set_xlim([0,1])
    ax4.plot(((data['timei'])/period) %1,data['magi'],'o',color='purple')
    if len(fittingParameters['i'])>1:
        timeForModel=np.arange(data['timei'][0],data['timei'][0]+2*period,0.01)
        magModelForPlot=modelToFit(timeForModel,fittingParameters['i'])
        ax4.plot((timeForModel/period)%1,magModelForPlot,'.',color='black')
    ax4.set_xlabel('phase')
    ax5 = plt.subplot2grid((3,2), (1,1)) # topleft    
    ax5.set_ylabel('z (mag)')
    ax5.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax5.invert_yaxis()
    ax5.set_xlim([0,1])
    ax5.plot(((data['timez'])/period) %1,data['magz'],'o',color='purple')
    if len(fittingParameters['z'])>1:
        timeForModel=np.arange(data['timez'][0],data['timez'][0]+2*period,0.01)
        magModelForPlot=modelToFit(timeForModel,fittingParameters['z'])
        ax5.plot((timeForModel/period)%1,magModelForPlot,'.',color='black')
    ax5.set_xlabel('phase')
    ax6 = plt.subplot2grid((3,2), (2,1)) # topleft    
    ax6.set_ylabel('y (mag)')
    ax6.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax6.invert_yaxis()
    ax6.set_xlim([0,1])
    ax6.plot(((data['timey'])/period) %1,data['magy'],'o',color='purple')
    if len(fittingParameters['y'])>1:
        timeForModel=np.arange(data['timey'][0],data['timey'][0]+2*period,0.01)
        magModelForPlot=modelToFit(timeForModel,fittingParameters['y'])
        ax6.plot((timeForModel/period)%1,magModelForPlot,'.',color='black')
    ax6.set_xlabel('phase')
    
    plt.savefig(str(outDir)+'/Fit'+str(yearfin)+'.pdf')
def computation(data,period,numberOfHarmonics,dataTeo,yearfin,index,outDir):
    print('fitting...')
    fitting=computingLcModel(data,period,numberOfHarmonics,index,outDir)
    timeForModel=np.arange(data['timeu'][0],data['timeu'][0]+2*period,0.01)
    #computing the magModelFromFit
    if len(fitting['u'])>1:
        magModelFromFit_u=modelToFit(timeForModel,fitting['u']) 
        ampl_u=max(magModelFromFit_u)-min(magModelFromFit_u)
    else:
        magModelFromFit_u=[9999.]
        ampl_u=9999.
    timeForModel=np.arange(data['timeg'][0],data['timeg'][0]+2*period,0.01)
    if len(fitting['g'])>1:
        #magModelFromFit_g=modelToFit(data['timeg'],fitting['g']) 
        magModelFromFit_g=modelToFit(timeForModel,fitting['g'])
        ampl_g=max(magModelFromFit_g)-min(magModelFromFit_g)
    else:
        magModelFromFit_g=[9999.]
        ampl_g=9999.
    timeForModel=np.arange(data['timer'][0],data['timer'][0]+2*period,0.01)
    if len(fitting['r'])>1:
        #magModelFromFit_r=modelToFit(data['timer'],fitting['r']) 
        magModelFromFit_r=modelToFit(timeForModel,fitting['r'])
        ampl_r=max(magModelFromFit_r)-min(magModelFromFit_r)
    else:
        magModelFromFit_r=[9999.]
        ampl_r=9999.
    timeForModel=np.arange(data['timei'][0],data['timei'][0]+2*period,0.01)
    
    if len(fitting['i'])>1:
        
        magModelFromFit_i=modelToFit(timeForModel,fitting['i'])
        
        if len(magModelFromFit_i)>0:
            ampl_i=max(magModelFromFit_i)-min(magModelFromFit_i)
        else:
            ampl_i=9999.
    else:
        magModelFromFit_i=[9999.]
        ampl_i=9999.
    timeForModel=np.arange(data['timez'][0],data['timez'][0]+2*period,0.01)    
    if len(fitting['z'])>1:
        
        magModelFromFit_z=modelToFit(timeForModel,fitting['z'])
        ampl_z=max(magModelFromFit_z)-min(magModelFromFit_z)
    else:
        magModelFromFit_z=[9999.]
        ampl_z=9999.
    timeForModel=np.arange(data['timey'][0],data['timey'][0]+2*period,0.01)    
    if len(fitting['y'])>1:
        
        magModelFromFit_y=modelToFit(timeForModel,fitting['y']) 
        ampl_y=max(magModelFromFit_y)-min(magModelFromFit_y)
    else:
        magModelFromFit_y=[9999.]
        ampl_y=9999.
    
   
    
    meanMag_u=meanmag_antilog(magModelFromFit_u)
    meanMag_g=meanmag_antilog(magModelFromFit_g)
    meanMag_r=meanmag_antilog(magModelFromFit_r)
    meanMag_i=meanmag_antilog(magModelFromFit_i)
    meanMag_z=meanmag_antilog(magModelFromFit_z)
    meanMag_y=meanmag_antilog(magModelFromFit_y)
    
   
    
    ampl_u=max(magModelFromFit_u)-min(magModelFromFit_u)
    ampl_g=max(magModelFromFit_g)-min(magModelFromFit_g)
    ampl_r=max(magModelFromFit_r)-min(magModelFromFit_r)
    ampl_i=max(magModelFromFit_i)-min(magModelFromFit_i)
    ampl_z=max(magModelFromFit_z)-min(magModelFromFit_z)
    ampl_y=max(magModelFromFit_y)-min(magModelFromFit_y)
    
    skew_u=stat.skew(magModelFromFit_u)#?
    skew_g=stat.skew(magModelFromFit_g)#?
    skew_r=stat.skew(magModelFromFit_r)#?
    skew_i=stat.skew(magModelFromFit_i)#?
    skew_z=stat.skew(magModelFromFit_z)#?
    skew_y=stat.skew(magModelFromFit_y)#?
    
    DeltaMean_u= meanMag_u-dataTeo['meanu']
    DeltaMean_g= meanMag_g-dataTeo['meang']
    DeltaMean_r= meanMag_r-dataTeo['meanr']
    DeltaMean_i= meanMag_i-dataTeo['meani']
    DeltaMean_z= meanMag_z-dataTeo['meanz']
    DeltaMean_y= meanMag_y-dataTeo['meany']
    
    DeltaAmpl_u=ampl_u-dataTeo['amplu']
    DeltaAmpl_g=ampl_g-dataTeo['amplg']
    DeltaAmpl_r=ampl_r-dataTeo['amplr']
    DeltaAmpl_i=ampl_i-dataTeo['ampli']
    DeltaAmpl_z=ampl_z-dataTeo['amplz']
    DeltaAmpl_y=ampl_y-dataTeo['amply']
    
    DeltaSkew_u=skew_u-skew_u #da cambiare quando calcolata la skew per bene
    DeltaSkew_g=skew_g-skew_g #da cambiare quando calcolata la skew per bene
    DeltaSkew_r=skew_r-skew_r #da cambiare quando calcolata la skew per bene
    DeltaSkew_i=skew_i-skew_i #da cambiare quando calcolata la skew per bene
    DeltaSkew_z=skew_z-skew_z #da cambiare quando calcolata la skew per bene
    DeltaSkew_y=skew_y-skew_y #da cambiare quando calcolata la skew per bene
    
    
    plotting(data,fitting,period,yearfin,outDir)
    
    
    finalResult={'DeltaMean_u':DeltaMean_u,'DeltaMean_g':DeltaMean_g,'DeltaMean_r':DeltaMean_r,
                 'DeltaMean_i':DeltaMean_i,'DeltaMean_z':DeltaMean_z,'DeltaMean_y':DeltaMean_y,
                 'DeltaAmpl_u':DeltaAmpl_u,'DeltaAmpl_g':DeltaAmpl_g,'DeltaAmpl_r':DeltaAmpl_r,
                 'DeltaAmpl_i':DeltaAmpl_i,'DeltaAmpl_z':DeltaAmpl_z,'DeltaAmpl_y':DeltaAmpl_y,
                 'DeltaSkew_u':DeltaSkew_u,'DeltaSkew_g':DeltaSkew_g,'DeltaSkew_r':DeltaSkew_r,
                 'DeltaSkew_i':DeltaSkew_i,'DeltaSkew_z':DeltaSkew_z,'DeltaSkew_y':DeltaSkew_y,
                 'chi_u':fitting['chi_u'],'chi_g':fitting['chi_g'],'chi_r':fitting['chi_r'],
                 'chi_i':fitting['chi_i'],'chi_z':fitting['chi_z'],'chi_y':fitting['chi_y'],
                  'fittingParametersAllband':fitting}
    
    return finalResult




