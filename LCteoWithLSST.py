#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 10:03:31 2021

@author: silvioleccia
"""

#from numpy import random
import numpy as np
from scipy.interpolate import interp1d


def generateLC(time_lsst,filters_lsst,output_ReadLCTeo,period_true=-99,
               ampl_true=1.,phase_true=0.,do_normalize=False):
    
    u_model=output_ReadLCTeo['u']
    g_model=output_ReadLCTeo['g']
    r_model=output_ReadLCTeo['r']
    i_model=output_ReadLCTeo['i']
    z_model=output_ReadLCTeo['z']
    y_model=output_ReadLCTeo['y']
    phase_model=output_ReadLCTeo['phase']
    
    #Se diamo un periodo arbitrario, usare quello (Caso generico), altrimenti
    #si usa il periodo del modello di Marcella.
    if period_true < -90:
        period_final=(output_ReadLCTeo['period'])/86400.
    else:
        period_final=period_true

    if do_normalize:
        meanu=output_ReadLCTeo['meanu']
        meang=output_ReadLCTeo['meang']
        meanr=output_ReadLCTeo['meanr']
        meani=output_ReadLCTeo['meani']
        meanz=output_ReadLCTeo['meanz']
        meany=output_ReadLCTeo['meany']
        amplg=output_ReadLCTeo['amplg']

        #normalizzo ad ampiezza g
        meanmags=[meanu,meang,meanr,meani,meanz,meany]
        mags_model_norm=[ [(u_model-meanu)/amplg],
                     [(g_model-meang)/amplg],
                     [(r_model-meanr)/amplg],
                     [(i_model-meani)/amplg],
                     [(z_model-meanz)/amplg],
                     [(y_model-meany)/amplg] ]
        model_u=interp1d(phase_model,mags_model_norm[0])
        model_g=interp1d(phase_model,mags_model_norm[1])
        model_r=interp1d(phase_model,mags_model_norm[2])
        model_i=interp1d(phase_model,mags_model_norm[3])
        model_z=interp1d(phase_model,mags_model_norm[4])
        model_y=interp1d(phase_model,mags_model_norm[5])
    else:
        model_u=interp1d(phase_model,u_model)
        model_g=interp1d(phase_model,g_model)
        model_r=interp1d(phase_model,r_model)
        model_i=interp1d(phase_model,i_model)
        model_z=interp1d(phase_model,z_model)
        model_y=interp1d(phase_model,y_model)
        meanmags=[0.,0.,0.,0.,0.,0.]

    allmodels=[ model_u,model_g,model_r,model_i,model_z,model_y ]
    
    t_time_0=min(time_lsst)
    
    ind_u=(np.where(filters_lsst == 'u'))[0]
    ind_g=(np.where(filters_lsst == 'g'))[0]
    ind_r=(np.where(filters_lsst == 'r'))[0]
    ind_i=(np.where(filters_lsst == 'i'))[0]
    ind_z=(np.where(filters_lsst == 'z'))[0]
    ind_y=(np.where(filters_lsst == 'y'))[0]

    timeLSSTu=time_lsst[ind_u]
    timeLSSTg=time_lsst[ind_g]
    timeLSSTr=time_lsst[ind_r]
    timeLSSTi=time_lsst[ind_i]
    timeLSSTz=time_lsst[ind_z]
    timeLSSTy=time_lsst[ind_y]
    
    magLSSTu=np.empty(len(ind_u))
    magLSSTg=np.empty(len(ind_g))
    magLSSTr=np.empty(len(ind_r))
    magLSSTi=np.empty(len(ind_i))
    magLSSTz=np.empty(len(ind_z))
    magLSSTy=np.empty(len(ind_y))
    
    def interpola(timeLSST,meanmags,model):
        magLSST=np.empty(len(timeLSST))
        phaselsst=np.empty(len(timeLSST))
        for i in np.arange(len(timeLSST)):
            phase_lsst_temp=((timeLSST[i]-t_time_0)/period_final) % 1.
            phaselsst[i]=phase_lsst_temp
            magLSST[i]=meanmags+ampl_true*model(phase_lsst_temp)
        return phaselsst,magLSST

    phaselsst_u,magLSSTu=interpola(timeLSSTu,meanmags[0],model_u)
    phaselsst_g,magLSSTg=interpola(timeLSSTg,meanmags[1],model_g)
    phaselsst_r,magLSSTr=interpola(timeLSSTr,meanmags[2],model_r)
    phaselsst_i,magLSSTi=interpola(timeLSSTi,meanmags[3],model_i)
    phaselsst_z,magLSSTz=interpola(timeLSSTz,meanmags[4],model_z)
    phaselsst_y,magLSSTy=interpola(timeLSSTy,meanmags[5],model_y)
    
    mag_all=np.empty(len(time_lsst))
    phase_all=np.empty(len(time_lsst))
    time_all=np.empty(len(time_lsst))

    #mag_all è ordinato come time_lsst    
    mag_all[ind_u]=magLSSTu
    mag_all[ind_g]=magLSSTg
    mag_all[ind_r]=magLSSTr
    mag_all[ind_i]=magLSSTi
    mag_all[ind_z]=magLSSTz
    mag_all[ind_y]=magLSSTy

    phase_all[ind_u]=phaselsst_u
    phase_all[ind_g]=phaselsst_g
    phase_all[ind_r]=phaselsst_r
    phase_all[ind_i]=phaselsst_i
    phase_all[ind_z]=phaselsst_z
    phase_all[ind_y]=phaselsst_y

    time_all[ind_u]=timeLSSTu
    time_all[ind_g]=timeLSSTg
    time_all[ind_r]=timeLSSTr
    time_all[ind_i]=timeLSSTi
    time_all[ind_z]=timeLSSTz
    time_all[ind_y]=timeLSSTy
        
    return {'timeu':timeLSSTu,'timeg':timeLSSTg,
                'timer':timeLSSTr,'timei':timeLSSTi,'timez':timeLSSTz,
                'timey':timeLSSTy,'magu':magLSSTu,'magg':magLSSTg,
                'magr':magLSSTr,'magi':magLSSTi,
                'magz':magLSSTz,'magy':magLSSTy,
           'phaseu':phaselsst_u,'phaseg':phaselsst_g,'phaser':phaselsst_r,
           'phasei':phaselsst_i,'phasez':phaselsst_z,'phasey':phaselsst_y,
           'mag_all':mag_all,'phase_all':phase_all,'time_all':time_all,
           'ind_u':ind_u,'ind_g':ind_g,'ind_r':ind_r,
           'ind_i':ind_i,'ind_z':ind_z,'ind_y':ind_y}

    
def noising(LcTeoLSST,snr,sigma):
    
#noising 
    def noisingBand(timeLSSTteo,magLSSTteo,snr,sigma):
        magNoised=[]
        for j in range(len(timeLSSTteo)):            
            dmag = 2.5*np.log10(1.+1./snr[j])
            noise = np.random.randint(-sigma,sigma)*dmag
            magNoisedComp=magLSSTteo[j]+noise
            magNoised.append(magNoisedComp)
    
        return magNoised, noise ,dmag

    magNoisedu,noiseu,dmagu=noisingBand(LcTeoLSST['timeu'],LcTeoLSST['magu'],snr['u'],sigma)
    magNoisedg,noiseg,dmagg=noisingBand(LcTeoLSST['timeg'],LcTeoLSST['magg'],snr['g'],sigma)
    magNoisedr,noiser,dmagr=noisingBand(LcTeoLSST['timer'],LcTeoLSST['magr'],snr['r'],sigma)
    magNoisedi,noisei,dmagi=noisingBand(LcTeoLSST['timei'],LcTeoLSST['magi'],snr['i'],sigma)
    magNoisedz,noisez,dmagz=noisingBand(LcTeoLSST['timez'],LcTeoLSST['magz'],snr['z'],sigma)
    magNoisedy,noisey,dmagy=noisingBand(LcTeoLSST['timey'],LcTeoLSST['magy'],snr['y'],sigma)
    
    #mag_all è ordinato come time_lsst
    mag_all=np.empty(len(LcTeoLSST['mag_all']))
    mag_all[LcTeoLSST['ind_u']]=magNoisedu
    mag_all[LcTeoLSST['ind_g']]=magNoisedg
    mag_all[LcTeoLSST['ind_r']]=magNoisedr
    mag_all[LcTeoLSST['ind_i']]=magNoisedi
    mag_all[LcTeoLSST['ind_z']]=magNoisedz
    mag_all[LcTeoLSST['ind_y']]=magNoisedy
    
    
    #noise_all è ordinato come time_lsst
    noise_all=np.empty(len(LcTeoLSST['mag_all']))
    noise_all[LcTeoLSST['ind_u']]=noiseu
    noise_all[LcTeoLSST['ind_g']]=noiseg
    noise_all[LcTeoLSST['ind_r']]=noiser
    noise_all[LcTeoLSST['ind_i']]=noisei
    noise_all[LcTeoLSST['ind_z']]=noisez
    noise_all[LcTeoLSST['ind_y']]=noisey
    
   #mag_all è ordinato come time_lsst
    dmag_all=np.empty(len(LcTeoLSST['mag_all']))
    dmag_all[LcTeoLSST['ind_u']]=dmagu
    dmag_all[LcTeoLSST['ind_g']]=dmagg
    dmag_all[LcTeoLSST['ind_r']]=dmagr
    dmag_all[LcTeoLSST['ind_i']]=dmagi
    dmag_all[LcTeoLSST['ind_z']]=dmagz
    dmag_all[LcTeoLSST['ind_y']]=dmagy
    
    output2={'timeu':LcTeoLSST['timeu'],'timeg':LcTeoLSST['timeg'],
                'timer':LcTeoLSST['timer'],'timei':LcTeoLSST['timei'],
                'timez':LcTeoLSST['timez'],
                'timey':LcTeoLSST['timey'],'magu':np.asarray(magNoisedu),'magg':np.asarray(magNoisedg),
                'magr':np.asarray(magNoisedr),'magi':np.asarray(magNoisedi),
                'magz':np.asarray(magNoisedz),'magy':np.asarray(magNoisedy),
            'phaseu':LcTeoLSST['phaseu'],'phaseg':LcTeoLSST['phaseg'],'phaser':LcTeoLSST['phaser'],
           'phasei':LcTeoLSST['phasei'],'phasez':LcTeoLSST['phasez'],'phasey':LcTeoLSST['phasey'],
            'mag_all':mag_all, 'time_all':LcTeoLSST['time_all'],'noise_all':noise_all,'dmag_all':dmag_all}
    return output2

