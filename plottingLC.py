#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 12:32:33 2021

@author: silvioleccia
"""

from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import numpy as np
from pandas import DataFrame
import os
from os import path
from astropy.table import Table
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker


def magAndPhaseMod2(zp, frequency,tref,amplitudes, phases):
    
    
    timet=range(int(tref),int(tref*10))
    #finalPhases=np.empty(len(modelPhases))
    #for i in range(len(modelPhases)):
    #    finalPhases[i]=modelPhases[i]/1000
    
    phase=[]
    for i in range(len(timet)):
        phase.append((timet[i]-tref)*frequency-int((timet[i]-tref)*frequency))
    magModel=[]
    for j in range(len(timet)):
        magModelComp=zp
        for i in range(len(amplitudes)):
            magModelComp=magModelComp+amplitudes[i]*np.cos(2*np.pi*(i+1)*frequency*(timet[j]-tref)+phases[i])
           
        #print(magModelComp, j)    
        magModel.append(magModelComp)
   
    
    
    return phase,magModel 

def magAndPhaseMod3(time, zp, frequency,tref,amplitudes, phases):
    
    
    timet=time
    #finalPhases=np.empty(len(modelPhases))
    #for i in range(len(modelPhases)):
    #    finalPhases[i]=modelPhases[i]/1000
    
    phase=[]
    for i in range(len(timet)):
        phase.append((timet[i]-tref)*frequency-int((timet[i]-tref)*frequency))
        
    magModel=[]
    for j in range(len(timet)):
        magModelComp=zp
        for i in range(len(amplitudes)):
            magModelComp=magModelComp+amplitudes[i]*np.cos(2*np.pi*(i+1)*frequency*(timet[j]-tref)+phases[i])
           
        #print(magModelComp, j)    
        magModel.append(magModelComp)
   
    
    
    return phase,magModel 

def plotting_LSST(title,phaseModelu,magModelu,phaseModelg,magModelg,
              phaseModelr,magModelr,phaseModeli,magModeli,phaseModelz,magModelz,phaseModely,magModely,
              phaseLSSTu,magLSSTu,phaseLSSTg,magLSSTg,phaseLSSTr,magLSSTr,
              phaseLSSTi,magLSSTi,phaseLSSTz,magLSSTz,phaseLSSTy,magLSSTy,
                  snru,snrg,snrr,snri,snrz,snry):
    
#    phaseLSSTAll=np.concatenate((phaseLSSTu,phaseLSSTg,phaseLSSTr,phaseLSSTi,phaseLSSTz,phaseLSSTy))
#    magLSSTAll=np.concatenate((magLSSTu,magLSSTg,magLSSTr,magLSSTi,magLSSTz,magLSSTy))
    
    if len(magLSSTu)!=0: 
        deltau=max(magLSSTu)-min(magLSSTu)
        meanu=np.mean(magModelu)
    else:
        deltau=.1
        meanu=10
        
    if len(magLSSTg)!=0: 
        deltag=max(magLSSTg)-min(magLSSTg)
        meang=np.mean(magModelg)
    else:
        deltag=.1
        meang=10
        
    if len(magLSSTr)!=0: 
        deltar=max(magLSSTr)-min(magLSSTr)
        meanr=np.mean(magModelr)
    else:
        deltar=.1
        meanr=10
        
    if len(magLSSTi)!=0: 
        deltai=max(magLSSTi)-min(magLSSTi)
        meani=np.mean(magModeli)
    else:
        deltai=.1
        meani=10
        
    if len(magLSSTz)!=0: 
        deltaz=max(magLSSTz)-min(magLSSTz)
        meanz=np.mean(magModelz)
    else:
        deltaz=.1
        meanz=10
        
    if len(magLSSTy)!=0: 
        deltay=max(magLSSTy)-min(magLSSTy)
        meany=np.mean(magModely)
    else:
        deltay=.1
        meany=10
    
    deltamax=max([deltau,deltag,deltar,deltai,deltaz,deltay])


#    ax2 = plt.GridSpec(5, 1, top=0.8)
#    gs_base = plt.GridSpec(5, 1, hspace=0)
    
    fig=plt.figure(figsize=(10,16), dpi=80)
    plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95
                        , left = 0.1, hspace = 0.08, wspace = 0.2)
    fig.suptitle(title)
   

    ax2 = plt.subplot2grid((4,2), (0,0),  colspan=2, rowspan=1) # topleft    
    ax2.set_ylabel('ALLbands (mag)')
    ax2.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax2.invert_yaxis()
    ax2.set_xlim([0,1])
    #ax2.plot(phaseGAIA1,magGAIA1,'o',color='red')
    ax2.plot(phaseLSSTu,magLSSTu,'o',color='purple')
    ax2.plot(phaseLSSTg,magLSSTg,'o',color='g')
    ax2.plot(phaseLSSTr,magLSSTr,'o',color='r')
    ax2.plot(phaseLSSTi,magLSSTi,'o',color='k')
    ax2.plot(phaseLSSTz,magLSSTz,'o',color='magenta')
    ax2.plot(phaseLSSTy,magLSSTy,'o',color='y')
    ax2.set_xlabel('phase')

    
    ax3 = plt.subplot2grid((4,2), (1,0),  colspan=1, rowspan=1) # topleft    
    ax3.set_ylabel('u (mag)')
    ax3.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax3.invert_yaxis()
    ax3.set_xlim([0,1])
    ax3.set_ylim([meanu+.7*deltamax,meanu-.7*deltamax])
    ax3.set_xticklabels([])
    #ax3.set_ylabel('G (mag)                    ')
    #ax3.set_xlabel('phase')
    ax3.plot(phaseModelu, magModelu,'.')
    #ax3.plot(phaseGAIA1,magGAIA1,'o',color='red')
    ax3.plot(phaseLSSTu,magLSSTu,'o',color='r')

    ax4 = plt.subplot2grid((4,2), (2,0),  colspan=1, rowspan=1) # topleft    
    ax4.set_ylabel('g (mag)')
    #ax4.set_ylabel('G (mag)')
    ax4.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax4.invert_yaxis()
    ax4.set_xlim([0,1])
    ax4.set_ylim([meang+.7*deltamax,meang-.7*deltamax])
    ax4.set_xticklabels([])
    #ax4.set_xlabel('phase')
    ax4.plot(phaseModelg, magModelg,'.')
    #ax4.plot(phaseGAIA1,magGAIA1,'o',color='red')
    ax4.plot(phaseLSSTg,magLSSTg,'o',color='r')
    
    ax5 = plt.subplot2grid((4,2), (3,0),  colspan=1, rowspan=1) # topleft    
    ax5.set_ylabel('r (mag)')
    ax5.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax5.invert_yaxis()
    ax5.set_xlim([0,1])
    ax5.set_ylim([meanr+.7*deltamax,meanr-.7*deltamax])
    ax5.set_xlabel('phase')
    ax5.plot(phaseModelr, magModelr,'.')
    #ax5.plot(phaseGAIA1,magGAIA1,'o',color='red')
    ax5.plot(phaseLSSTr,magLSSTr,'o',color='r')
    
    ax6 = plt.subplot2grid((4,2), (1,1),  colspan=1, rowspan=1) # topleft    
    ax6.set_ylabel('i (mag)')
    ax6.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax6.invert_yaxis()
    ax6.set_xticklabels([])
    ax6.set_xlim([0,1])
    ax6.set_ylim([meani+.7*deltamax,meani-.7*deltamax])
    #ax6.set_xlabel('phase')
    ax6.plot(phaseModeli, magModeli,'.')
    #ax6.plot(phaseGAIA1,magGAIA1,'o',color='red')
    ax6.plot(phaseLSSTi,magLSSTi,'o',color='r')
    
    ax7 = plt.subplot2grid((4,2), (2,1),  colspan=1, rowspan=1) # topleft    
    ax7.set_ylabel('z (mag)')
    ax7.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax7.invert_yaxis()
    ax7.set_xlim([0,1])
    ax7.set_ylim([meanz+.7*deltamax,meanz-.7*deltamax])
    #ax7.set_xlabel('phase')
    ax7.plot(phaseModelz, magModelz,'.')
    #ax7.plot(phaseGAIA1,magGAIA1,'o',color='red')
    ax7.plot(phaseLSSTz,magLSSTz,'o',color='r')
    ax7.set_xticklabels([])
    
    ax8 = plt.subplot2grid((4,2), (3,1),  colspan=1, rowspan=1) # topleft    
    ax8.set_ylabel('y (mag)')
    ax8.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax8.invert_yaxis()
    ax8.set_xlim([0,1])
    ax8.set_ylim([meany+.7*deltamax,meany-.7*deltamax])
#    ax8.set_xlabel('phase')
    ax8.plot(phaseModely, magModely,'.')
    #ax8.plot(phaseGAIA1,magGAIA1,'o',color='red')
    ax8.plot(phaseLSSTy,magLSSTy,'o',color='r')
    ax8.set_xlabel('phase')


#    plt.savefig('LC_interpolata_noised_allband.pdf')

def plotting_LSST_saturation(title,phaseModelu,magModelu,phaseModelg,magModelg,
              phaseModelr,magModelr,phaseModeli,magModeli,phaseModelz,magModelz,phaseModely,magModely,
              phaseLSSTu,magLSSTu,phaseLSSTg,magLSSTg,phaseLSSTr,magLSSTr,
              phaseLSSTi,magLSSTi,phaseLSSTz,magLSSTz,phaseLSSTy,magLSSTy,
#              phaseLSSTu_sat,magLSSTu_sat,phaseLSSTg_sat,magLSSTg_sat,phaseLSSTr_sat,magLSSTr_sat,
#              phaseLSSTi_sat,magLSSTi_sat,phaseLSSTz_sat,magLSSTz_sat,phaseLSSTy_sat,magLSSTy_sat,
              phaseLSSTu_satlevel,magLSSTu_satlevel,phaseLSSTg_satlevel,magLSSTg_satlevel,phaseLSSTr_satlevel,magLSSTr_satlevel,
              phaseLSSTi_satlevel,magLSSTi_satlevel,phaseLSSTz_satlevel,magLSSTz_satlevel,phaseLSSTy_satlevel,magLSSTy_satlevel,
                  snru,snrg,snrr,snri,snrz,snry):
    
#    phaseLSSTAll=np.concatenate((phaseLSSTu,phaseLSSTg,phaseLSSTr,phaseLSSTi,phaseLSSTz,phaseLSSTy))
#    magLSSTAll=np.concatenate((magLSSTu,magLSSTg,magLSSTr,magLSSTi,magLSSTz,magLSSTy))
    
    if len(magLSSTu)!=0: 
        deltau=max(magLSSTu)-min(magLSSTu)
        meanu=np.mean(magModelu)
    else:
        deltau=.1
        meanu=10
        
    if len(magLSSTg)!=0: 
        deltag=max(magLSSTg)-min(magLSSTg)
        meang=np.mean(magModelg)
    else:
        deltag=.1
        meang=10
        
    if len(magLSSTr)!=0: 
        deltar=max(magLSSTr)-min(magLSSTr)
        meanr=np.mean(magModelr)
    else:
        deltar=.1
        meanr=10
        
    if len(magLSSTi)!=0: 
        deltai=max(magLSSTi)-min(magLSSTi)
        meani=np.mean(magModeli)
    else:
        deltai=.1
        meani=10
        
    if len(magLSSTz)!=0: 
        deltaz=max(magLSSTz)-min(magLSSTz)
        meanz=np.mean(magModelz)
    else:
        deltaz=.1
        meanz=10
        
    if len(magLSSTy)!=0: 
        deltay=max(magLSSTy)-min(magLSSTy)
        meany=np.mean(magModely)
    else:
        deltay=.1
        meany=10
    
    deltamax=max([deltau,deltag,deltar,deltai,deltaz,deltay])


#    ax2 = plt.GridSpec(5, 1, top=0.8)
#    gs_base = plt.GridSpec(5, 1, hspace=0)
    
    fig=plt.figure(figsize=(10,16), dpi=80)
    plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95
                        , left = 0.1, hspace = 0.08, wspace = 0.2)
    fig.suptitle(title)
   

    ax2 = plt.subplot2grid((4,2), (0,0),  colspan=2, rowspan=1) # topleft    
    ax2.set_ylabel('ALLbands (mag)')
    ax2.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax2.invert_yaxis()
    ax2.set_xlim([0,1])
    #ax2.plot(phaseGAIA1,magGAIA1,'o',color='red')
    ax2.plot(phaseLSSTu,magLSSTu,'o',color='purple')
    ax2.plot(phaseLSSTg,magLSSTg,'o',color='g')
    ax2.plot(phaseLSSTr,magLSSTr,'o',color='r')
    ax2.plot(phaseLSSTi,magLSSTi,'o',color='k')
    ax2.plot(phaseLSSTz,magLSSTz,'o',color='magenta')
    ax2.plot(phaseLSSTy,magLSSTy,'o',color='y')
    ax2.set_xlabel('phase')

    
    ax3 = plt.subplot2grid((4,2), (1,0),  colspan=1, rowspan=1) # topleft    
    ax3.set_ylabel('u (mag)')
    ax3.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax3.invert_yaxis()
    ax3.set_xlim([0,1])
    ax3.set_ylim([meanu+.7*deltamax,meanu-.7*deltamax])
    ax3.set_xticklabels([])
    #ax3.set_ylabel('G (mag)                    ')
    #ax3.set_xlabel('phase')
    ax3.plot(phaseModelu, magModelu,'.')
    #ax3.plot(phaseGAIA1,magGAIA1,'o',color='red')
    ax3.plot(phaseLSSTu,magLSSTu,'o',color='r')

    ax4 = plt.subplot2grid((4,2), (2,0),  colspan=1, rowspan=1) # topleft    
    ax4.set_ylabel('g (mag)')
    #ax4.set_ylabel('G (mag)')
    ax4.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax4.invert_yaxis()
    ax4.set_xlim([0,1])
    ax4.set_ylim([meang+.7*deltamax,meang-.7*deltamax])
    ax4.set_xticklabels([])
    #ax4.set_xlabel('phase')
    ax4.plot(phaseModelg, magModelg,'.')
    #ax4.plot(phaseGAIA1,magGAIA1,'o',color='red')
    ax4.plot(phaseLSSTg,magLSSTg,'o',color='r')
    
    ax5 = plt.subplot2grid((4,2), (3,0),  colspan=1, rowspan=1) # topleft    
    ax5.set_ylabel('r (mag)')
    ax5.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax5.invert_yaxis()
    ax5.set_xlim([0,1])
    ax5.set_ylim([meanr+.7*deltamax,meanr-.7*deltamax])
    ax5.set_xlabel('phase')
    ax5.plot(phaseModelr, magModelr,'.')
    #ax5.plot(phaseGAIA1,magGAIA1,'o',color='red')
    ax5.plot(phaseLSSTr,magLSSTr,'o',color='r')
    
    ax6 = plt.subplot2grid((4,2), (1,1),  colspan=1, rowspan=1) # topleft    
    ax6.set_ylabel('i (mag)')
    ax6.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax6.invert_yaxis()
    ax6.set_xticklabels([])
    ax6.set_xlim([0,1])
    ax6.set_ylim([meani+.7*deltamax,meani-.7*deltamax])
    #ax6.set_xlabel('phase')
    ax6.plot(phaseModeli, magModeli,'.')
    #ax6.plot(phaseGAIA1,magGAIA1,'o',color='red')
    ax6.plot(phaseLSSTi,magLSSTi,'o',color='r')
    
    ax7 = plt.subplot2grid((4,2), (2,1),  colspan=1, rowspan=1) # topleft    
    ax7.set_ylabel('z (mag)')
    ax7.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax7.invert_yaxis()
    ax7.set_xlim([0,1])
    ax7.set_ylim([meanz+.7*deltamax,meanz-.7*deltamax])
    #ax7.set_xlabel('phase')
    ax7.plot(phaseModelz, magModelz,'.')
    #ax7.plot(phaseGAIA1,magGAIA1,'o',color='red')
    ax7.plot(phaseLSSTz,magLSSTz,'o',color='r')
    ax7.set_xticklabels([])
    
    ax8 = plt.subplot2grid((4,2), (3,1),  colspan=1, rowspan=1) # topleft    
    ax8.set_ylabel('y (mag)')
    ax8.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in',which='major')
    ax8.invert_yaxis()
    ax8.set_xlim([0,1])
    ax8.set_ylim([meany+.7*deltamax,meany-.7*deltamax])
#    ax8.set_xlabel('phase')
    ax8.plot(phaseModely, magModely,'.')
    #ax8.plot(phaseGAIA1,magGAIA1,'o',color='red')
    ax8.plot(phaseLSSTy,magLSSTy,'o',color='r')
    ax8.set_xlabel('phase')


#    plt.savefig('LC_interpolata_noised_allband_saturation.pdf')

