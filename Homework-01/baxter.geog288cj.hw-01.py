#!/usr/bin/env python
# coding: utf-8

# # __GEOG 288CJ Homework-01__
# ## Ian Baxter

# In[1]:


# import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline
import metpy.calc as mpcalc
from metpy.units import units


# ## 1) Review the definition of adiabatic lapse rate and static stability. A good reference is Wallace and Hobbs Atmospheric Science: an introductory survey

# Adiabatic lapse rate: The rate of temperature change for a parcel of air that is being raised or lowered in the atmosphere without losing or gaining heat.
# 
# Static stability: Refers to the restoring force of the atmosphere, in which vertical mixing is inhibited by greater stability. Static stability is commonly computed as change in potential temperature over change in height.

# ## 2C)

# <img src="./question2.png">

# The diameter of the major eddies are (based on Taylor's hypothesis)...
# 
# P = lambda / M 
# 
# lambda = P * M
# 
# lambda = 60 s * 5 m/s = 300 m

# ## 3) Would there be a boundary layer on a planet that had an atmosphere. but that did not experience a diurnal variation of net radiation at the ground?

# Yes, there would be a boundary layer on a planet that has an atmosphere. There will be a boundary layer as long as there is a an atmosphere because there will still be a distribution of momentum and heat fluxes associated with the surface. As long as there is heat coming off of the surface or winds passing over the surface (and being affected by friction), there will be a boundary layer. Even at night on Earth there exists a boundary layer.

# ## 4) There are four radiosonde profiles in Gaucho Space named A, B, C and D. In each file, columns have time (s), Pressure level (hPa), Air Temperature (C), Relative Humidity (%), Wind Speed (knots), Wind Direction (degrees), Altitude (m) and Geopotential Height (m)

# ## 4a) For each sounding, compute the water vapor mixing ratio, potential temperature and virtual potential temperature. Use surface pressure the value at t=0

# In[2]:


# Calculations

# Read in radiosonde data
def read_radiosonde(fname):
    names = ['Time','P','T','U','Wsp','Wdir','Altitude','Geo_Pot']
    return pd.read_csv(fname,skiprows=[0],names=names,delimiter='\t')

# Compute saturation vapor pressure
def calc_saturation_vapor_pressure(T):
    return 0.6112 * np.exp((17.67 * T) / (T + 243.5))

def calc_saturation_mixing_ratio(es,P):
    es = es * 10
    return 0.622 * es / (P - es)

# Compute water vapor mixing ratio
def calc_water_vapor_mixing_ratio(U,rs):
    return (U / 100) * rs

# Compute potential temperature (dry)
def calc_potential_temperature(T,Po,P):
    return (T + 273.15) * (Po / P) ** (2/7)

# Compute virtual potential temperature
def calc_virtual_potential_temperature(theta,r,r_l):
    #return theta * (1 + 0.61 * r - r_l)
    return theta * (r + 0.622) / (0.622 * (1 + r))

def calc(df):
    es = calc_saturation_vapor_pressure(df['T'])
    rs = calc_saturation_mixing_ratio(es,df['P'])
    r = calc_water_vapor_mixing_ratio(df['U'],rs)
    theta = calc_potential_temperature(df['T'],df['P'][0],df['P'])
    theta_v = calc_virtual_potential_temperature(theta,r,0)
    return r,theta,theta_v


# In[3]:


#---Plot
fig = plt.figure(figsize=(20,8))
left = -0.2
bottom = 1.1

fnum = ['A','B','C','D']
labs_list = {'A':'a','B':'b','C':'c','D':'d'}

for idx,f in enumerate(fnum):
    df = read_radiosonde('Radiosonde-'+f+'.txt')
    r,theta,theta_v = calc(df)
    
    #---Plots
    ax = fig.add_subplot(1,4,idx+1)
    #ax.set_title('Water Vapor Mixing Ratio',loc='left')
    ax.set_xlim([0,10])
    ax.set_ylim([0,3])
    ax.set_xlabel('Mixing ratio [g/kg]',c='dodgerblue')
    ax.grid()
    if idx == 0:
        ax.set_ylabel('Altitude [km]')

    #---Labels
    ax.text(left, bottom, labs_list[f],
             horizontalalignment='left',
             verticalalignment='bottom',
             weight='bold',
             fontsize='12',
             transform=ax.transAxes)

    ax.text(left+0.1, bottom, 'Radiosonde-'+str(f),
             horizontalalignment='left',
             verticalalignment='bottom',
             fontsize='12',
             transform=ax.transAxes)
    
    #---Mixing ratio
    ln1 = ax.plot(r*1000,df['Altitude']/1000, label='Mixing Ratio [g/kg]',
                 c='dodgerblue')
    
    #---Wind speed
    ax2 = ax.twiny()
    ax2.xaxis.set_ticks_position('bottom')
    rspine = ax2.spines['bottom']
    rspine.set_position(('axes', -0.08))
    ax2.set_xlabel('Wind speed [m/s]',c='fuchsia')
    ax2.xaxis.set_label_position('bottom')
    ax2.set_xlim([0,30])
    
    ln2 = ax2.plot(df['Wsp'],df['Altitude']/1000, label='Wind Speed [m/s]',
                 c='fuchsia')
    
    #---Potential Temperature & Virtual Potential Temperature
    ax3 = ax.twiny()
    ax3.set_xlim([285,310])
    ax3.set_ylim([0,3])
    ax3.set_xlabel('Potential Temperature [K]')
    
    ln3 = ax3.plot(theta,df['Altitude']/1000, label='Potential Temperature [K]',
                  c='red')
    ln4 = ax3.plot(theta_v,df['Altitude']/1000, label='Virtual Potential Temperature [K]',
                  c='green', linestyle='--')

    #---Layer lines
    #---Radiosonde-A Layers
    if idx == 0:
        ax.axhline(0.1,c='black',linestyle='--',linewidth=1)
        ax.text(1.04, 0.005, 'SL',
                ha='left',
                va='bottom',
                fontsize='12',
                transform=ax.transAxes)
        ax.axhline(2.25,c='black',linestyle='--',linewidth=1)
        ax.text(1.04, 0.47, 'ML',
                ha='left',
                va='bottom',
                fontsize='12',
                transform=ax.transAxes)
        ax.text(1.04, 0.83, 'FA',
             ha='left',
             va='bottom',
             fontsize='12',
             transform=ax.transAxes)
    #---Radiosonde-B Layers
    elif idx == 1:
        ax.axhline(0.1,c='black',linestyle='--',linewidth=1)
        ax.text(1.04, 0.0, 'SBL',
                ha='left',
                va='bottom',
                fontsize='12',
                transform=ax.transAxes)
        ax.axhline(0.65,c='black',linestyle='--',linewidth=1)
        ax.text(1.04, 0.1, 'RL',
                ha='left',
                va='bottom',
                fontsize='12',
                transform=ax.transAxes)
        ax.text(1.04, 0.5, 'FA',
                ha='left',
                va='bottom',
                fontsize='12',
                transform=ax.transAxes)
    #---Radiosonde-C Layers
    elif idx == 2:
        ax.axhline(0.18,c='black',linestyle='--',linewidth=1)
        ax.text(1.04, 0.01, 'SBL',
                ha='left',
                va='bottom',
                fontsize='12',
                transform=ax.transAxes)
        ax.axhline(0.28,c='black',linestyle='--',linewidth=1)
        ax.text(1.04, 0.06, 'RL',
                ha='left',
                va='bottom',
                fontsize='12',
                transform=ax.transAxes)
        ax.text(1.04, 0.5, 'FA',
                ha='left',
                va='bottom',
                fontsize='12',
                transform=ax.transAxes)
    #---Radiosonde-D Layers
    elif idx == 3:
        ax.axhline(0.365,c='black',linestyle='--',linewidth=1)
        ax.axhline(1.20,c='black',linestyle='--',linewidth=1)
        ax.text(1.04, 0.05, 'SL',
             ha='left',
             va='bottom',
             fontsize='12',
             transform=ax.transAxes)
    #---Legend
    if idx == 3:
        lns = ln1+ln2+ln3+ln4
        labs = [l.get_label() for l in lns]
        ax.legend(lns, labs, bbox_to_anchor=(1.1, 0.5), frameon=True)

#---Figure spacing & Title
fig.subplots_adjust(wspace=1.0)
fig.suptitle('4b) ', x=0.05, y=0.95,
             ha='left',
             va='bottom',
             fontsize='24');


# ## 4b) For each sounding, make vertical profile plots of mixing ratio, temperature, virtual potential temperature and wind speeds. Limit the plots to elevations 0-3000 m so that variations can be better seen.

# ## 4c) Discuss each profile and identify (if possible) surface layer (SL), mixed layer (ML), stable boundary layer (SBL) and free atmosphere (FA)

# Radionsonde-A has a very well mixed layer, where potential temperature and virtual potential temperature barely changes vertically for 1.5 km.
# 
# Radiosonde-B has a very stable atmospheric structure. 
# 
# Radiosonde-C 
# 
# Radiosonde-D is very unstable near the surface, to the point of being super adiabatic, beneath a well mixed layer. The rapid decrease in potential temperature with height, is indicative of this instability and suggests that this may be an afternoon profile. It may be that solar radiation over the duration of the day and warms the surface, the near surface air becomes more buoyant and creates turbulence.

# ## 4d) Estimate the top of the BL in each sounding

# Radiosonde-A:
# 
# Radiosonde-B: 625 m
# 
# Radiosonde-C: 

# ## 4e) Associate each sounding with early morning, noon, afternoon, late afternoon and night. Discuss and justify your answers.
# 
# Radiosonde-A: Afternoon
# 
# Radiosonde-B: Early Morning
# 
# Radiosonde-C: Noon
# 
# Radiosonde-D: Late Afternoon

# In[ ]:




