#!/usr/bin/env python
# coding: utf-8

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

# ## History
# this code have been translate from the matlab code of julien brondex used in Brondex, 2019 by B. Urruty and cleaned by P. Mathiot (2022).
# 
# ## Usage
# It take the data from the inversion to compute the parameter of different friction

# ### Modules

# In[1]:


# load module
import numpy as np
import xarray as xr


# ## Files names

# In[2]:


cfilein='ANT50.GL1-ISMIP6_5.restart_weertman_only.nc'
cfileout='ANT50.GL1-ISMIP6_5.restart_beta_coulomg_reg.nc'


# ### Constants

# In[3]:


# Ice sheet parameters
m = 1.0 / 3.0  # Exposant loi Weertman NL ou Schoof
Cmax = 0.4     # Cmax loi Schoof
hth=75         # [m]   treshold on high above flotation 
u0=300         # [m/a]  

# Physical constant
z_sl=0         # [m] sea level heigh
yrs=31556926.0 # [s]
g = 9.81 * yrs**2
rho_ice = 917.0  * 1.0e-6 /(yrs**2)
rho_sea = 1027.0 * 1.0e-6 /(yrs**2)


# ### Load input files

# In[4]:


# load data
ds_weert=xr.open_dataset(cfilein)
ds_weert

ds_cr=ds_weert.copy(deep=True)

# In[5]:


# geometry
da_haf   = ds_weert['haf']
da_gmask = ds_weert['groundedmask'] 

# velocities
da_u     = ds_weert['ssavelocity 1']
da_v     = ds_weert['ssavelocity 2']

# beta
da_wbeta = ds_weert['beta']


# ### Retreive taub in Weertmann linear case

# In[6]:


# Weertman linear
# Retreive friction parameter
C = 10**da_wbeta

# Compute |ub|
u_norm=np.sqrt(da_u**2 + da_v**2)

# Compute Weertmann taub
taub = C * u_norm


# ### Convert beta weertmann to beta Coulomb

# In[7]:


coef=(u_norm/(u_norm+u0))

Beta_j=taub/(coef**m)

# ### Restore Beta
# In order to have a generic file, we need to increase beta when haf below a treshold (hth). The correction is based on haf.

# In[9]:


# coefficient reducteur lambda
# haf > hth => no change
# haf < 0   => floating part 
lbd=da_haf/hth
lbd.values[lbd.values>1]=1
lbd.values[lbd.values<0]=1

# définition du Beta qui sera injecté dans le modèle
Beta_out1=Beta_j*1.0/lbd

# define a background value for floating part at 100kPa
Beta_out1.values[da_gmask.values == -1] = 0.1

# add removed attributes
Beta_out1.attrs=da_wbeta.attrs

# ## Save data

# In[10]:

ds_cr['beta_cr']=Beta_out1
ds_cr.attrs['History']=('Convertion of the Weertmann beta parameter'
                              'from file '+cfilein+' to Coulomb regularised'
                              'beta parameter using Beta_Weertmann2Joughin_nc script')
ds_cr.to_netcdf(cfileout)

