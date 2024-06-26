import os
import glob
import arviz as az
import numpy as np
import pandas as pd
from pytensor.tensor import TensorVariable
from typing import Optional, Tuple
import pytensor.tensor as pt

import pymc as pm
import jax
import seaborn as sns
# import aesara
# import aesara.tensor as at
# import pytensor.tensor as pt
import matplotlib.pyplot as plt
import xarray as xr
import re
from scipy import stats
from scipy.special import logit
from scipy.special import expit as logist
from matplotlib.axes import Axes
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from typing import List
from pandas import DataFrame
import sys
import pickle

jax.default_backend()
jax.devices()

RANDOM_SEED = 8927
rng = np.random.default_rng(RANDOM_SEED)

path_out = './pymc_model'
with open(os.path.join('pymc_model','brain_16mm_vshd.pkl'),'rb') as f:
    inputBrain16um = pickle.load(f)
df = inputBrain16um['dsg']
inf_aver = inputBrain16um['wsf']

g_gf = inf_aver.to_numpy()
d_sg = df.to_numpy()
ngens = len(inf_aver)
nspost = len(df)
nfeat = len(inf_aver.columns)

with pm.Model() as simNB_model:
    mu_m = pm.Gamma("mu_m",1,1)
    o_m = pm.Exponential("o_m",3)
    a_m = pm.Deterministic("a_m",1/o_m**2)
    m_g = pm.Gamma("m_g",a_m,a_m/mu_m,shape = (1,ngens))
    
    a_y = pm.Gamma("a_y",10,10/20)
    y_s = pm.Gamma("y_s",10,10/30,shape = (nspost,1))
    
    b0 = pm.Gamma("b0",9,3,shape = 2)
    o_g = pm.Exponential("o_g",b0[0], shape = (nspost,1))
    a_g = pm.Deterministic("a_g",1/o_g**2)
    
    o_es = pm.Exponential("o_es",b0[1])
    a_es = pm.Deterministic("a_es",1/o_es**2)
    mu_es = pm.Gamma("mu_es",1,100)
    s_eg = pm.Gamma("s_eg",a_es,a_es/mu_es,shape = (1,ngens))
    
    
    k_r  = pm.Gamma("k_r",7,1)
    x_rf  = pm.Gamma("x_rf",k_r/50,k_r,shape = (50,nfeat))
    B_s = pm.Gamma("B_s",7,1) 
    N_s =pm.Gamma("N_s",4*5,5) #4 cells per spot vn =5 acording to seggestion of section 2.1  
    z_sr  = pm.Gamma("z_sr",B_s/50,1/(N_s/B_s),shape = (nspost,50))
    mu_sf = pm.Deterministic("mu_sf",pm.math.dot(z_sr,x_rf))
    w_sf = pm.Gamma("w_sf",alpha = mu_sf.eval()*10, beta =10)
    
    musg1_p1 = pm.math.dot(w_sf,g_gf.T)
    musg1_p2 = m_g * musg1_p1
    musg1_p3 = musg1_p2+s_eg
    musg1_p4 = pm.Deterministic("musg1_p4",musg1_p3*y_s)
    
    count = pm.NegativeBinomial("count", mu=musg1_p4, alpha=a_g, observed=d_sg)
    trace = pm.sample(draws=100,tune = 10, target_accept=0.95, random_seed=RANDOM_SEED, chains = 10)

with open(os.path.join(path_out, 'trace_visHD.pkl'), 'wb') as buff:
    pickle.dump({'trace':trace})
