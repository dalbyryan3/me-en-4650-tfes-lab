# %% 
# ME EN 4650  Lab10:Heat Exchanger Lab  Ryan Dalby    
import numpy as np
from numpy import random
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from pandas.plotting import table
import os
import matlab.engine # Must install from matlab installation direcotry using: cd "matlabroot\extern\engines\python" && python setup.py install 

# %%
# Define useful functions
mmhg_to_pa = lambda mmhg : mmhg*133.3 # Lambda function to covert mmHg to Pa
inh2o_to_pa = lambda inh2o : inh2o*249.04 # Lambda function to covert inH20 to Pa
in_to_m = lambda inch : inch*0.0254 # Lambda function to convert in to m
cm_to_m = lambda cm : cm/100 # Lambda function to convert cm to m
mm_to_m = lambda mm : mm/1000 # Lambda function to convert mm to m
degCToK = lambda degC : degC+273.15 # Lambda function to covert degC to K

# Generate EffectivenessNTU figure
# Plot effectiveness versus NTU and Cr for crossflow heat exchanger

Cr=np.array([0,.25,.5,.75])
NTU=np.linspace(0,1,100)
e=np.zeros((100,5))

plt.figure(figsize=(10,5))
plt.minorticks_on()
plt.grid(b=True, which='major', axis='both', linestyle='-')
plt.grid(b=True, which='minor', axis='both', linestyle='--')
plt.xlim((0,1))
plt.ylim((0,0.7))
plt.title('Effectiveness versus NTU for a crossflow heat exchanger', fontsize=16)
plt.xlabel('NTU', fontsize=16)
plt.ylabel('$\epsilon$', fontsize=16)
for k in range(Cr.size):
    e[:,k]=(1-np.exp(-NTU*(1-Cr[k])))/(1-Cr[k]*np.exp(-NTU*(1-Cr[k])))
    plt.plot(NTU,e[:,k],'k-')
    plt.text(NTU[-1]*1.02,e[-1,k]*.98,str(Cr[k]), fontsize=16)
k=k+1
e[:,k]=NTU/(1+NTU)
plt.plot(NTU,e[:,k],'k-')
plt.text(NTU[-1]*1.02,e[-1,k]*.98,'1', fontsize=16)
plt.text(NTU[-1]*1.02,e[-1,1]*1.11,'$C_r$', fontsize=16)
plt.show()