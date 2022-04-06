import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
%matplotlib inline
from scipy.stats import norm
import time

L = 4
M = 5
T = 250
X0 = 100
K = 80
r = 0.06
sigma = 0.4
N0 = 10**4

sum = 0

for l in range(1,L):
    Nl = int(N * (pow(pow(M,l-1),-0.5) + pow(pow(M,l),-0.5))/np.sqrt(pow(M,l-1)+pow(M,l)))
    if l == 1:
        X    = np.zeros((T,Nl), np.float64)
        X[0] = X0
        for t in range(1,T):
            rand  = np.random.standard_normal(Nl)
            X[t] = X[t-1]*(np.exp((r-0.5*sigma**2)*1/T + sigma*np.sqrt(1/T)*rand))
    else:
        Nl    = M**l
        Ml    = int(pow(M,L+1)/Nl) * Nl
        Xf    = np.zeros((int(T*M/Nl), Ml), np.float64)
        Xf[0] = X0
        for t in range(1,int(T*M/Nl)):
            rand  = np.random.standard_normal(Ml)
            Xf[t] = Xf[t-1]*np.exp((r-0.5*sigma**2)*1/(Nl*T) + sigma*np.sqrt(1/(Nl*T))*rand)
        Xc    = np.zeros((int(T/Nl), Ml), np.float64)
        Xc[0] = X0
        for t in range(1,int(T/Nl)):
            Xc[t] = Xc[t-1]*np.exp((r-0.5*sigma**2)*M/(Nl*T) + sigma*np.sqrt(M/(Nl*T))*(rand[t] + rand[t+M]))
        sum += np.exp(-r*T)*(np.maximum(Xf[-1,l]-K,0) - np.maximum(Xc[-1,l]-K,0))/Nl

# plt.figure(figsize=(12,5))
# plt.plot(Xf,linewidth=1)