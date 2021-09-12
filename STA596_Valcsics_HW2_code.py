# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 20:31:38 2021

@author: angya
"""
import  matplotlib.pyplot  as  pyplot
from pydataset import data
import numpy as np

dat = data('chickwts')
data('chickwts', show_doc=(True))

wt = dat['weight']
feed = dat['feed']
np.unique(feed)

# A bootstrap CI for single mean
tot = len(wt)
xbar = np.mean(wt)

B = 10000
xbar_boot = np.zeros(B)

for k in range(B):
    xbar_boot[k] = np.mean(np.random.choice(wt,tot, replace = True))

alp = 0.05
[lo, hi] = np.quantile(xbar_boot, [alp/2, 1 - alp/2])

np.round([lo, hi], 3) # array([243.549, 279.099])
pyplot.hist(xbar_boot)
pyplot.axvline(xbar,color='black')
pyplot.axvline(lo, color = 'blue')
pyplot.axvline(hi, color = 'blue')

# Hypothesis test for single mean
mu0 = 250
# shifted data with mean mu0
x_star = wt - (xbar - mu0)
B = 10000
xbar_boot = np.zeros(B)
for k in range(B):
    xbar_boot[k] = np.mean(np.random.choice(x_star, tot, replace = True))
# xbar_boot is now centered at mu0 = 250
# calculate statistics
Swt2 = np.var(wt)**2 # var of x_star is equal to var of wt
T = (xbar - mu0)/np.sqrt(Swt2 / tot)
T_star = (xbar_boot - mu0)/np.sqrt(Swt2 / tot)

alp = 0.05
Tcrit = np.quantile(np.abs(T_star), 1 - alp)
np.abs(T) > Tcrit
# We fail to reject the null hypothesis

# graphically
pyplot.hist(T_star)
pyplot.axvline(T,color='black')
pyplot.axvline(Tcrit,color='red')

# A bootstrap CI for difference of two means
meatmeal_wt = dat['weight'][dat['feed']=='meatmeal']
soybean_wt = dat['weight'][dat['feed']=='soybean']
n = len(meatmeal_wt)
m = len(soybean_wt)
xbar_meat = np.mean(meatmeal_wt)
xbar_soy = np.mean(soybean_wt)
d = xbar_meat - xbar_soy

B = 10000
xbarm_boot = np.zeros(B)
xbars_boot = np.zeros(B)
for k in range(B):
    xbarm_boot[k] = np.mean(np.random.choice(meatmeal_wt, n, replace = True))
    xbars_boot[k] = np.mean(np.random.choice(soybean_wt, m, replace = True))
    
d_boot = xbarm_boot - xbars_boot

alp = 0.05
[lo, hi] = np.quantile(d_boot, [alp/2, 1 - alp/2])
np.round([lo, hi], 3)

pyplot.hist(d_boot)
pyplot.axvline(d,color='black')
pyplot.axvline(lo, color = 'blue')
pyplot.axvline(hi, color = 'blue')

# Hypothesis test for difference of two means
B = 10000
# initialize the list for the test statistic replicate
d_boot = np.zeros(B)
# iterate for the specified number of iterations
for k in range(B):
    # concatenate the two samples into one
    samples = np.append(meatmeal_wt, soybean_wt)
    # permute the entire appended set
    samples_perm = np.random.permutation(samples)
    # create the hypothesized samples
    sample_1_hyp = samples_perm[:n]
    #  and the rest is the second sample
    sample_2_hyp = samples_perm[n:]
    # compute the test statistic replicate and append it to the list of permutation replicates
    d_boot[k] = np.mean(sample_1_hyp)-np.mean(sample_2_hyp)

pyplot.hist(d_boot)
pyplot.axvline(d,color='black')

S1_sqrd = np.var(meatmeal_wt)**2
S2_sqrd = np.var(soybean_wt)**2
test_stat = d/np.sqrt(S1_sqrd/n + S2_sqrd/m)
test_stat_boot = d_boot/np.sqrt(S1_sqrd/n + S2_sqrd/m)
pyplot.hist(test_stat_boot)
pyplot.axvline(test_stat ,color='black')

alp = 0.05
Tcrit = np.quantile(np.abs(test_stat_boot), 1 - alp)
np.abs(test_stat) > Tcrit

pyplot.hist(test_stat_boot)
pyplot.axvline(test_stat,color='black')
pyplot.axvline(Tcrit,color='red')


