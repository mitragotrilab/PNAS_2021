#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt


# In[ ]:


time = np.array(range(24*10+1))


# In[ ]:


# define parameters for SEAPORT
D = np.array([6e7, 6e7, 6e7, 6e7, 1.5e7])   # 4 priming doses, subsequent doses (mg APO)
V = 5000   # volume of distribution (total blood volume in human: 5L)
p = 0.9746     # parameters determined from fitting
f = 0.5742
s = 0.001392
e = 10.51

# PK model equation for SEAPORT
# seaport = (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*time) + f*(1-p)*(e-s)*np.exp(-f*time) + s*p*(e-f)*np.exp(-s*time))


# In[ ]:


### simulation over 10 days
d = int(10)     # 10 days
n_prime = int(1) # initial priming dosing everyday for the first 2 days
freq_prime = 4 # priming dose frequency
n = int(2)      # dosing every two days
freq = freq_prime + int((d-freq_prime*n_prime)/n)   # total dosing frequency

# create an empty array for [APO] 
seaport = np.empty(0, dtype=int)

for i in np.array(range(freq)):  # dosing number index
    if i == 0:   # if priming dose
        D_s = D[i]    # inject 60 mg
        seaport = np.append(seaport, (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*time) + f*(1-p)*(e-s)*np.exp(-f*time) + s*p*(e-f)*np.exp(-s*time)))
        # update seaport array with PK from the first dose
          # add to PK from the first dose at the corresponding timepoint

    if i != 0 and i < freq_prime:     # if priming dose
        D_s = D[i]      # inject 60 mg
        for j in np.array(range(24*(d-n_prime*i)+1)):    # dosing time
            dose = (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*j) + f*(1-p)*(e-s)*np.exp(-f*j) + s*p*(e-f)*np.exp(-s*j))
            # [APO] at each hr from the second dose injection solely from the second dose
            seaport[j+24*n_prime*i] = seaport[j+24*n_prime*i] + dose   
            # add to PK from the first dose at the corresponding timepoint

    if i == 4:     # dose beyond priming dose
        D_s = D[len(D)-1]      # repeat injection of 15 mg
        for j in np.array(range(24*(d-n_prime*freq_prime-n*(i-4))+1)):    # dosing time
            dose = (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*j) + f*(1-p)*(e-s)*np.exp(-f*j) + s*p*(e-f)*np.exp(-s*j))
            # [APO] at each hr from the injection solely from the current dose            
            seaport[j+24*(n_prime*freq_prime+n*(i-4))] = seaport[j+24*(n_prime*freq_prime+n*(i-4))] + dose
            # add to PK from the first dose at the corresponding timepoint

    if i == 5:     # dose beyond priming dose
        D_s = D[len(D)-1]      # repeat injection of 15 mg
        for j in np.array(range(24*(d-n_prime*freq_prime-n*(i-4))+1)):    # dosing time
            dose = (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*j) + f*(1-p)*(e-s)*np.exp(-f*j) + s*p*(e-f)*np.exp(-s*j))
            # [APO] at each hr from the injection solely from the current dose            
            seaport[j+24*(n_prime*freq_prime+n*(i-4))] = seaport[j+24*(n_prime*freq_prime+n*(i-4))] + dose
            # add to PK from the first dose at the corresponding timepoint

    if i == 6:     # dose beyond priming dose
        D_s = D[len(D)-1]      # repeat injection of 15 mg
        for j in np.array(range(24*(d-n_prime*freq_prime-n*(i-4))+1)):    # dosing time
            dose = (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*j) + f*(1-p)*(e-s)*np.exp(-f*j) + s*p*(e-f)*np.exp(-s*j))
            # [APO] at each hr from the injection solely from the current dose            
            seaport[j+24*(n_prime*freq_prime+n*(i-4))] = seaport[j+24*(n_prime*freq_prime+n*(i-4))] + dose
            # add to PK from the first dose at the corresponding timepoint


# In[ ]:


# maximum tolerated dose and minimum effective concentration
MTD = np.repeat(10,len(time))    # MTD: 10 ng/mL
MEC = np.repeat(4,len(time))     # MEC: 4 ng/mL

# plot
plt.figure(figsize=(6, 3))
plt.subplot(111).spines["top"].set_visible(False)
plt.subplot(111).spines["right"].set_visible(False)
plt.subplot(111).get_xaxis().tick_bottom() 
plt.subplot(111).get_yaxis().tick_left() 
plt.yticks(fontsize=14, fontweight='bold', fontname = "Arial")  
plt.xticks(np.arange(0, 24*(d+1), 24), fontsize=14, fontweight='bold', fontname = "Arial")
plt.subplot(111).set_xlabel("Time (hr)", fontsize=14, fontweight='bold', fontname = "Arial")
plt.subplot(111).set_ylabel("Plasma [APO] (ng/mL)", fontsize=14, fontweight='bold', fontname = "Arial")

plt.plot(time, seaport, lw=3, color="red")
plt.plot(time, MEC, lw=3, color='darkslategrey', linestyle='dashed')
plt.plot(time, MTD, lw=3, color='darkslategrey', linestyle='dashed')


# In[ ]:


# define parameters for rApokyn
D_a = 1e7   # mg APO dose
V = 5000    # volume of distribution (total blood volume in human: 5L)
f_a = 0.5318    # parameters determined from fitting
e_a = 14.06

# PK model equation for rApokyn
# rapokyn = (D_a/(V*(e_a-f_a)*e_a))*((-f_a*e_a)*np.exp(-e_a*time) + f_a*e_a*np.exp(-f_a*time))


# In[ ]:


### simulation over 10 days
d = int(10)     # 10 days
n = int(3)      # dosing three times a day (30 times per 10 days)
t = int(24/n)   # dosing interval (every 8 hours)
freq = int(d*n)    # total dosing frequency

# create an empty array for [APO] 
rapokyn = np.empty(0, dtype=int)

for i in np.array(range(freq)):  # dosing number index
    if i == 0:   # if first dose
        rapokyn = np.append(rapokyn, ((D_a/(V*(e_a-f_a)*e_a))*((-f_a*e_a)*np.exp(-e_a*time) + f_a*e_a*np.exp(-f_a*time))))
        # update seaport array with PK from the first dose

    else:     # if second dose and beyond
        i = int(i)
        for j in np.array(range(24*d-t*i+1)):    # dosing time
            dose = (D_a/(V*(e_a-f_a)*e_a))*((-f_a*e_a)*np.exp(-e_a*j) + f_a*e_a*np.exp(-f_a*j))
            # [APO] at each hr from the second dose injection solely from the second dose
            rapokyn[j+t*i] = rapokyn[j+t*i] + dose   
            # add to PK from the first dose at the corresponding timepoint


# In[ ]:


# maximum tolerated dose and minimum effective concentration
MTD = np.repeat(10,len(time))
MEC = np.repeat(4,len(time))

# plot
plt.figure(figsize=(9, 3))
plt.subplot(111).spines["top"].set_visible(False)
plt.subplot(111).spines["right"].set_visible(False)
plt.subplot(111).get_xaxis().tick_bottom() 
plt.subplot(111).get_yaxis().tick_left() 
plt.yticks(fontsize=14, fontweight='bold', fontname = "Arial")  
plt.xticks(np.arange(0, 24*(d+1), 24), fontsize=14, fontweight='bold', fontname = "Arial")
plt.subplot(111).set_xlabel("Time (hr)", fontsize=14, fontweight='bold', fontname = "Arial")
plt.subplot(111).set_ylabel("Plasma [APO] (ng/mL)", fontsize=14, fontweight='bold', fontname = "Arial")

plt.plot(time, rapokyn, lw=3, color="black")
plt.plot(time, MEC, lw=3, color='darkorange', linestyle='dashed')
plt.plot(time, MTD, lw=3, color='darkorange', linestyle='dashed')


# In[ ]:


### SEAPORT simulation over 30 days
time = np.array(range(24*30+1))

d = int(30)     # 30 days
n_prime = int(1) # initial priming dosing everyday for the first 2 days
freq_prime = 4 # priming dose frequency
n = int(2)      # dosing every two days
freq = freq_prime + int((d-freq_prime*n_prime)/n)   # total dosing frequency

# create an empty array for [APO] 
seaport = np.empty(0, dtype=int)

for i in np.array(range(freq)):  # dosing number index
    if i == 0:   # if priming dose
        D_s = D[i]    # inject 60 mg
        seaport = np.append(seaport, (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*time) + f*(1-p)*(e-s)*np.exp(-f*time) + s*p*(e-f)*np.exp(-s*time)))
        # update seaport array with PK from the first dose
          # add to PK from the first dose at the corresponding timepoint

    if i != 0 and i < freq_prime:     # if priming dose
        D_s = D[i]      # inject 60 mg
        for j in np.array(range(24*(d-n_prime*i)+1)):    # dosing time
            dose = (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*j) + f*(1-p)*(e-s)*np.exp(-f*j) + s*p*(e-f)*np.exp(-s*j))
            # [APO] at each hr from the second dose injection solely from the second dose
            seaport[j+24*n_prime*i] = seaport[j+24*n_prime*i] + dose   
            # add to PK from the first dose at the corresponding timepoint

    if i == 4:     # dose beyond priming dose
        D_s = D[len(D)-1]      # repeat injection of 15 mg
        for j in np.array(range(24*(d-n_prime*freq_prime-n*(i-4))+1)):    # dosing time
            dose = (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*j) + f*(1-p)*(e-s)*np.exp(-f*j) + s*p*(e-f)*np.exp(-s*j))
            # [APO] at each hr from the injection solely from the current dose            
            seaport[j+24*(n_prime*freq_prime+n*(i-4))] = seaport[j+24*(n_prime*freq_prime+n*(i-4))] + dose
            # add to PK from the first dose at the corresponding timepoint

    if i == 5:     # dose beyond priming dose
        D_s = D[len(D)-1]      # repeat injection of 15 mg
        for j in np.array(range(24*(d-n_prime*freq_prime-n*(i-4))+1)):    # dosing time
            dose = (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*j) + f*(1-p)*(e-s)*np.exp(-f*j) + s*p*(e-f)*np.exp(-s*j))
            # [APO] at each hr from the injection solely from the current dose            
            seaport[j+24*(n_prime*freq_prime+n*(i-4))] = seaport[j+24*(n_prime*freq_prime+n*(i-4))] + dose
            # add to PK from the first dose at the corresponding timepoint

    if i == 6:     # dose beyond priming dose
        D_s = D[len(D)-1]      # repeat injection of 15 mg
        for j in np.array(range(24*(d-n_prime*freq_prime-n*(i-4))+1)):    # dosing time
            dose = (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*j) + f*(1-p)*(e-s)*np.exp(-f*j) + s*p*(e-f)*np.exp(-s*j))
            # [APO] at each hr from the injection solely from the current dose            
            seaport[j+24*(n_prime*freq_prime+n*(i-4))] = seaport[j+24*(n_prime*freq_prime+n*(i-4))] + dose
            # add to PK from the first dose at the corresponding timepoint

    if i == 7:     # dose beyond priming dose
        D_s = D[len(D)-1]      # repeat injection of 15 mg
        for j in np.array(range(24*(d-n_prime*freq_prime-n*(i-4))+1)):    # dosing time
            dose = (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*j) + f*(1-p)*(e-s)*np.exp(-f*j) + s*p*(e-f)*np.exp(-s*j))
            # [APO] at each hr from the injection solely from the current dose            
            seaport[j+24*(n_prime*freq_prime+n*(i-4))] = seaport[j+24*(n_prime*freq_prime+n*(i-4))] + dose
            # add to PK from the first dose at the corresponding timepoint

    if i == 8:     # dose beyond priming dose
        D_s = D[len(D)-1]      # repeat injection of 15 mg
        for j in np.array(range(24*(d-n_prime*freq_prime-n*(i-4))+1)):    # dosing time
            dose = (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*j) + f*(1-p)*(e-s)*np.exp(-f*j) + s*p*(e-f)*np.exp(-s*j))
            # [APO] at each hr from the injection solely from the current dose            
            seaport[j+24*(n_prime*freq_prime+n*(i-4))] = seaport[j+24*(n_prime*freq_prime+n*(i-4))] + dose
            # add to PK from the first dose at the corresponding timepoint

    if i == 9:     # dose beyond priming dose
        D_s = D[len(D)-1]      # repeat injection of 15 mg
        for j in np.array(range(24*(d-n_prime*freq_prime-n*(i-4))+1)):    # dosing time
            dose = (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*j) + f*(1-p)*(e-s)*np.exp(-f*j) + s*p*(e-f)*np.exp(-s*j))
            # [APO] at each hr from the injection solely from the current dose            
            seaport[j+24*(n_prime*freq_prime+n*(i-4))] = seaport[j+24*(n_prime*freq_prime+n*(i-4))] + dose
            # add to PK from the first dose at the corresponding timepoint

    if i == 10:     # dose beyond priming dose
        D_s = D[len(D)-1]      # repeat injection of 15 mg
        for j in np.array(range(24*(d-n_prime*freq_prime-n*(i-4))+1)):    # dosing time
            dose = (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*j) + f*(1-p)*(e-s)*np.exp(-f*j) + s*p*(e-f)*np.exp(-s*j))
            # [APO] at each hr from the injection solely from the current dose            
            seaport[j+24*(n_prime*freq_prime+n*(i-4))] = seaport[j+24*(n_prime*freq_prime+n*(i-4))] + dose
            # add to PK from the first dose at the corresponding timepoint

    if i == 11:     # dose beyond priming dose
        D_s = D[len(D)-1]      # repeat injection of 15 mg
        for j in np.array(range(24*(d-n_prime*freq_prime-n*(i-4))+1)):    # dosing time
            dose = (D_s/(V*(e-f)*(e-s)))*(((p*f-f)*(e-s)-(p*s)*(e-f))*np.exp(-e*j) + f*(1-p)*(e-s)*np.exp(-f*j) + s*p*(e-f)*np.exp(-s*j))
            # [APO] at each hr from the injection solely from the current dose            
            seaport[j+24*(n_prime*freq_prime+n*(i-4))] = seaport[j+24*(n_prime*freq_prime+n*(i-4))] + dose
            # add to PK from the first dose at the corresponding timepoint

### rApokyn simulation over 30 days
d = int(30)     # 30 days
n = int(3)      # dosing three times a day (3 times per day)
t = int(24/n)   # dosing interval (every 8 hours)
freq = int(d*n)    # total dosing frequency

# create an empty array for [APO] 
rapokyn = np.empty(0, dtype=int)

for i in np.array(range(freq)):  # dosing number index
    if i == 0:   # if first dose
        rapokyn = np.append(rapokyn, ((D_a/(V*(e_a-f_a)*e_a))*((-f_a*e_a)*np.exp(-e_a*time) + f_a*e_a*np.exp(-f_a*time))))
        # update seaport array with PK from the first dose

    else:     # if second dose and beyond
        i = int(i)
        for j in np.array(range(24*d-t*i+1)):    # dosing time
            dose = (D_a/(V*(e_a-f_a)*e_a))*((-f_a*e_a)*np.exp(-e_a*j) + f_a*e_a*np.exp(-f_a*j))
            # [APO] at each hr from the second dose injection solely from the second dose
            rapokyn[j+t*i] = rapokyn[j+t*i] + dose   
            # add to PK from the first dose at the corresponding timepoint


# In[ ]:


# Time spent between MTC and MEC
time_30d = np.array(range(24*30-96+1))
t_btw_seaport = np.zeros(24*30-96+1, dtype=int)
t_btw_rapokyn = np.zeros(24*30-96+1, dtype=int)

for i in range(len(t_btw_seaport)-1):
    if 4 <= seaport[96+i] <= 10 and 4 <= seaport[96+i+1] <= 10:
        t_btw_seaport[i+1] = t_btw_seaport[i] + 1
    else:
        t_btw_seaport[i+1] = t_btw_seaport[i]
        
for i in range(len(t_btw_rapokyn)-1):
    if 4 <= rapokyn[96+i] <= 10 and 4 <= rapokyn[96+i+1] <= 10:
        t_btw_rapokyn[i+1] = t_btw_rapokyn[i] + 1
    else:
        t_btw_rapokyn[i+1] = t_btw_rapokyn[i]

plt.subplot(111).spines["top"].set_visible(False)
plt.subplot(111).spines["right"].set_visible(False)
plt.subplot(111).get_xaxis().tick_bottom() 
plt.subplot(111).get_yaxis().tick_left() 
plt.yticks(np.arange(0, 24*30, 100), fontsize=14, fontweight='bold', fontname = "Arial")  
plt.xticks(np.arange(0, 24*30, 100), fontsize=14, fontweight='bold', fontname = "Arial")
plt.subplot(111).set_xlabel("Time from 96 hr post first injection (hr)", fontsize=14, fontweight='bold', fontname = "Arial")
plt.subplot(111).set_ylabel("Time spent between MTC and MEC (hr)", fontsize=14, fontweight='bold', fontname = "Arial")

plt.plot(time_30d, t_btw_rapokyn, lw=5, color="black", label = 'rApokyn')
plt.plot(time_30d, t_btw_seaport, lw=5, color="red", label = 'SEAPORT')
leg = plt.legend()


# In[ ]:


# Time spent below MEC
time_30d = np.array(range(24*30-96+1))
t_below_seaport = np.zeros(24*30-96+1, dtype=int)
t_below_rapokyn = np.zeros(24*30-96+1, dtype=int)

for i in range(len(t_below_seaport)-1):
        if 4 >= seaport[96+i]:
            t_below_seaport[i+1] = t_below_seaport[i] + 1
        else:
            t_below_seaport[i+1] = t_below_seaport[i]
        
for i in range(len(t_below_rapokyn)-1):
        if 4 >= rapokyn[96+i]:
            t_below_rapokyn[i+1] = t_below_rapokyn[i] + 1
        else:
            t_below_rapokyn[i+1] = t_below_rapokyn[i]

plt.subplot(111).spines["top"].set_visible(False)
plt.subplot(111).spines["right"].set_visible(False)
plt.subplot(111).get_xaxis().tick_bottom() 
plt.subplot(111).get_yaxis().tick_left() 
plt.yticks(np.arange(0, 24*30, 30), fontsize=14, fontweight='bold', fontname = "Arial")  
plt.xticks(np.arange(0, 24*30, 100), fontsize=14, fontweight='bold', fontname = "Arial")
plt.subplot(111).set_xlabel("Time from 96 hr post first injection (hr)", fontsize=14, fontweight='bold', fontname = "Arial")
plt.subplot(111).set_ylabel("Time spent below MEC (hr)", fontsize=14, fontweight='bold', fontname = "Arial")

plt.plot(time_30d, t_below_rapokyn, lw=5, color="black", label = 'rApokyn')
plt.plot(time_30d, t_below_seaport, lw=5, color="red", label = 'SEAPORT')
leg = plt.legend()


# In[ ]:


# Time spent above MTC
time_30d = np.array(range(24*30-96+1))
t_above_seaport = np.zeros(24*30-96+1, dtype=int)
t_above_rapokyn = np.zeros(24*30-96+1, dtype=int)

for i in range(len(t_above_seaport)-1):
        if 10 <= seaport[96+i]:
            t_above_seaport[i+1] = t_above_seaport[i] + 1
        else:
            t_above_seaport[i+1] = t_above_seaport[i]
        
for i in range(len(t_above_rapokyn)-1):
        if 10 <= rapokyn[96+i]:
            t_above_rapokyn[i+1] = t_above_rapokyn[i] + 1
        else:
            t_above_rapokyn[i+1] = t_above_rapokyn[i]
            
plt.subplot(111).spines["top"].set_visible(False)
plt.subplot(111).spines["right"].set_visible(False)
plt.subplot(111).get_xaxis().tick_bottom() 
plt.subplot(111).get_yaxis().tick_left() 
plt.yticks(np.arange(0, 24*30, 30), fontsize=14, fontweight='bold', fontname = "Arial")  
plt.xticks(np.arange(0, 24*30, 100), fontsize=14, fontweight='bold', fontname = "Arial")
plt.subplot(111).set_xlabel("Time from 96 hr post first injection (hr)", fontsize=14, fontweight='bold', fontname = "Arial")
plt.subplot(111).set_ylabel("Time spent above MTC (hr)", fontsize=14, fontweight='bold', fontname = "Arial")

plt.plot(time_30d, t_above_rapokyn, lw=5, color="black", label = 'rApokyn')
plt.plot(time_30d, t_above_seaport, lw=5, color="red", label = 'SEAPORT')
leg = plt.legend()


# In[ ]:




