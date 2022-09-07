##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Name: Danh-Tai HOANG       Email: hoangdanhtai@gmail.com
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utils import *

##===========================================================================
result_folder = "analysis_results/"

x = np.loadtxt("%sR_mean_patient.txt"%result_folder, dtype="str")
genes = x[:,0]
coef = x[:,1].astype(float)
p = x[:,3].astype(float)
print("genes.shape:", genes.shape)

## find coef_min with p < 0.05:
i = p < 0.05
p0 = p[i]
R0 = coef[i]

n0 = R0.shape[0]
R0_min = np.min(R0)
print("n_genes with p<0.05:", n0)
print("R0_min:", R0_min)

np.savetxt("%sR0_min.txt"%result_folder, np.array([n0, R0_min]), fmt="%s")

#labels = np.load("labels_all.npy")
#print(labels.shape)

#preds = np.load("preds_all.npy")
#print(preds.shape)

#coef, slope, p = compute_coef_slope_p(labels, preds)

#print(sum(coef > 0.4),sum(coef > 0.45),sum(coef > 0.5), sum(coef > 0.55), sum(coef > 0.6))
#print(sum(slope > 0.3),sum(slope > 0.4),sum(slope > 0.45), sum(slope > 0.5), sum(slope > 0.55))
#n_predictable_genes = np.sum(p<0.05)
#print("number of predictable genes:", n_predictable_genes)

#np.savetxt("R_slope_p.txt" , np.array((genes, coef, slope, p)).T, fmt="%s %s %s %s")

#np.savetxt("n_genes_predicted.txt", np.array([n_predictable_genes]), fmt="%s")

##===========================================================================
#coef0= np.linspace(0.4, 0.6, 11, endpoint=True)
#n = np.zeros(len(coef0))
#for i, t in enumerate(coef0):
#    n[i] = sum(coef > t)
#    print(i, t, n[i])

#np.savetxt("n_genes_threshold_patient.txt", np.array((coef0, n)).T, fmt="%s %s")
##===========================================================================
## plot
bins1 = np.linspace(-0.2,0.8,11, endpoint=True)
bins2 = np.linspace(0,0.5,11, endpoint=True)

nx,ny = 2,1
fig, ax = plt.subplots(ny,nx,figsize=(nx*3.5,ny*3))

ax[0].hist(coef,bins=bins1,histtype='bar',rwidth=0.9)
ax[0].set_xlabel("Correlation")
ax[0].set_ylabel("Number of genes")
ax[0].plot([R0_min, R0_min],[0,7000], "--", label="p-value=0.05")
ax[0].set_xticks([-0.2,0,0.2,0.4,0.6,0.8])
ax[0].set_ylim(0,7000)
ax[0].set_yticks([2000,4000,6000])

ax[1].hist(p,bins=bins2,histtype='bar',rwidth=0.95)
ax[1].set_xlabel("p-value")
ax[1].set_ylabel("Number of genes")
ax[1].set_xticks([0,0.1,0.2,0.3,0.4,0.5])

ax[0].legend()
plt.tight_layout(h_pad=1, w_pad= 1.5)
plt.savefig("%sR_p_patient.pdf"%result_folder, format='pdf', dpi=50)

##===========================================================================

bins1 = np.linspace(-0.2,0.8,11, endpoint=True)

nx,ny = 1,1
fig, ax = plt.subplots(ny,nx,figsize=(nx*3.5,ny*3))

ax.hist(coef,bins=bins1,histtype='bar',color="lightblue",edgecolor="black",rwidth=0.85)
ax.plot([R0_min, R0_min],[0,7000], "--", color = "red", label="p-value=0.05")
ax.set_xlabel("Correlation")
ax.set_ylabel("Number of genes")
ax.set_xticks([-0.2,0,0.2,0.4,0.6,0.8])
ax.set_ylim(0,7000)
ax.set_yticks([2000,4000,6000])

ax.legend()
plt.tight_layout(h_pad=1, w_pad= 1.5)
plt.savefig("%sR_patient.pdf"%result_folder, format='pdf', dpi=50)

##===========================================================================

print(" --- completed ---")




