import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Read the data for <111>{110} stacking fault energy.
sfe_110_dft = pd.read_csv("../REF_DATA/sfe_110_dft.csv",delimiter=",",decimal=".",header=0,names=["dis","energyeV"])
sfe_110_dft2 = pd.read_csv("../REF_DATA/110_111.csv",delimiter=",",decimal=".",header=0,names=["dis","energyJ"])
sfe_110_gap21 = pd.read_csv("../data/sfe_110.csv",delimiter=" ",decimal=".",header=0,names=["dis","energyeV","energyJ"])

# Change energy units from eV/Å^2 to J/m^2.
# Read the data for <112>{110} stacking fault energy.
sfe_112_dft1 = pd.read_csv("../REF_DATA/dft2_112_111.csv",delimiter=",",decimal=".",header=0,names=["dis","energyeV"])
sfe_112_dft2 = pd.read_csv("../REF_DATA/112_111.csv",delimiter=",",decimal=".",header=0,names=["dis","energyJ"])
sfe_112_gap21 = pd.read_csv("../data/sfe_112.csv",delimiter=" ",decimal=".",header=0,names=["dis","energyeV","energyJ"])

# plot data
fig, (ax, bx) = plt.subplots(nrows=1, ncols=2, figsize=(15,8))
plt.rcParams['font.size'] = '20'

# Plot for <111>{110}.
ax.grid(c='gainsboro',ls='--',lw=0.7)
ax.set_title('<111>{110}',fontsize=20)
ax.scatter(sfe_110_dft["dis"]*2.834*np.sqrt(3)/2, sfe_110_dft["energyeV"] * 0.001 * 16.0217733, s=60, marker='s', c='red', label='DFT')
ax.scatter(sfe_110_dft2["dis"]*2.834*np.sqrt(3)/2, sfe_110_dft2["energyJ"], s=60, marker='v', c='blue', label='DFT2')
ax.scatter(sfe_110_gap21["dis"], sfe_110_gap21["energyJ"], s=60, marker='o', c='orange', label='ML-IAP')
ax.plot(sfe_110_gap21["dis"], sfe_110_gap21["energyJ"], c='orange',lw=2,ls='--')

ax.set_xlim([0,2.45])
ax.set_xlabel('Displacements, [$\AA$]', fontsize=20)
ax.set_ylabel('Stacking Fault Energy, [J/m$^2$])', fontsize=20)
ax.legend(loc='upper right',fontsize=20)
ax.text(0.6, 0.2, 'IAP='+str(np.round(max(sfe_110_gap21["energyJ"]),4)), transform=ax.transAxes, fontsize='large', horizontalalignment='right',
            verticalalignment='bottom')

# Plot for <111>{112}
bx.grid(c='gainsboro',ls='--',lw=0.7)
bx.set_title('<111>{112}',fontsize=20)
bx.scatter(sfe_112_dft1["dis"]*2.834*np.sqrt(3)/2, sfe_112_dft1["energyeV"]* 0.001 * 16.0217733, s=60, marker='s', c='red', label='DFT1')
bx.scatter(sfe_112_dft2["dis"]*2.834*np.sqrt(3)/2, sfe_112_dft2["energyJ"], s=60, marker='v', c='blue', label='DFT2')
bx.scatter(sfe_112_gap21["dis"], sfe_112_gap21["energyJ"], s=60, marker='o', c='orange', label='ML-IAP')
bx.plot(sfe_112_gap21["dis"], sfe_112_gap21["energyJ"],c='orange',lw=1.5,ls='--')
bx.text(0.6, 0.2, 'IAP='+str(np.round(max(sfe_112_gap21["energyJ"]),4)), transform=bx.transAxes, fontsize='large', horizontalalignment='right',
            verticalalignment='bottom')

bx.set_xlim([0,2.45])
bx.set_xlabel('Displacements, [$\AA$]', fontsize=20)
bx.set_ylabel('Stacking Fault Energy, [J/m$^2$])', fontsize=20)
bx.legend(loc='upper right',fontsize=20)

f = open('../data/results.txt','a')
f.write('=========================================\n')
f.write('Max (110) SF energy\n')
f.write(str(max(sfe_110_gap21["energyJ"]))+'\n')
f.write('Max (112) SF energy\n')
f.write(str(max(sfe_112_gap21["energyJ"]))+'\n')
#plt.show()
fig.savefig('sfe.png',dpi=300)
