import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

eos_dft = pd.read_csv("./eos_dft.csv",delimiter=";",decimal=",",header=None,names=["volume","energy"])
eos_gap_lei = pd.read_csv("./eos_gap_lei.csv",delimiter=" ",decimal=".",header=None,names=["volume","energy"])

eos_gap_lei["energy"] = eos_gap_lei["energy"] * 1000

# plot data
fig = plt.figure(figsize=(15.5,7.5))
ax = plt.subplot()

l1=ax.scatter(eos_dft["volume"], eos_dft["energy"], s=40, marker='o', c='blue', alpha=0.7)
l2=ax.scatter(eos_gap_lei["volume"], eos_gap_lei["energy"], s=40, marker='s', c='red', alpha=0.7)

#ax.set_xticks([0.8,1.0,1.2,1.4,1.6,1.8,2.0])

ax.set_xlabel('volume $\AA^3$', fontsize='x-large')
ax.set_ylabel('Energy (meV)', fontsize='x-large')

plt.legend(handles=[l1,l2],labels=['Dragoni_DFT','Lei_GAP'],loc='upper center')
#plt.show()
fig.savefig('eos',dpi=300)
