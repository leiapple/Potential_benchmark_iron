import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pprint import pprint

eV_to_Ry=1/13.60569301
Ry_to_eV= 13.60569301

eos_dft = pd.read_csv("../REF_DATA/eos_dft.csv",delimiter=";",decimal=",",header=None,names=["volume","energy"])
eos_gap18 = pd.read_csv("../REF_DATA/eos_gap18.csv",delimiter=" ",decimal=".",header=None,names=["volume","energy"])
eos_gap21 = pd.read_csv("../data/eos_mlip.csv",delimiter=" ",decimal=".",header=None,names=["volume","energy"])
min_gap18=min(eos_gap18["energy"])
min_gap21=min(eos_gap21["energy"])
eos_dft["energy"] = eos_dft["energy"] -min(eos_dft["energy"]) 
eos_gap18["energy"] = (eos_gap18["energy"] -min(eos_gap18["energy"])) * 1000
eos_gap21["energy"] = (eos_gap21["energy"] -min(eos_gap21["energy"])) * 1000

bp_dft = pd.read_csv("../REF_DATA/BainPath_DFT.csv",delimiter=",",decimal=".",header=None,names=["ratio","energy"])
bp_gap18 = pd.read_csv("../REF_DATA/BainPath_GAP.csv",delimiter=",",decimal=".",header=None,names=["ratio","energy"])
bp_gap21 = pd.read_csv("../data/bain_path.csv",delimiter=" ",decimal=".",header=None,names=["ratio","energy","ratio2","lx","ly","lz"])
bp_dft["energy"] = bp_dft["energy"] -min(bp_dft["energy"])
bp_gap18["energy"] = bp_gap18["energy"] -min(bp_gap18["energy"])
bp_gap21["energy"] = (bp_gap21["energy"] -min(bp_gap21["energy"]))*1000/2

#db1 = pd.read_csv("./db1_ca_ratio.dat",delimiter=" ",decimal=".",header=None,names=["i","ratio","energy"])
#db9 = pd.read_csv("./db9_ca_ratio.dat",delimiter=" ",decimal=".",header=None,names=["i","ratio","energy"])
#db1["energy"] = (db1["energy"] - min_gap18)*1000
#db9["energy"] = (db9["energy"] - min_gap21)*1000

fig, (ax, bx) = plt.subplots(nrows=1, ncols=2, figsize=(15,7))
plt.rcParams['font.size'] = '20'
ax.grid(c='gainsboro',ls='--',lw=0.7)
bx.grid(c='gainsboro',ls='--',lw=0.7)

# plot EV curve
ax.set_title('Energy volume curve',fontsize=20)
ax.scatter(eos_dft["volume"], eos_dft["energy"], s=30, marker='o', c='dodgerblue', alpha=1,label='DFT')
ax.scatter(eos_gap21["volume"], eos_gap21["energy"],  s=30, marker='s', c = 'red', alpha=1,label='New_IAP')
ax.plot(eos_dft["volume"], eos_dft["energy"],lw=1.5,alpha=1)
ax.plot(eos_gap21["volume"], eos_gap21["energy"],lw=1.5,alpha=1)

ax.set_xlabel('Volume, [$\AA^3$]',fontsize='20')
ax.set_ylabel('Energy, [meV]',fontsize='20')
ax.legend(loc='upper center',fontsize=20,markerscale=1.5)


# Plot right side
bx.set_title('Bain path',fontsize=20)
bx.scatter(bp_dft["ratio"],bp_dft["energy"],s=30, marker='o', c='dodgerblue', alpha=1,label='DFT')
bx.scatter(bp_gap21["ratio"], bp_gap21["energy"],  s=30, marker='s', c = 'red', alpha=1,label='New_IAP')
bx.plot(bp_dft["ratio"], bp_dft["energy"],lw=1.5,alpha=1)
bx.plot(bp_gap21["ratio"], bp_gap21["energy"],lw=1.5,alpha=1)
#bx.scatter(db1["ratio"],db1["energy"], s=5, marker = '+', c='darkgrey',alpha=0.4) #,label='DB1' )
#bx.scatter(db9["ratio"],db9["energy"], s=5, marker = '*', c='pink',alpha=0.4) #,label='DB10' )

plt.ylim(0,400)

bx.set_xlabel('c/a ratio',fontsize='20')
bx.legend(loc='upper center',fontsize=20,markerscale=1.5)

for label in (ax.get_xticklabels() + ax.get_yticklabels() + bx.get_xticklabels() + bx.get_yticklabels()):
	label.set_fontsize(20)

fig.savefig('eos_bp.png',dpi=300)
