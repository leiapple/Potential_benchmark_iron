import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

bainpath_dft = pd.read_csv("./BainPath_DFT.csv",delimiter=",",decimal=".",header=None,names=["ratio","energy"])
bainpath_gap = pd.read_csv("./BainPath_GAP.csv",delimiter=",",decimal=".",header=None,names=["ratio","energy"])
bainpath_gap_lei = pd.read_csv("./BainPath_GAP_lei.csv",delimiter=" ",decimal=".",header=None,names=["ratio","energy","ratio2","lx","ly","lz"])

bainpath_dft["energy"] = bainpath_dft["energy"]-min(bainpath_dft["energy"])
bainpath_gap["energy"] = bainpath_gap["energy"]-min(bainpath_gap["energy"])
bainpath_gap_lei["energy"] = (bainpath_gap_lei["energy"]-min(bainpath_gap_lei["energy"]))*1000 /2 

"""
# read my gap calculation data
fname = 'BainPath_GAP_lei.csv'
f = open(fname, 'rt')

# read data from file
print()
print("Data read from file:")
ratio = []
ene = []
while True:
    line = f.readline().strip()
    if line == '': break
    if line[0] == '#' or line[0] == '!': continue
    v, e = [float(x) for x in line.split()[:2]]
    
    ratio.append(v)
    ene.append(float(e))
    print(v, e)
print()
f.close()

ratio = np.array(ratio)
ene = np.array(ene)

ene = (ene -min(ene))*1000/2
"""
# plot data
fig = plt.figure(figsize=(15.5,7.5))
ax = plt.subplot()

l1=ax.scatter(bainpath_dft["ratio"], bainpath_dft["energy"], s=40, marker='o', c='gold', alpha=0.7)
l2=ax.scatter(bainpath_gap["ratio"], bainpath_gap["energy"], s=40, marker='h', c='blue', alpha=0.7)
l3=ax.scatter(bainpath_gap_lei["ratio"], bainpath_gap_lei["energy"], s=40, marker='s', c='red', alpha=0.7)

ax.set_xticks([0.8,1.0,1.2,1.4,1.6,1.8,2.0])

ax.set_xlabel('c/a ratio', fontsize='x-large')
ax.set_ylabel('Energy (meV)', fontsize='x-large')

plt.legend(handles=[l1,l2,l3],labels=['Dragoni_DFT','Dragoni_GAP','Lei_GAP'],loc='upper center')
#plt.show()
fig.savefig('db1_thesis',dpi=300)
