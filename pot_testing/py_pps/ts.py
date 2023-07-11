import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

eVA32GPa = 160.2176621

# load dft data
dft_100 = np.genfromtxt('../REF_DATA/ts_100_dft.csv')
dft_110 = np.genfromtxt('../REF_DATA/ts_110_dft.csv')

ts_100 = np.genfromtxt('../data/ts_100.csv')
ts_110 = np.genfromtxt('../data/ts_110.csv')

fig, ((ax,bx)) = plt.subplots(nrows=1, ncols=2, figsize=(16,8))
c1 = '#c1272d'
c2 = '#0000a7'
c3 = '#eecc16'
c4 = '#008176'
c5 = '#b3b3b3'
import matplotlib
font = {'family' : 'helvetica',
'weight' : 'normal',
'size' : 20}
matplotlib.rc('font', **font)

ax.plot(ts_100[1:,0]-0.025, eVA32GPa * (ts_100[1:,2]-ts_100[:-1,2])/(0.05*ts_100[1:,3]),label='ML-IAP', c=c1, lw=3)
bx.plot(ts_110[1:,0]-0.025, eVA32GPa * (ts_110[1:,2]-ts_110[:-1,2])/(0.05*ts_110[1:,3]),label='ML-IAP', c=c1, lw=3)

ax.scatter(dft_100[:,0], dft_100[:,1], facecolors='none', edgecolors=c3, s=80)
ax.plot(dft_100[:,0], dft_100[:,1], c=c3, ls="--", lw=2.5,label='DFT')

bx.scatter(dft_110[:,0], dft_110[:,1], facecolors='none', edgecolors=c3,s=80)
bx.plot(dft_110[:,0], dft_110[:,1], c=c3, ls="--", lw=2.5,label='DFT')

ax.legend()
ax.hlines(0,0,5,ls=':',color='grey',lw=2)
bx.hlines(0,0,5,ls=':',color='grey',lw=2)

for pplot in [ax, bx]:
    pplot.set_xlim(0,5)
    pplot.set_ylim(-2.5,37.5)
ax.set_xlabel('Separation distance, (Å)', weight='bold' )
bx.set_xlabel('Separation distance, (Å)', weight='bold' )
ax.set_ylabel('Normal stress, (GPa)',weight='bold' )
bx.set_ylabel('Normal stress, [GPa]', weight='bold' )

ax.annotate("(a)", xy=(-0.17, 1), weight='bold', xycoords="axes fraction", fontname='Helvetica',fontsize=26)
bx.annotate("(b)", xy=(-0.17, 1), weight='bold', xycoords="axes fraction", fontname='Helvetica',fontsize=26)

#-----
f = open('../data/results.txt','a')
f.write('=========================================\n')
f.write('Max traction along (100) (GPa)\n')
f.write(str(max(eVA32GPa * (ts_100[1:,2]-ts_100[:-1,2])/(0.05*ts_100[1:,3])))+'\n')
f.write('Max traction along (110) (GPa)\n')
f.write(str(max(eVA32GPa * (ts_110[1:,2]-ts_110[:-1,2])/(0.05*ts_110[1:,3])))+'\n')

plt.savefig("ts.png", format='png', dpi=600)

