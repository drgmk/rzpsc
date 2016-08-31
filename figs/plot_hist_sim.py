import sys
root = '/Users/grant/astro/doc/rz-psc/'
sys.path.append(root+'rz-psc')
from funcs import *

# grab data and exclude some sources
data = ascii.read(root+"figs/all-lc.txt")
reject = (data['obs'] != 'WASP') & (data['obs'] != 'KELT')
data = data[np.invert(reject)]
data = data[data['flux'] < 1.15]

nbins = 160
bins = 160.0 * np.arange(nbins)/(nbins-1)

fig,ax = plt.subplots(3,1,sharex=True,figsize=(8,6))
fig.subplots_adjust(hspace=0)
ax[2].set_xlabel('Period (days)')
ax[0].set_ylim([0,34])
#ax[0].set_ylabel('number')
ax[1].set_ylim([0,34])
ax[1].set_ylabel('number')
ax[2].set_ylim([0,34])
#ax[2].set_ylabel('number')

histreal,dtall,tall,minall,gapall = get_power(data,bins)
ax[0].plot(bins,histreal,alpha=1,c='blue')
ax[1].plot(bins,histreal,alpha=1,c='blue')
ax[2].plot(bins,histreal,alpha=1,c='blue')

nrun = 500

hmean,hstd = simulate_single_dips(data,bins,nrun=nrun,ndip=295,flux=0.1,width=1.0,dwidth=3.0)
ax[0].plot(bins,hmean,'--',alpha=1,c='black')
ax[0].plot(bins,hmean+hstd,'--',alpha=1,c='black')
ax[0].plot(bins,hmean-hstd,'--',alpha=1,c='black')

hmean,hstd = simulate_rep_dips(data,bins,nrun=nrun,ndip=100,plow=60.0,phigh=70.0,nrep=3,flux=0.1,width=1.0,dwidth=3.0)
ax[1].plot(bins,hmean,'--',alpha=1,c='black')
ax[1].plot(bins,hmean+hstd,'--',alpha=1,c='black')
ax[1].plot(bins,hmean-hstd,'--',alpha=1,c='black')

hmean,hstd = simulate_rep_dips(data,bins,nrun=nrun,ndip=100,plow=64.0,phigh=65.0,nrep=3,flux=0.1,width=1.0,dwidth=3.0)
ax[2].plot(bins,hmean,'--',alpha=1,c='black')
ax[2].plot(bins,hmean+hstd,'--',alpha=1,c='black')
ax[2].plot(bins,hmean-hstd,'--',alpha=1,c='black')


fig.savefig('hist_single.png',bbox_inches = 'tight')
