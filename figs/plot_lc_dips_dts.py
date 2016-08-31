import sys
root = '/Users/grant/astro/doc/rz-psc/'
sys.path.append(root+'rz-psc')
from funcs import *

# grab data and exclude some sources
data = ascii.read(root+"figs/all-lc.txt")
reject = (data['obs'] != 'WASP') & (data['obs'] != 'KELT')
#reject = (data['obs'] == 'AAVSO') & (data['JD'] < 2456413.5)
data = data[np.invert(reject)]
data = data[data['flux'] < 1.15]
#data = data[data['JD'] > 2448012.5]

nbins = 1000
bins = 4000.0 * np.arange(nbins)/(nbins-1)
histreal,dtall,tall,minall,gapall = get_power(data,bins,plot=1)

def fmin(dt,alpha=1.,k=1.):
    return 1.-k/dt**alpha

wid = tall[:,2]-tall[:,0]
ok = (wid > 0.5)
wid = wid[ok]
minall = minall[ok]
gapall = gapall[ok]

fig,ax = plt.subplots()
ax.scatter(wid,minall,20*(10/gapall))
ok = (gapall < 2.)
ax.plot(wid[ok],minall[ok],'ro')
x = np.arange(3,22,0.1)
ax.plot(x,fmin(x,0.5,np.sqrt(3)))
ax.plot(x,fmin(x,1,3.))
ax.plot(x,fmin(x,2.,9.))
ax.set_xlim([0,22])
ax.set_ylim([-0.1,1])
ax.set_xlabel('Event duration (days)')
ax.set_ylabel('Minimum flux')
fig.savefig('fmin-wid.png',bbox_inches = 'tight')
