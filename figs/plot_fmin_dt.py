import sys
root = '/Users/grant/astro/doc/rz-psc/'
sys.path.append(root+'rz-psc')
from funcs import *

# grab data and exclude some sources
data = ascii.read(root+"figs/all-lc.txt")
reject = (data['obs'] != 'WASP') & (data['obs'] != 'KELT')
data = data[np.invert(reject)]
data = data[data['flux'] < 1.15]

dt,df = dfdt(data)

fig,ax = plt.subplots(figsize=(8,6))
fig1,ax1 = plt.subplots(10,1,figsize=(8,16),sharex=True)
fig1.subplots_adjust(hspace=0)
year = 2000
i = 0
while year < 2015:
    # get this year's data
    t0 = Time(datetime(year, 5, 1, 0, 0, 0)) # start of decent data
    t0.format = 'jd'
    ok = (data['JD'] > t0.value) & (data['JD'] < t0.value+365.25)
    if np.any(ok) == True:
        dt,df = dfdt(data[ok])
        ax.plot(dt,df,'.')
        ax1[i].plot(dt,df,'.')
        ax1[i].set_ylim([0.,0.9])
        i += 1
    year += 1

ax.set_xlabel('time between measurements below $F_{min}$')
ax.set_ylabel('$F_{min}$')
fig.savefig('fmin_dt.png',bbox_inches = 'tight')

ax1[9].set_xlabel('time between measurements below $F_{min}$')
fig1.savefig('fmin_dt_panel.png',bbox_inches = 'tight')
