'''fit light curve with the curtain model, use mcmc to show degeneracies
    and what is acutally constrained, if anything'''

import sys
root = '/Users/grant/astro/doc/rz-psc/'
sys.path.append(root+'rz-psc')
import subprocess
import emcee
import corner
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt

# Define the probability function as likelihood * prior.
def lnprior(par):
    t0, g1, g2, tau, dudt = par
    if t0 > 0 and 0.1 < g1 and 0.1 < g2 and (g1<10 or g2<10):
        return 0.0
#    t0, dur1, dur2, tau = par
#    if t0 > 0 and 0.05 < dur1 and 0.05 < dur2 and 0 < tau < 1:
#        return 0.0
    return -np.inf

def getlc(par1, x):
    par = par1.copy()
    par[4] = 10**par[4]
#    t0, dur1, dur2, tau = par1.copy()
#    dudt = 10.
#    g1 = dur1 * dudt
#    g2 = dur2 * dudt
#    par = np.array([t0,g1,g2,tau,dudt])
    out = subprocess.check_output('/Users/grant/astro/code/c/dippers/dip_gauss1d'+' '+str(len(par)/5)+' '+' '.join(par.astype(np.str))+' '+' '.join(x.astype(np.str)),shell=True)
    out1 = np.fromstring(out,dtype=float,sep=' ')
    model = out1.astype(np.float)
    return model
    
def lnlike(par, x, y, yerr):
    model = getlc(par,x)
    return -0.5*(np.sum( ((y-model)/yerr)**2 ) )/len(x)

def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

def plotchain(chain,fname):
    nwalkers,nrun,ndim = chain.shape
    plt.figure()
    fig,ax = plt.subplots(ndim,sharex=True)
    fig.subplots_adjust(hspace=0)
    for i in range(ndim):
        for j in range(nwalkers):
            ax[i].plot(chain[j,:,i])
    plt.savefig(fname)

# Get data
# grab data and exclude some sources
data = np.genfromtxt(root+"figs/all-lc.txt",delimiter=',',usecols=[0,1],skip_header=1)
date = 2454739.25 #2455894.
wid = 6
ok = (data[:,0] > date-wid) & (data[:,0] < date+wid) & (data[:,1] < 1.05)
x = data[ok,0]
y = data[ok,1]
yerr = np.zeros(len(y))+0.02

# Choose the "true" parameters.
t0_ini = x[np.argmin(y)]
tau_ini = 0.75
g1_ini = 1.
g2_ini = 5.
dudt_ini = np.log10(10)
par = np.array([t0_ini,g1_ini,g2_ini,tau_ini,dudt_ini])
#dur1_ini = 0.1
#dur2_ini = 0.7
#par = np.array([t0_ini,dur1_ini,dur2_ini,tau_ini])

# Find the maximum likelihood value.
chi2 = lambda *args: -2 * lnlike(*args)
result = op.minimize(chi2, par, args=(x, y, yerr))
#t0_ml, g1_ml, g2_ml, tau_ml, dudt_ml = result["x"]
print result["x"]
#print("""Maximum likelihood result:
#    t0   = {0} (truth: {1})
#    g1   = {2} (truth: {3})
#    g2   = {4} (truth: {5})
#    tau  = {6} (truth: {7})
#    dudt = {8} (truth: {9})
#""".format(t0_ml, t0_ini, g1_ml, g1_ini, g2_ml, g2_ini, tau_ml, tau_ini, dudt_ml, dudt_ini))
like = lnlike(result["x"],x,y,yerr)
print("lnlike: {0}",format(like))

# Plot the dataset and the first best fit model.
xmod = np.linspace(x.min(),x.max(),num=500)
ymod = getlc(result["x"],xmod)
fig = plt.plot(x,y,'.')
plt.plot(xmod,ymod)
plt.savefig('mcmc.pdf')

# Set up the sampler.
ndim, nwalkers = len(par), 40
pos = [result["x"] + 1e-2*np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr), threads=20)

# Clear and run the chains.
print("Running MCMC (burnin)...")
pos, prob, state = sampler.run_mcmc(pos, 1000)
sampler.reset()
print("Running MCMC (for real)...")
sampler.run_mcmc(pos, 2000, rstate0=state)
print("Done.")

plotchain(sampler.chain,'mcmc-time.pdf')

# Make the triangle plot.
samples = sampler.chain[:, :, :].reshape((-1, ndim))
fig = corner.corner(samples, labels=["$t_0$", "$HWHM_{in}$", "$HWHM_{out}$", "$tau$", "$log(du/dt)$"])
#fig = corner.corner(samples, labels=["$t_0$", "dur1", "dur2", "$tau$"])
fig.savefig("mcmc-triangle.pdf")

fig = plt.figure(frameon=False)
plt.plot(x, y, '.')
for par in samples[np.random.randint(len(samples), size=100)]:
    plt.plot(xmod, getlc(par,xmod), color="k", alpha=0.1)
fig.savefig("mcmc-samples.pdf")

#exit()

# modify samples to physical parameters
au2km = 1.496e8 # km
mstar = 1.
rstar = 1.
rsun  = 7e5 # Rsun in km
rstar2au = rstar*rsun/au2km
vearth = 2*3.14159265359/365.25 # au/day

psampler = sampler
# du/dt to semi-major axis in au = (v/v_Earth)^2 M_Sun/M_star
psampler.chain[:,:,4] = np.log10( ( vearth/(rstar2au*10**psampler.chain[:,:,4]) )**2 * mstar )
# HWHMs to au
psampler.chain[:,:,1] = psampler.chain[:,:,1] * rstar2au
psampler.chain[:,:,2] = psampler.chain[:,:,2] * rstar2au

psamples = psampler.chain[:, :, :].reshape((-1, ndim))

keep = np.where(10**psamples.T[4]*3.14159/psamples.T[1]> 1.0)

fig = corner.corner(psamples[keep], labels=["$t_0$", "$HWHM_{in} (au)$", "$HWHM_{out} (au)$", "$tau$", "$log(a) (au)$"])#,range=[1.,1,1,1,1])
fig.savefig("mcmc-triangle-phys.pdf")
