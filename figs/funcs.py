from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.stats import sigma_clip
from datetime import datetime

"""Function to take lowest pointi in a light curve and collect
    all the points around it as part of a single occultation.
    Return a new light curve with the occultation excluded, and
    the lower, minimum, and upper times for the event"""
def clipit(x,y,sig=0.1, mean=1.0, xsig=1.0):
#    plt.scatter(x, y)
    
    # find the lowest point
    ind = y.argmin(axis=0)
#    plt.scatter(x[ind], y[ind], c='r')

    # shrink to the lower side, until you get flux within (mean - sig)
    # or hits the zeroth index. imin is the index of the highest value
    # we will keep from the lower half
    imin = ind
    while(y[imin] < (mean-xsig*sig)):
        imin = imin - 1
        if imin < 0:
            break

    # grow to the upper side, until you get flux within (mean - sig)
    # or hits the array size. imax is the index of the lowest value
    # we will keep from the upper half
    imax = ind
    while(y[imax] < (mean-xsig*sig)):
        imax = imax + 1
        if imax == y.size:
            break

#    print(ind)
#    print(imin)
#    print(imax)
#    plt.scatter(x[imax],y[imax],c='orange')

    # data to keep, lower and upper sides
    indip = np.ones(len(x),dtype=bool)
    #     indip[imin+1:imax] = True # data in the dip
    #     new_x = x[np.invert(indip)]
    #     new_y = y[np.invert(indip)]
    #     indip[imin] = True # add end points
    #     indip[imax] = True
    #     in_x = x[indip]
    if imin >= 0:
        low_x = x[:imin+1]
        low_y = y[:imin+1]
        indip[:imin] = False # keep last of lower chunk
    else:
        low_x = []
        low_y = []
        indip[0] = False
    if imax < y.size:
        hig_x = x[imax:]
        hig_y = y[imax:]
        indip[imax+1:] = False # keep first of upper chunk
    else:
        hig_x = []
        hig_y = []
        indip[-1] = False
    new_x = np.concatenate((low_x,hig_x))
    new_y = np.concatenate((low_y,hig_y))
#    plt.scatter(new_x,new_y,c='yellow')
    if imin < 0:
        xmin = np.nan
    else:
        xmin = x[imin]
    if imax == y.size:
        xmax = np.nan
    else:
        xmax = x[imax]

    in_x = x[indip]
    gap = 0
    if len(in_x) > 1:
        dt = in_x[1:]-in_x[:len(in_x)-1]
        gap = np.max(dt)

    return(new_x, new_y, xmin, x[ind], xmax, y[ind], gap)


def dbox(ax1, tmin, tmid, tmax, mid = 1.0, wid = 0.1, alpha=0.1):
    import matplotlib.patches as patches
    ax1.add_patch(
                  patches.Rectangle(
                  (tmin, mid-wid), # (x,y)
                  (tmax-tmin),     # width
                  2*wid,           # height
                  alpha=alpha, facecolor='red', edgecolor='none')
                  )
    ax1.scatter(tmid, mid, color='red')

"""Function to figure mean and std, and then to iteratively get dips until there aren't any left
   threshold is the number of sigma data need to be below to be considered a dip"""
def get_dips(xx,yy,threshold=4.,xsig=1.0):
    
    # mean and std
    ok = np.where((yy > 0.7) & (yy < 1.2)) # remove obvious outliers
    clip_masked = sigma_clip(yy[ok],sigma=2,iters=10,cenfunc=np.mean)
    std = np.std(yy[ok][~clip_masked.mask])
    mean = np.mean(yy[ok][~clip_masked.mask])
    
    tdip = []
    ydip = []
    gap = []
    while min(yy) < mean-threshold*std:
        
        (nx, ny, tmin, tmid, tmax, ymin, gap1) = clipit(xx,yy,mean=mean,sig=std,xsig=xsig)
        xx = nx
        yy = ny
        # make or add to array of dips
        ydip = np.append(ydip,ymin)
        gap = np.append(gap,gap1)
        try:
            tdip = np.vstack((tdip,[tmin, tmid, tmax]))
        except ValueError:
            tdip = np.array([[tmin, tmid, tmax]])

    # if there were dips, sort into increasing order of central time
    if len(tdip) == 0:
        return (np.nan,np.nan,np.nan,mean,std)
    elif len(tdip) == 1:
        return (tdip,ydip,gap,mean,std)
    else:
        srt = np.argsort(tdip[:,1])
        return (tdip[srt],ydip[srt],gap[srt],mean,std)

# fuction to get time differences between dips
def get_dtdip(tdip):
    
    if not (type(tdip) == np.ndarray):
        return np.nan
    
    # sort into increasing order of central time
    srt = np.argsort(tdip[:,1])
    tdip = tdip[srt]

    dtdip = []
    for i in range(len(tdip)-1):
        for j in range(len(tdip)):
            if j <= i:
                continue
            dtmin = tdip[j,0]-tdip[i,2] # [i,2] to include width of lower dip
            dtmid = tdip[j,1]-tdip[i,1] # [i,1]   (all [i,1] to ignore)
            dtmax = tdip[j,2]-tdip[i,0] # [i,0]
            try:
                dtdip = np.vstack((dtdip,[dtmin, dtmid, dtmax]))
            except ValueError:
                dtdip = np.array([[dtmin, dtmid, dtmax]])

    # if there were dips, sort into increasing order of central time
    if len(dtdip) == 0:
        return np.nan
    elif len(dtdip) == 1:
        return dtdip
    else:
        srt = np.argsort(dtdip[:,1])
        return dtdip[srt]

def get_power(data,bins,plot=0):
    
    if plot:
        fig, ax = plt.subplots(figsize=(16,10))
    year = 1900
    j = 11
    dtall = []
    tall = []
    minall = []
    gapall = []
    while year < 2015:
        
        # get this year's data
        t0 = Time(datetime(year, 5, 1, 0, 0, 0)) # start of decent data
        t0.format = 'jd'
        ok = np.where((data['JD'] > t0.value) & (data['JD'] < t0.value+365.25))
        if len(ok[0]) == 0:
            year += 1
            continue
        t = data['JD'][ok] - t0.value
        f = data['flux'][ok]
        
        tdip,ydip,gap,mean,std = get_dips(t,f,threshold=6.0,xsig=1.0)
        dtdip = get_dtdip(tdip)
        
        if type(tdip) == np.ndarray:
            minall = np.append(minall,ydip)
            gapall = np.append(gapall,gap)
            try:
                tall = np.vstack((tall,tdip+t0.value))
            except ValueError:
                tall = tdip+t0.value

        # show light curves and detected dips, time differences between dips
        if plot:
            off = 160.
            offd = 0.
            tdip += off
            ax.plot(t+off,0.5*f+0.5+(j-1),'.')
            ax.text(off+10,j-0.05,year)
            if type(tdip) == np.ndarray:
                for i in range(len(tdip)):
                    dbox(ax, tdip[i,0], tdip[i,1], tdip[i,2], mid=j, wid=0.3)

            if type(dtdip) == np.ndarray:
                for i in range(len(dtdip)):
                    dbox(ax, dtdip[i,0]+offd, dtdip[i,1]+offd, dtdip[i,2]+offd, mid=j, wid=0.4, alpha=0.1)
                    dbox(ax, dtdip[i,0]+offd, dtdip[i,1]+offd, dtdip[i,2]+offd, mid=1, wid=0.4, alpha=0.04)
            
        if type(dtdip) == np.ndarray:
            try:
                dtall = np.vstack((dtall,dtdip))
            except ValueError:
                dtall = dtdip
    
        j -= 1
        year += 1

#    dtall = get_dtdip(tall) # dts for all, including inter-year dips

    nbin = len(bins)
    hist = np.zeros(nbin)
    for i in range(len(dtall)):
        ini = np.where((bins > dtall[i,0]) & (bins < dtall[i,2]))
        hist[ini] += 1

    if plot:
        ax.text(off+10,0.95,'all')
        ax.yaxis.set_visible(False)
        ax.set_xlabel('time (days)')
        ax.set_ylim([0.5,11.5])
        ax.set_xlim([0,475])
        fig.savefig('lc_dips_dts.pdf',bbox_inches = 'tight')

    return (hist,dtall,tall,minall,gapall)

'''functions to simulate dips'''
def add_1dip(datain,time,flux=0.1,width=1.0):
    data = datain
    ini = np.where((data['JD'] > time-width/2.) & (data['JD'] < time+width/2.))
    if len(ini[0]) > 0:
        data['flux'][ini] = flux
    return data

def add_rep_dip(datain,time,period=80.0,flux=0.1,nrep=3,width=1.0):
    data = datain
    for i in np.arange(nrep):
        data = add_1dip(data,time+i*period,flux=flux,width=width)
    return data

def add_many_single_dips(datain,ndip=100,flux=0.1,width=1.0,dwidth=2.0):
    data = datain
    data['flux'] = np.random.normal(loc=1,scale=0.05,size=len(data))
    for i in np.arange(ndip):
        time = np.random.uniform(low=np.min(data['JD']),high=np.max(data['JD']))
        rwidth = np.random.uniform(low=width,high=width+dwidth)
        data = add_1dip(data,time,flux=flux,width=rwidth)
    return data

def add_many_rep_dips(datain,plow=60.0,phigh=80.0,ndip=100,nrep=3,flux=0.1,width=1.0,dwidth=2.0):
    data = datain
    data['flux'] = np.random.normal(loc=1,scale=0.05,size=len(data))
    for i in np.arange(ndip):
        time = np.random.uniform(low=np.min(data['JD']),high=np.max(data['JD']))
        period = np.random.uniform(low=plow,high=phigh)
        rwidth = np.random.uniform(low=width,high=width+dwidth)
        data = add_rep_dip(data,time,width=rwidth,period=period,nrep=nrep)
    return data

def simulate_single_dips(datain,bins,nrun=50,ndip=1000,flux=0.1,width=1.0,dwidth=2.0):
    data = datain
    histall = []
    for j in np.arange(nrun):
        fake = add_many_single_dips(data,ndip=ndip,flux=flux,width=width,dwidth=dwidth)
        hist,dtall,tall,minall,gapall = get_power(fake,bins,plot=0)
        try:
            histall = np.vstack((histall,hist))
        except ValueError:
            histall = hist
    histmean = np.mean(histall,axis=0)
    histstd = np.std(histall,axis=0)
    return (histmean,histstd)

def simulate_rep_dips(datain,bins,nrun=50,ndip=1000,flux=0.1,width=1.0,dwidth=2.0,plow=60.0,phigh=80.0,nrep=3):
    data = datain
    histall = []
    for j in np.arange(nrun):
        fake = add_many_rep_dips(data,plow=plow,phigh=phigh,nrep=nrep,ndip=ndip,flux=flux,width=width,dwidth=dwidth)
        hist,dtall,tall,minall,gapall = get_power(fake,bins,plot=0)
        try:
            histall = np.vstack((histall,hist))
        except ValueError:
            histall = hist
    histmean = np.mean(histall,axis=0)
    histstd = np.std(histall,axis=0)
    return (histmean,histstd)


'''Function to return Fmin and dt, similar approach to dm-dt'''
def dfdt(data):
    year = 1850
    df = []
    dt = []
    while year < 2015:
        # get this year's data
        t0 = Time(datetime(year, 5, 1, 0, 0, 0)) # start of decent data
        t0.format = 'jd'
        for thresh in 0.9*np.arange(50)/49.:
            ok = (data['JD'] > t0.value) & (data['JD'] < t0.value+365.25) & (data['flux'] < thresh)
            if np.any(ok) == False:
                continue
            t = data['JD'][ok]
            f = data['flux'][ok]
            n = len(t)
            for i in range(n):
                for j in 1+i+np.arange(n-i-1):
                    df.append( thresh )
                    dt.append( np.fabs(t[i]-t[j]) )
        year += 1
    return np.array(dt),np.array(df)
