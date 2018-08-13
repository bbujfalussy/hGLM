import numpy as np
import matplotlib.pyplot as plt

#----------------------------------------------------------------------------
# Plotting functions
# plot only the dendritic and somatic response - different trials on the same subplot
def plotTraces(data,xmax=1000):
    plt.figure(figsize=(16,8))
    cols = ['b', 'g', 'r', 'c']
    ax = plt.subplot(2, 1, 1)
    nReps = len(data.vdata)
    for i_rep in np.arange(0,nReps):
        trace = data.vdata[i_rep]
        ax.plot(data.taxis, trace, color=cols[i_rep % 4])
        plt.xlim(0,xmax)
#        plt.ylim(-80,0)

    ax = plt.subplot(2, 1, 2)
    if data.recordDend:
        for i_rep in np.arange(0,nReps):
            trace = data.vDdata[i_rep]
            Dtrace = trace[0]
            ax.plot(data.taxis, Dtrace, color=cols[i_rep % 4])
            plt.xlim(0,xmax)
 #           plt.ylim(-80,0)
    plt.show(block=False)


# plot the stimulus + somatic and dendritic response
# each trial is a different row
def plotResp(data,xmax=1000):
    plt.figure(figsize=(16,8))
    cols = ['b', 'g', 'r', 'c', 'm']
    # if np.ndim(data.vdata) == 2:
    #     for trace in data.vdata:
    #         plt.plot(data.taxis, trace)
    # else:
    nReps = len(data.vdata)
    nPlots = nReps
    nsyn = sum(data.Ensyn)
    if (data.GABA == True): nIsyn = sum(data.Insyn)

    mat_plots = np.arange(3*nPlots).reshape(3, nPlots)

    for i_rep in np.arange(0,nReps):
        i_plot = i_rep 
        ax = plt.subplot(3, nPlots, mat_plots[0,i_plot]+1)
        ind_etimes = (data.stim[:,0] == i_rep + 1)
        if (sum(ind_etimes) > 0) :
            etimes = np.array(data.stim[ind_etimes,:])
            if (len(etimes.shape) == 1) :
                etimes = np.array([data.stim[ind_etimes,:]])
            
            etimes = etimes[:,(1,2)]
            et = etimes[etimes[:,1] < xmax,:]
            if (len(et) > 0):
                ax.plot(et[:,1], et[:,0] + 0.5, '|', color='r')

        ind_itimes = (data.stim[:,0] == -1 * i_rep - 1)
        if (sum(ind_itimes) > 0) :
            itimes = np.array(data.stim[ind_itimes,:])
            if (len(itimes.shape) == 1):
                itimes = np.array([data.stim[ind_itimes,:]])
            itimes = itimes[:,(1,2)]
            it = itimes[itimes[:,1] < xmax,:]
            if (len(it) > 0) :
                ax.plot(it[:,1], it[:,0] + nsyn + 0.5, '|', color='b')

        plt.xlim(0,xmax)
        ymax = nsyn
        if (data.GABA == True): ymax = ymax + nIsyn
        plt.ylim(0,ymax)


        
        ax = plt.subplot(3, nPlots, mat_plots[1,i_plot]+1)
        trace = data.vdata[i_rep]
        ax.plot(data.taxis, trace)
        plt.xlim(0,xmax)
        plt.ylim(-80,0)

        if data.recordDend:
            trace = data.vDdata[i_rep]
            nDends = len(trace)
            for i_dend in np.arange(0,nDends):
                Dtrace = trace[i_dend]
                if (i_dend == 0):
                    ax = plt.subplot(3, nPlots, mat_plots[2,i_plot]+1)
                ax.plot(data.taxis, Dtrace, color=cols[i_dend])
            plt.xlim(0,xmax)
            plt.ylim(-80,0)
    plt.show(block=False)
