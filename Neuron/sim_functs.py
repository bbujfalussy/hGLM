#----------------------------------------------------------------------------
# Functions and Classes
import struct
import sys
from numpy import loadtxt

def initOnsetSpikes():
    model.ncAMPAlist[0].event(data.st_onset * 1000)

def initSpikes():
    if (len(data.etimes)>0):
        for s in data.etimes:
            model.ncAMPAlist[int(s[0])].event(float(s[1]))
            if data.NMDA: model.ncNMDAlist[int(s[0])].event(float(s[1]))
    if data.GABA == True:
        if (len(data.itimes)>0):
            for s in data.itimes:
                model.ncGABAlist[int(s[0])].event(float(s[1]))

def storeSimOutput(v,vD):
        data.vdata.append(v)
        data.vDdata.append(vD)

def storeSimInputOutput(v,vD,Et,It):
        data.vdata.append(v)
        data.vDdata.append(vD)
        if (len(data.stim) == 0):
            nrep = 1
        else :
            nrep = max(data.stim[:,0]) + 1
            
        if (len(Et) > 0):
            Et = np.column_stack((nrep * np.ones(len(Et[:,0])), Et))
            if (len(data.stim) == 0):
                data.stim = Et
            else :
                data.stim = np.row_stack((data.stim, Et))
        if (len(It) > 0) :
            It = np.column_stack((-1 * nrep * np.ones(len(It[:,0])), It))
            if (len(data.stim) == 0):
                data.stim = It
            else :
                data.stim = np.row_stack((data.stim, It))

#-----------------------------------------------
# Synapse location functions
#-----------------------------------------------

# def gen1dendLocs(dend, nsyn, spread):
#     locs = []
#     isd = (spread[1]-spread[0])/(nsyn)
#     pos = np.arange(0, spread+isd/10., isd) + 0.5-spread/2
#     pos = np.arange(spread[0], spread[1], isd)[0:nsyn]
#     for p in pos:
#         locs.append([dend, p])
#     return locs

        

def genDendLocs(dends, nsyn, spread):
    locs = []
    n_dends = len(dends)
    if len(nsyn) != len(dends):
        nsyn = np.repeat(nsyn[0], len(dends))
    for i_dend in np.arange(0,n_dends):
        dend = dends[i_dend]
        nsyn_dend = nsyn[i_dend]

        isd = (spread[1]-spread[0])/(nsyn_dend)
        pos = np.arange(spread[0], spread[1], isd)[0:nsyn_dend] 
        if (len(pos) != nsyn_dend):
            print 'error: synapse number mismatch, stop simulation! dend:', i_dend, 'created=', len(pos), '!=', nsyn_dend
            sys.exit(1)
        for p in pos:
            locs.append([dend, p])

    return locs

def genDendSomaLocs(dends, nsyn, spread):
    locs = []
    n_dends = len(dends)
    nsyn_soma = nsyn[0]
    nsyn_dends = np.delete(nsyn, 0)
    for p in np.arange(0,nsyn_soma):
        locs.append([-1, 0.5])# -1 is not an index - it means that the synapse is at the soma
    for i_dend in np.arange(0,n_dends):
        dend = dends[i_dend]
        nsyn_dend = nsyn_dends[i_dend]
        isd = (spread[1]-spread[0])/(nsyn_dend)
        pos = np.arange(spread[0], spread[1], isd)[0:nsyn_dend] 
        if (len(pos) != nsyn_dend):
            print 'error: synapse number mismatch, stop simulation! dend:', i_dend, 'created=', len(pos), '!=', nsyn_dend
            sys.exit(1)
        for p in pos:
            locs.append([dend, p])
    return locs

def genSomaLocs(nsyn):
    locs = []
    if len(nsyn) > 1:
        nsyn = sum(nsyn)
    for p in np.arange(0,nsyn):
        locs.append([-1, 0.5])# -1 is not an index - it means that the synapse is at the soma
    return locs


def genRandomLocs(nsyn):
    locs = []
    nsyn = sum(nsyn)
    for s in np.arange(0,nsyn):
        dend = np.random.randint(low=0, high=len(model.dends))
        pos = np.random.uniform()
        locs.append([dend, pos])
    return locs

def genAllLocs(isd=1):   # intersynaptic distance in microns 
    locs = []
    dend_n = 0
    for dend in model.dends:
        distance = 0
        while distance<dend.L:
            locs.append([dend_n, distance/dend.L])
            distance = distance + isd
        dend_n = dend_n + 1
    return locs # list of [dend distance] pairs

def readDendLocs(data):
    locs = []
    if data.NMDA:
        data.Ntype='NMDA'
    else:
        data.Ntype='AMPA'
    if data.ACTIVE:
        if data.ACTIVEdend:
            data.act='Adend'
        else :
            data.act='Asoma'
    else:
        data.act='passive'

    outdir = data.simType + '/4B_passive'
    filename=data.Ntype+'_ApN1_'+data.locBias+'_passive_'+str(data.ind_Rates)
    bfname = './'+outdir+'/synlocs_'+filename+'.bin'
    eilocs = np.fromfile(bfname,dtype=float)
    eilocs = np.delete(eilocs,0)
    eilocs = np.reshape(eilocs, (2,192)).T
    Elocs = eilocs[0:160,:]
    Ilocs = eilocs[160:192,:]
    Ilocs[:,1] = Ilocs[:,1] - 1
    return Elocs, Ilocs



#-----------------------------------------------
# Input generation functions
#-----------------------------------------------

def genPoissonPulse(nsyn, rate, duration, onset):
    times = np.array([])
    print 'Poisson train for', nsyn, 'cells with', rate, 'Hz rate and', duration, 'ms duration'
    while times.shape[0]<2: # at least two spikes are required
        P =  br.OfflinePoissonGroup(nsyn, rate, duration * br.ms)
        times = np.array(P.spiketimes)    
    times[:,1] = times[:,1] * 1000 + onset
    return times


def genPoissonTrain(Ensyn, bgErate, fgErate, duration, omegaUP, omegaDown, Insyn=0, bgIrate = 0, fgIrate = 0, random_seed=1):
    ## rates are in Hz
    ## duration: total simulation time in s
    ## omegaUp: the rate of switching from 0 to 1, 1/s
    ## omegaDown: switching from Up to Down, 1/s
    ## times: two columns: cell id, spike times (ms)
    np.random.seed(random_seed)
    state = np.int(0)
    t_last = 0
    Etimes = np.array([])
    Ensyn = sum(Ensyn)
    Insyn = sum(Insyn)
    
    if (Insyn > 0) : Itimes = np.array([])
    while (Etimes.shape[0]<2): # at least two spikes are required
        Etimes = np.array([[0,0]])
        if (Insyn > 0) : Itimes = np.array([[0,0]])
        while (t_last < duration):
            if (state == 0):
                omega = omegaUP
            else:
                omega = omegaDown
            if (state == 0):
                Erate = bgErate
                if (Insyn > 0) : Irate = bgIrate
            else:
                Erate = fgErate
                if (Insyn > 0) : Irate = fgIrate
            t_next = t_last + np.random.exponential(1.0/omega)
            if (t_next > duration): t_next = duration
            dur_state = (t_next - t_last) * 1000
            PE =  br.OfflinePoissonGroup(Ensyn, Erate, dur_state * br.ms)
            Ettimes = np.array(PE.spiketimes)
            if (Ettimes.shape[0] > 0):
                Ettimes[:,1] = Ettimes[:,1] * 1000 + t_last * 1000
                Etimes = np.vstack((Etimes,Ettimes))
            if ((Insyn > 0)&(Irate > 0)) :
                PI =  br.OfflinePoissonGroup(Insyn, Irate, dur_state * br.ms)
                Ittimes = np.array(PI.spiketimes)
                if (Ittimes.shape[0] > 0):
                    Ittimes[:,1] = Ittimes[:,1] * 1000 + t_last * 1000
                    Itimes = np.vstack((Itimes,Ittimes))
            state = np.int(1-state)
            t_last = t_next
#            print t_last
            
        Etimes = np.delete(Etimes,0,0)
        if (Insyn > 0) :
            Itimes = np.delete(Itimes,0,0)
            times = [Etimes, Itimes]
        else :
            times = Etimes
    return times

def readOriTrain(duration, ir, rep, Ensyn, bgErate, fgErate, Insyn, bgIrate, fgIrate):

    ## eitimes = readPoissonTrain(duration, ir, rep, Ensyn, bgErate, fgErate, Insyn, bgIrate, fgIrate)
    ## rates are in Hz
    ## duration: total simulation time in s
    ## times: two columns: cell id, spike times (ms)
    Etimes = np.array([])
    Ensyn = sum(Ensyn)
    Insyn = sum(Insyn)
    duration = duration*1000
    
    fname = './allDend/orientations/Espikes_d'+str(duration)+'_r'+str(ir+1)+'_rep'+str(rep+1)+'_Ne'+str(Ensyn)+'_e'+str(bgErate)+'_E'+str(fgErate)+'.dat'
    Etimes = loadtxt(fname, comments="#", delimiter=" ", unpack=False)    
    n_sp = len(Etimes)

    
    if (Insyn > 0) :
        Itimes = np.array([])
        fname = './allDend/orientations/Ispikes_d'+str(duration)+'_r'+str(ir+1)+'_rep'+str(rep+1)+'_Ni'+str(Insyn)+'_i'+str(bgIrate)+'_I'+str(fgIrate)+'.dat'
        Itimes = loadtxt(fname, comments="#", delimiter=" ", unpack=False)    
        times = [Etimes, Itimes]
    else :
        times = Etimes

    return times

def readPoissonTrain(duration, ir, ori, Ensyn, bgErate, fgErate, Insyn, bgIrate, fgIrate):

    ## eitimes = readPoissonTrain(duration, ir, ori, Ensyn, bgErate, fgErate, Insyn, bgIrate, fgIrate)
    ## rates are in Hz
    ## duration: total simulation time in s
    ## times: two columns: cell id, spike times (ms)
    Etimes = np.array([])
    Ensyn = sum(Ensyn)
    Insyn = sum(Insyn)
    duration = duration*1000
    
    fname = './allDend/orientations/spE_d'+str(duration)+'_r'+str(ir+1)+'_o'+str(ori+1)+'_Ne'+str(Ensyn)+'_e'+str(bgErate)+'_E'+str(fgErate)+'.dat'
    Etimes = loadtxt(fname, comments="#", delimiter=" ", unpack=False)    
    n_sp = len(Etimes)

    
    if (Insyn > 0) :
        Itimes = np.array([])
        fname = './allDend/orientations/spI_d'+str(duration)+'_r'+str(ir+1)+'_o'+str(ori+1)+'_Ni'+str(Insyn)+'_i'+str(bgIrate)+'_I'+str(fgIrate)+'.dat'
        Itimes = loadtxt(fname, comments="#", delimiter=" ", unpack=False)    
        times = [Etimes, Itimes]
    else :
        times = Etimes

    return times


def genDSinput(nsyn, tInterval, onset, direction):
    times = np.zeros([nsyn, 2])
    if (direction=='OUT'):
        times[:,0] = np.arange(0, nsyn)
    else:
        times[:,0] = np.arange(nsyn-1, -1, -1)
    times[:,1] = np.arange(0, nsyn*tInterval, tInterval) + onset
    return times

def addBground(times1, times2):
    syns = np.concatenate([times1[:,0], times2[:,0]])
    times = np.concatenate([times1[:,1], times2[:,1]])
    inds = np.argsort(times)
    et = np.array([syns[inds], times[inds]])
    times = np.transpose(et)
    return times

def gen_bursts(Erate, minN, maxN, Tmax, nsyns):
    # generate a series of presynaptic bursts of minN < N < maxN presynaptic E cells
    # Erate is the average firing rate of the cells - the burst rate is calculated fomr 
    # from Erate = (minN + maxN)/2 * burst_rate / nsyns
    # as:
    burst_rate = round(Erate * nsyns * 2 / float(minN + maxN), 2)
    burst_times = genPoissonPulse(nsyn=1, rate=burst_rate, duration=Tmax * 1000, onset=0)
    nBursts = burst_times.shape[0]
    for iburst in range(nBursts):
        tBurst = burst_times[iburst, 1]
        nInput = int(np.random.randint(minN, maxN+1, 1)) 
        start_syn = np.random.randint(0, nsyns - nInput + 1, 1)
        active_syns = np.arange(start_syn, start_syn + nInput)
        t_spikes = tBurst + np.arange(0, nInput) / float(1)
        burst_i = np.transpose(np.array([active_syns, t_spikes]))
        if (iburst == 0):
            burst = burst_i
        else:
            burst = np.vstack((burst, burst_i))
        # if (iburst == 25):
        # print(burst_i)
        # print tBurst, nInput
        # print burst.shape 
        # print(data.etimes)
    print('number of bursts:', nBursts, ', number of spikes in bursts:', burst.shape[0])
    print('expected firing rate:', Erate, ', real firing rate:', round(burst.shape[0] / float(nsyns) / float(Tmax), 2))
    return(burst)

#--------------------------------------------------
# Simulation functions
#--------------------------------------------------

def sim_nsynPoiss(Ensyn, Erate, bRate=0, minN=5, maxN=20):
    # simulates the respons to synchron bursts embedded in a Poisson background input
    # simulates only excitatory cells
    # bRate is the rate of bursts, Erate is the rate for the Poisson background
    # minN and maxN are the number of synapses participating in a single burst

    Ensyn = sum(Ensyn)
    # print Ensyn, Erate, data.TSTOP * 1000
    data.etimes = genPoissonPulse(Ensyn, Erate, data.TSTOP * 1000, 0)
    data.itimes = []
    print 'background spikes:', data.etimes.shape[0]

    if (bRate>0): # bursts
        fetimes = genPoissonPulse(1, bRate, data.TSTOP * 1000, 0)
        nBursts = fetimes.shape[0]
        for iburst in np.arange(1,nBursts):
            tBurst = fetimes[iburst, 1]
            nInput = int(np.random.randint(minN, maxN+1, 1)) 
            synsInput = np.random.choice(Ensyn, nInput, replace=False)
            burst_i = np.transpose(np.array([synsInput, np.tile(tBurst, nInput)]))
            if (iburst == 0):
                burst = burst_i
            else:
                burst = np.vstack((burst, burst_i))
            # print tBurst, nInput
            # print burst.shape 
            # print(data.etimes)
    print 'total spikes:', data.etimes.shape[0]
    data.etimes = addBground(burst, data.etimes)

    # Run
    fih = lb.h.FInitializeHandler(1, initSpikes)
    taxis, v, vD = lb.simulate(model, t_stop=data.TSTOP * 1000,
                                     NMDA=data.NMDA, recDend=data.recordDend,
                                     i_recDend=data.locDendRec)
    Et, It = data.etimes, data.itimes
    return taxis, v, vD, Et, It


def sim_oneRandomInput(Ensyn, Insyn, Erate, Irate, bErate=0, bIrate=5):
    # simulates response to a single presynaptic up state
    # Ensyn and Insyn are lists - but may contain only one element
    Ensyn = sum(Ensyn)
    Insyn = sum(Insyn)
    print Ensyn, Erate, data.st_duration * 1000, data.st_onset * 1000
    data.etimes = genPoissonPulse(Ensyn, Erate, data.st_duration * 1000, data.st_onset * 1000)
#        print(data.etimes.shape[0])
    if (bErate>0):
        data.betimes = genPoissonPulse(Ensyn, bErate, data.TSTOP * 1000, 0)
#            print(data.betimes.shape[0])
        data.etimes = addBground(data.betimes, data.etimes)
    if (Irate>0):
        data.itimes = genPoissonPulse(Insyn, Irate, data.st_duration * 1000, data.st_onset * 1000)
    if (bIrate>0):
        data.bitimes = genPoissonPulse(Insyn, bIrate, data.TSTOP * 1000, 0)
        if (Irate>0):
            data.itimes = addBground(data.bitimes, data.itimes)
        else :
            data.itimes = data.bitimes
    if ((Irate==0) & (bIrate==0)): data.GABA=False 


    # Run
    fih = lb.h.FInitializeHandler(1, initSpikes)
    taxis, v, vD = lb.simulate(model, t_stop=data.TSTOP * 1000,
                                     NMDA=data.NMDA, recDend=data.recordDend,
                                     i_recDend=data.locDendRec)
    Et, It = data.etimes, data.itimes
    return taxis, v, vD, Et, It



def sim_balanceInput(Ensyn, Insyn, fgErate, fgIrate, bgErate, bgIrate, duration, omega, rand_seed=1):
    # generates balanced E and I inputs for a number of presynaptic assemblies
    # Ensyn is a list containing vectors with the synapses in a given assembly
    nden = len(Ensyn)
    Ensyn_sig = [Ensyn[0]]
    Insyn = [sum(Insyn)]
    bgIr = bgIrate / float(nden)
    fgIr = fgIrate / float(nden)
    eitimes = genPoissonTrain(Ensyn_sig, bgErate, fgErate, duration, omega[0], omega[1], Insyn, bgIr, fgIr, rand_seed)
    etimes = eitimes[0]
    itimes = eitimes[1]
    print 'number of dendrites stimulated: ', nden
    print len(etimes[:,0]), 'E spikes generated for dendrite 1'
    if (nden>1):
        for iden in np.arange(1,nden):
            Ensyn_sig = [Ensyn[iden]]
            beitimes = genPoissonTrain(Ensyn_sig, bgErate, fgErate, duration, omega[0], omega[1], Insyn, bgIr, fgIr, rand_seed+iden)
            betimes = beitimes[0]
            bitimes = beitimes[1]
            betimes[:,0] = betimes[:,0] + sum(Ensyn[0:iden])
            etimes = addBground(betimes, etimes)
            itimes = addBground(bitimes, itimes)
            print len(betimes[:,0]), 'E spikes generated for dendrite', iden + 1
    data.etimes = etimes
    data.itimes = itimes

    # Run
    fih = lb.h.FInitializeHandler(1, initSpikes)
    taxis, v, vD = lb.simulate(model, t_stop=data.TSTOP * 1000,
                                     NMDA=data.NMDA, recDend=data.recordDend,
                                     i_recDend=data.locDendRec)

    return taxis, v, vD, etimes, itimes


def sim_oriInput(Ensyn, Insyn, fgErate, fgIrate, bgErate, bgIrate, duration, ori, ir):
# simulates response to orientation selective inputs read from file
    nden = len(Ensyn)
    Ensyn = [sum(Ensyn)]
    Insyn = [sum(Insyn)]
    bgIr = bgIrate / float(nden)
    fgIr = fgIrate / float(nden)

    eitimes = readPoissonTrain(duration, ir, ori, Ensyn, bgErate, fgErate, Insyn, bgIrate, fgIrate)
    etimes = eitimes[0]
    itimes = eitimes[1]
    print 'all of the dendrites are stimulated'
    print len(etimes[:,0]), 'E and ',  len(itimes[:,0]),  'I spikes read from file'
    data.etimes = etimes
    data.itimes = itimes
 
    # Run
    fih = lb.h.FInitializeHandler(1, initSpikes)
    taxis, v, vD = lb.simulate(model, t_stop=data.TSTOP * 1000,
                                     NMDA=data.NMDA, recDend=data.recordDend,
                                     i_recDend=data.locDendRec)

    return taxis, v, vD, etimes, itimes


def sim_mixedoriInput(Ensyn, Insyn, fgErate, fgIrate, bgErate, bgIrate, duration, rep, ir, seed=12):
# simulates response to gratings with mixed orientations. Input spikes are read from file.
    nden = len(Ensyn)
    Ensyn = [sum(Ensyn)]
    Insyn = [sum(Insyn)]
    random_seed = seed

    eitimes = readOriTrain(data.TSTOP, data.ind_Rates, rep, Ensyn, bgErate, fgErate, Insyn, bgIrate, fgIrate)
    etimes = eitimes[0]
    itimes = eitimes[1]
    print 'all of the dendrites are stimulated'
    print len(etimes[:,0]), 'E and ',  len(itimes[:,0]),  'I spikes read from file'

    if (data.inburstrate > 0):
        np.random.seed(random_seed)
        etimes_burst = gen_bursts(Erate=data.inburstrate, minN=5, maxN=30, Tmax=data.TSTOP, nsyns=sum(Ensyn))
        etimes = addBground(etimes, etimes_burst)
    data.etimes = etimes
    data.itimes = itimes
 
    # Run
    fih = lb.h.FInitializeHandler(1, initSpikes)
    taxis, v, vD = lb.simulate(model, t_stop=data.TSTOP * 1000,
                                     NMDA=data.NMDA, recDend=data.recordDend,
                                     i_recDend=data.locDendRec)

    return taxis, v, vD, etimes, itimes


#--------------------------------------------------
# Iteration functions
#--------------------------------------------------

def SIM_balanceIteration(EIrates):
    fgErate, fgIrate, bgErate, bgIrate = EIrates
    for iter in np.arange(data.nIter):
        r_seed = data.r_seed * 100 + 4 * iter
        print 'Running E bg. rate', bgErate, 'I bg. rate', bgIrate, 'E rate', fgErate, 'I rate', fgIrate, 'iteration:', iter

        # sim_balanceInput(Ensyn=data.Ensyn, Insyn=data.Insyn, 
        #                                               fgErate=fgErate, fgIrate=fgIrate,
        #                                               bgErate=bgErate, bgIrate=bgIrate, 
        #                                               duration=data.TSTOP, omega=data.omega, rand_seed=r_seed)
        data.taxis, v, vD, Et, It = sim_balanceInput(Ensyn=data.Ensyn, Insyn=data.Insyn, 
                                                      fgErate=fgErate, fgIrate=fgIrate,
                                                      bgErate=bgErate, bgIrate=bgIrate, 
                                                      duration=data.TSTOP, omega=data.omega, rand_seed=r_seed)
        storeSimInputOutput(v,vD,Et,It)


def SIM_rateIteration(ERate, IRate):
    nIter = len(ERate)
    for iter in np.arange(nIter):
        print 'Running E rate', ERate[iter], 'iteration:', iter
        np.random.seed(data.r_seed * 100 + iter)
        data.taxis, v, vD, Et, It = sim_oneRandomInput(Ensyn=data.Ensyn, Insyn=data.Insyn, 
                                                        Erate=ERate[iter], Irate=IRate[iter],
                                                        bErate=data.EbGroundRate,
                                                        bIrate=data.IbGroundRate)
        storeSimInputOutput(v,vD,Et,It)        


def SIM_oriIteration(EIrates):
    fgErate, fgIrate, bgErate, bgIrate = EIrates
    for iter in np.arange(data.nIter):
        r_seed = data.r_seed * 100 + 4 * iter
        print 'Running E bg. rate', bgErate, 'I bg. rate', bgIrate, 'E rate', fgErate, 'I rate', fgIrate, 'iteration:', iter

        data.taxis, v, vD, Et, It = sim_oriInput(Ensyn=data.Ensyn, Insyn=data.Insyn, 
                                                      fgErate=fgErate, fgIrate=fgIrate,
                                                      bgErate=bgErate, bgIrate=bgIrate, 
                                                      duration=data.TSTOP, ori=iter, ir=data.ind_Rates)
        storeSimInputOutput(v,vD,Et,It)



def SIM_mixedoriIteration(EIrates):
    fgErate, fgIrate, bgErate, bgIrate = EIrates
    for iter in np.arange(data.nIter):
        r_seed = data.r_seed * 100 + 4 * iter
        print 'Running E bg. rate', bgErate, 'I bg. rate', bgIrate, 'E rate', fgErate, 'I rate', fgIrate, 'iteration:', iter

        data.taxis, v, vD, Et, It = sim_mixedoriInput(Ensyn=data.Ensyn, Insyn=data.Insyn, 
                                                      fgErate=fgErate, fgIrate=fgIrate,
                                                      bgErate=bgErate, bgIrate=bgIrate, 
                                                      duration=data.TSTOP, rep=iter, ir=data.ind_Rates, seed=iter + 100)
        storeSimInputOutput(v,vD,Et,It)


        
def SIM_nsynIteration(tInterval, bGround=False):
#   ---------------------------------------
#   first each synapse is activated alone
    It = []
    for nsyn in np.arange(1, data.Ensyn[0]+1):
        print 'activating synapse ', nsyn, ' alone'
        data.etimes = np.array([[nsyn-1, data.st_onset * 1000]])
        fih = lb.h.FInitializeHandler(1, initSpikes)
        taxis, v, vD = lb.simulate(model, t_stop=data.TSTOP * 1000,
                                         NMDA=data.NMDA, recDend=data.recordDend,
                                         i_recDend=data.locDendRec)
        Et = data.etimes
        storeSimInputOutput(v, vD, Et, It)

#   ---------------------------------------
#   next, synapses are activated together
    for nsyn in np.arange(1, data.Ensyn[0]+1):
        print 'activating ', nsyn, ' synapses together'
        data.etimes = genDSinput(nsyn, tInterval, data.st_onset * 1000,'OUT')
        # if bGround:
        #     data.etimes = addBground(data.bEnsyn, data.Ensyn,
        #                              data.EbGroundRate, data.etimes)
        #     data.itimes = addBground(data.bInsyn, data.Insyn,
        #                              data.IbGroundRate, [0,0])
        fih = lb.h.FInitializeHandler(1, initSpikes)
        taxis, v, vD = lb.simulate(model, t_stop=data.TSTOP * 1000,
                                         NMDA=data.NMDA, recDend=data.recordDend,
                                         i_recDend=data.locDendRec)
        Et = data.etimes
        storeSimInputOutput(v, vD, Et, It)
    data.taxis = taxis


def SIM_burstIteration(burstRate, maxN):
    Ensyn = data.Ensyn
    Erate = data.Erate

    for iter in np.arange(data.nIter):
        r_seed = data.r_seed * 100 + 4 * iter
        np.random.seed(r_seed)
        print 'Running E bg. rate', Erate, 'burst rate', burstRate, 'iteration:', iter

        # sim_balanceInput(Ensyn=data.Ensyn, Insyn=data.Insyn, 
        #                                               fgErate=fgErate, fgIrate=fgIrate,
        #                                               bgErate=bgErate, bgIrate=bgIrate, 
        #                                               duration=data.TSTOP, omega=data.omega, rand_seed=r_seed)
        data.taxis, v, vD, Et, It = sim_nsynPoiss(Ensyn, Erate, bRate=burstRate, minN=5, maxN=maxN)
        storeSimInputOutput(v,vD,Et,It)



def SIM_currentSteps(iRange, bGround=False):
    Et, It = [], []
    for step in iRange:
        if bGround:
            if (data.EbGroundRate>0):
                data.etimes = genPoissonPulse(data.Ensyn, data.EbGroundRate, data.TSTOP * 1000, 0)
                Et = data.etimes
            if (data.IbGroundRate>0):
                data.itimes = genPoissonPulse(data.Insyn, data.IbGroundRate, data.TSTOP * 1000, 0)
                It = data.itimes
                fih = lb.h.FInitializeHandler(1, initSpikes)
        print 'Running input step current', step
        model.stim.amp = step
        taxis, v, vD = lb.simulate(model, t_stop=data.TSTOP * 1000, NMDA=data.NMDA, recDend=data.recordDend, i_recDend=data.locDendRec)

        storeSimInputOutput(v, vD, Et, It)
    data.taxis = taxis


def SIM_minisDistribution(data):
#   ---------------------------------------
#   first each excitatory synapse is activated alone
    It, Et, data.etimes, data.itimes = [], [], [], []
    Enum = sum(data.Ensyn)
    for nsyn in np.arange(0, Enum):
    # for nsyn in np.arange(0, 10):
        print 'activating excitatory synapse ', nsyn, ' alone'
        data.etimes = np.array([[nsyn, data.st_onset * 1000]])
        fih = lb.h.FInitializeHandler(1, initSpikes)
        taxis, v, vD = lb.simulate(model, t_stop=data.TSTOP * 1000,
                                         NMDA=data.NMDA, recDend=data.recordDend,
                                         i_recDend=data.locDendRec)
        Et = data.etimes
        storeSimInputOutput(v, vD, Et, It)

#   ---------------------------------------
#   next, each inhibitory synapse is activated alone
    It, Et, data.etimes, data.itimes = [], [], [], []
    Inum = sum(data.Insyn)
    for nsyn in np.arange(0, Inum):
        print 'activating inhibitory synapse ', nsyn, ' alone'
        data.itimes = np.array([[nsyn, data.st_onset * 1000]])
        fih = lb.h.FInitializeHandler(1, initSpikes)
        taxis, v, vD = lb.simulate(model, t_stop=data.TSTOP * 1000,
                                         NMDA=data.NMDA, recDend=data.recordDend,
                                         i_recDend=data.locDendRec)
        It = data.itimes
        storeSimInputOutput(v, vD, Et, It)

    data.taxis = taxis
