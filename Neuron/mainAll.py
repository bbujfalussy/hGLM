import numpy as np
import struct
import matplotlib.pyplot as plt
import time

import brian as br

import libcell as lb
import saveClass as sc

# this does not creates a new module, just executes the content of the file
# as if it was written in the same file!
execfile('sim_functs.py') 


#----------------------------------------------------------------------------
# Data saving object; # Data storage lists
data = sc.emptyObject()
data.vdata, data.vDdata, data.stim = [], [], []

#----------------------------------------------------------------------------
# Simulation CONTROL
# PAR1 = 0
# PAR2 = 2
data.model = 'L23'
data.simType = 'allDend' # allDend
data.stimType = 'mixedori' # allDend, orientations, mixedori, Poisson, minis
data.actType = 'Adend' # passive, aSoma, Adend

### controls the up-state transitions for allDend
iomega = 2
data.omegai = 1 # 1

data.ind_Rates = 0 #PAR1 # 0-5
iclust = 2 #PAR2 # 0-5; 1, 4, 13, 100, 200, 629
data.iAmp = 2 # 0-2 - the amplitude of the dendritically injected current in the orientation tuning experiments

data.r_seed = 1
noSynDend = []

data.inburstrate = 0

#---------------------------------------------------------------------------
# simtype sets the location of the synapses and the stimulus type
data.SHOWTRACES = True
data.SAVE = True
data.SHOWSYNS = True

### number of iterations for allDend; number of orientations for ori
data.nIter = 2 # 10 for allDend; 16 for ori and 10 for mixedori
### time parameters
data.st_onset = 0.11 # in seconds - only for minis
data.st_duration = 2.5 # doesn't matter
data.TSTOP = 2 # 41 for allDend, 18 for orientations and 48 for mixedori, 0.3 for minis

data.ACTIVE = True
data.ACTIVEaxonSoma = True
data.ACTIVEdend = False
data.ACTIVEdendNa = False
data.ACTIVEdendCa = False
data.ACTIVEhotSpot = False
if (data.actType == 'passive'):
   data.ACTIVE = False
   data.ACTIVEaxonSoma = False
if (data.actType == 'Adend'):
    data.ACTIVEdend = True
    data.ACTIVEdendNa = True
    data.ACTIVEdendCa = True
    data.ACTIVEhotSpot = True
    data.ACTIVEaxonSoma = False


data.SPINES = False

#----------------------------------------------------------------------------
# synapse parameters
#----------------------------------------------------------------------------
data.SYN = True
data.GABA = True    
data.NMDA = True
data.NMDAkinetic = False
data.ApN = 0.5 # AMPA per NMDA ratio

data.Egmax = 0.5 # nS - NMDA max conductance
data.Igmax = 1
data.Irev = -80
data.locBias = 'none'

execfile('init_params.py')

#----------------------------------------------------------------------------
# Simulation general parameters
data.dt = 0.2
lb.h.dt = data.dt
lb.h.steps_per_ms = 1.0/lb.h.dt
data.recordDend = True

#----------------------------------------------------------------------------
# Create neuron and add mechanisms
# if data.model == 'BS': model = lb.BS()

if data.model == 'L23': model = lb.L23()
if data.SPINES: lb.addSpines(model)
if data.ACTIVE: lb.init_active(model, axon=data.ACTIVEaxonSoma,
                             soma=data.ACTIVEaxonSoma, dend=data.ACTIVEdend,
                             dendNa=data.ACTIVEdendNa, dendCa=data.ACTIVEdendCa)
if data.ACTIVE:
    if data.ACTIVEhotSpot: lb.hotSpot(model)

if (data.ICLAMP):
    if data.iclampLoc[0]=='soma':
        lb.add_somaStim(model, data.iclampLoc[1], onset=data.iclampOnset,
                        dur=data.iclampDur, amp=data.iclampAmp)
    if data.iclampLoc[0]=='dend':
        lb.add_dendStim(model, data.iclampLoc[1], data.iclampLoc[2],
                 onset=data.iclampOnset, dur=data.iclampDur, amp=data.iclampAmp)


#----------------------------------------------------------------------------
# Generate synapse locations
#----------------------------------------------------------------------------

# Elocs: synapse locations - list, with elements [dendrite number, synapse location]
data.Elocs = genAllLocs(10) # distance between synapses
print 'number of excitatory synapses:', len(data.Elocs)

# inhibitory synapses - some are on the soma, others on the dendrites
#data.nIter = len(data.rateE) 
if (data.stimType == "allDend") :
 	# Nsomas = [420, 100, 100, 100, 420, 420, 420] # looks that this used to be 420 420 420 100 100 100
 	Nsomas = [420, 420, 420, 100, 100, 100] # looks that this used to be 420 420 420 100 100 100
 	nsyn_soma = Nsomas[data.ind_Rates]
else :
	nsyn_soma = 1	
if data.GABA:
    Isomalocs = []
    for p in np.arange(0,nsyn_soma):
        Isomalocs.append([-1, 0.5])# -1 is not an index - it means that the synapse is at the soma
    Idendlocs = genAllLocs(100)
    data.Ilocs = Isomalocs + Idendlocs
 #   data.Ilocs = Isomalocs
    print 'number of somatic synapses:', nsyn_soma
    print 'number of inhibitory synapses:', len(data.Ilocs)


#########################################################################
## number of synapses in a given cluster
## Ensyn: list with the number of E synapses in a particular cluster
## must sum to len(data.Elocs)
maxClusts = [100.,100.,100., 100., 4., 1.]
maxClust = maxClusts[iclust] # maximum number of synapses / cluster; 1-4-100
nsyn_tree = np.empty((0, 2))
elocs = np.asarray(data.Elocs)
for ibranch in np.arange(max(elocs[:,0])+1):
	nsyn_branch = float(sum(elocs[:,0] == ibranch))
	nclust_branch = np.ceil(nsyn_branch / maxClust)
	syn_branch = np.ones((nclust_branch, 2)) * maxClust
	if ((nsyn_branch % maxClust) > 0) :
		syn_branch[nclust_branch-1, 1] = nsyn_branch % maxClust
	syn_branch[:,0] = ibranch
	nsyn_tree = np.vstack((nsyn_tree, syn_branch))

data.Ensyn = list(nsyn_tree[:,1])
data.Insyn = [nsyn_soma, 11,11, 9,6,8,5,8, 12,11,13,6, 11,8] # use this if stimType == the orientation or allDend

# data.Ensyn = [48,58, 52,34,45,38,44, 68,50,62,31, 60,39] # 629 - use this if stimType == orientations or mixedori

## Insyn: number of I synapses in a given cluster
## must sum to len(data.Ilocs)
## otherwise its structure doesn't matter, since I synapses are not clustered

# data.Insyn = [120, 119] # 239 - used to use this for allDend

if (iclust == 0): 
	data.Ensyn = [629] # 629 - one cluster
	data.Insyn = [nsyn_soma, 119] # 119

if (iclust == 1): 
	data.Ensyn = [106, 213, 211, 99] # 629 - four clusters
	data.Insyn = [nsyn_soma, 22, 36, 42, 19] # 

if (iclust == 2): 
    data.Ensyn = [48,58, 52,34,45,38,44, 68,50,62,31, 60,39]
    data.Insyn = [nsyn_soma, 11,11, 9,6,8,5,8, 12,11,13,6, 11,8] # 239

#nClust=len(data.Ensyn)
#import cell_save as cs
#cs.save_syn(data, outdir="AllDend/passive", clust=nClust)

#######################################################################################
## which dendrites to use - if not All!
data.locDend = [11, 24, 70, 95] # 1dend: 11; 4dend: 11 - g2, 24 - g3, 70 - g10, 95 - g12
# best for stimulation - terminal dendrites with diameter  < 1 um
# 11 - dend1_1212 L = 121um, d = 0.6 um
# 24 - dend2_1122 L = 129um, d = 0.6 um
# 70 - dend3_12122221 L = 231um, d = 0.75 um
# 95 - dend4_1121 L = 137um, d = 0.5 um

## which dendrites to record from
data.locDendRec = [11, 24, 70, 95] 
# best for recording - non-terminal apical branches, 100 um from soma, diam > 1 um
# 38 - dend2_1212[0.5] L = 140 um, d=1.1 um, dist ~ 157 um
# 47 - dend2_1222[0.5] L=150, d=1.1, dist ~ 157 um
# 69 - dend3_1212222[0.5] L=47, d=1.4, dist ~ 135 um
# 70 - dend3_12122221[0.5] L = 231um, d = 0.75 um
# mainAll.model.dends[11].name()
# mainAll.lb.h.topology()

# np.random.shuffle(data.Elocs)

# Insert synapses
lb.add_AMPAsyns(model, locs=data.Elocs, gmax=data.ApN * data.Egmax, NoSynDends=noSynDend)
if (data.NMDAkinetic):
    lb.add_NMDAkin_syns(model, locs=data.Elocs, gmax=data.Egmax, NoSynDends=noSynDend)
else:
	if (data.NMDA):
	    lb.add_NMDAsyns(model, locs=data.Elocs, gmax=data.Egmax, NoSynDends=noSynDend)
	else:
	    lb.add_NMDAsyns(model, locs=data.Elocs, gmax=0, NoSynDends=noSynDend)
lb.add_GABAsyns(model, locs=data.Ilocs, gmax=data.Igmax, rev=data.Irev, NoSynDends=noSynDend)

# ----------------------------------------------------------------------------
# Run simulation
# ----------------------------------------------------------------------------
print 'simulation type: ', data.simType
if (data.ACTIVE == True):
    if (data.ACTIVEdend == True) :
        print 'dendrites are active, all branches are stimulated'
    else :
        print 'the soma is active, all branches are stimulated'
else :
    print 'the neuron is passive, all branches are stimulated'
if (data.stimType == "allDend") :
    SIM_balanceIteration(data.EIrates)
elif (data.stimType == "orientations"):
    SIM_oriIteration(data.EIrates)
elif (data.stimType == "mixedori"):
    SIM_mixedoriIteration(data.EIrates)
elif (data.stimType == "minis"):
	SIM_minisDistribution(data)


#----------------------------------------------------------------------------
# show the traces

if data.SHOWTRACES:
    print 'plotting traces ...'
    import cell_traces as ct
    if (data.stimType == 'minis'):
        ct.plotTraces(data, data.TSTOP * 1000)
    else :
        ct.plotResp(data, data.TSTOP * 1000)


#----------------------------------------------------------------------------
# Save data - number of clusters are encoded in the locBias parameter
data.locBias = str(len(data.Ensyn))
if (len(data.Ensyn) < 10 ):
   data.locBias = '0' + data.locBias

if data.SAVE:
    print 'saving data ...'
    if (data.stimType == "allDend") :
       outdir = data.simType + '/' + data.actType
    else:
       outdir = data.simType + '/' + data.stimType

    import cell_save as cs
    cs.save_sim(data, out_binary=True, out_vdend=True, out_pickle=False, 
    			outdir=outdir, rate=data.ind_Rates, omega=iomega, iAMP=data.iAmp)

#----------------------------------------------------------------------------
# show the synapses

if data.SHOWSYNS:
    print 'plotting the cell with the synapses ...'
    import cell_draw as cd
    if (data.simType=='nsyn') :
        cd.plot_syns(data, model, False)
    else : 
        cd.plot_syns(data, model, True)

# #----------------------------------------------------------------------------
# # inter-synapse matrix

# # print 'calculating the inter-synapse distance matrix'

# # DSyn = lb.synDist(model, data.Elocs)
# # imgplot = plt.imshow(DSyn)
# # plt.show(block=False)

# # print 'calculating the inter-synapse distance matrix'
# # DSynE = lb.synDistFromSoma(model, data.Elocs)
# # DSynI = lb.synDistFromSoma(model, data.Ilocs)

# # bfname='EdistFromSoma.bin'
# # binfile = file(bfname, 'wb')
# # # and write out two integers with the row and column dimension
# # header = struct.pack('2I', DSynE.shape[0], DSynE.shape[1])
# # binfile.write(header)
# # # then loop over columns and write each
# # for i in range(DSynE.shape[1]):
# #     ddata = struct.pack('%id' % DSynE.shape[0], *DSynE[:,i])
# #     binfile.write(ddata)
# # binfile.close()


# # bfname='IdistFromSoma.bin'
# # binfile = file(bfname, 'wb')
# # # and write out two integers with the row and column dimension
# # header = struct.pack('2I', DSynI.shape[0], DSynI.shape[1])
# # binfile.write(header)
# # # then loop over columns and write each
# # for i in range(DSynI.shape[1]):
# #     ddata = struct.pack('%id' % DSynI.shape[0], *DSynI[:,i])
# #     binfile.write(ddata)
# # binfile.close()
