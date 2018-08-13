# ----------------------------------------------------------
# Library of cell classes and functions
#
# Tiago Branco, MRC Laboratory of Molecular Biology, 2013
# email: tbranco@mrc-lmb.cam.ac.uk
# ----------------------------------------------------------

import numpy as np
import neuron

from neuron import h
from neuron import load_mechanisms
from neuron import gui

#load_mechanisms('/directory_where_mod_files_have_been_compiled')
h('objref nil')

# ----------------------------------------------------------
# MODELS
class L23(object):

    # Cell morphology is from cat, all lengths and diameters
    # are scaled to 70% to approximate it to mouse values

    def __init__(self):
        h('xopen("./L23.hoc")')
        props(self)
        self._geom()
        self._topol()
        self._changeLength()
        self._biophys()

    def _geom(self):
        self.axon = h.Section()
        # self.axon.L = 1
        self.axon.L = 300
        self.axon.diam = 1

    def _topol(self):            
        self.soma = h.soma
        self.dends = []
        for sec in h.allsec():
            self.dends.append(sec)
            sec.nseg = 7
        self.dends.pop()   # Remove soma from the list
        self.dends.pop()   # and the Axon
        for sec in self.dends:
            sec.diam = sec.diam * self.rescale
        self.axon.connect(self.soma,1,0)
    
    def _biophys(self):
        for sec in h.allsec():
            sec.cm = self.CM
            sec.insert('pas')
            sec.e_pas = self.E_PAS
            sec.g_pas = 1.0/self.RM
            sec.Ra = self.RA

    def _changeLength(self):
        for sec in h.allsec():
            sec.L = sec.L * self.rescale


# ----------------------------------------------------------
# INSTRUMENTATION FUNCTIONS
def props(model):

    # morphology - rescale factor
    model.rescale = 0.7 # 0.7
    
    # Passive properties
    model.CM = 1.0
    model.RM = 7000.0
    model.RA = 100.0 
    model.E_PAS = -75
    model.CELSIUS = 35

    # Active properties
    model.Ek = -90
    model.Ena = 60
    model.Eca = 140
    
    ## unit:  pS/um2 ; in comments: the values reported in Smith et al, 2013
    model.gna_axon = 1000 #  pS/um2 - missing
    model.gkv_axon = 100 #  pS/um2 - missing
    
    model.gna_soma = 1000 #  pS/um2 - 100 mS/cm2
    model.gkv_soma = 100 # pS/um2 - 10 mS/cm2
    model.gkm_soma = 2.2 #  pS/um2 
    model.gkca_soma = 3 #  pS/um2
    model.gca_soma = 0.5 # pS/um2 - 0.05 mS/cm2
    model.git_soma = 0.0003 #  S/cm2 - 0.0003 mS/cm2
    
    model.gna_dend = 80 #  pS/um2 - 60 mS/cm2
    model.gna_dend_hotSpot = 600 # no date
    model.gkv_dend = 3 # 0.3  pS/um2 - mS/cm2
    model.gkm_dend = 1 # pS/um2
    model.gkca_dend = 3 # pS/um2
    model.gca_dend = 0.5 # pS/um2 - 0.05 mS/cm2
    model.git_dend = 0.00015 # S/cm2 - 0.00015 ms/cm2
    model.gh_dend = 0

def init_active(model, axon=False, soma=False, dend=True, dendNa=False,
                dendCa=False):
    if axon:
        model.axon.insert('na'); model.axon.gbar_na = model.gna_axon
        model.axon.insert('kv'); model.axon.gbar_kv = model.gkv_axon
        model.axon.ena = model.Ena
        model.axon.ek = model.Ek
        print 'active conductances added in the axon'
        
    if soma:
        model.soma.insert('na'); model.soma.gbar_na = model.gna_soma
        model.soma.insert('kv'); model.soma.gbar_kv = model.gkv_soma
        model.soma.insert('km'); model.soma.gbar_km = model.gkm_soma
        model.soma.insert('kca'); model.soma.gbar_kca = model.gkca_soma
        model.soma.insert('ca'); model.soma.gbar_ca = model.gca_soma
        model.soma.insert('it'); model.soma.gbar_it = model.git_soma
        # model.soma.insert('cad');
        model.soma.ena = model.Ena
        model.soma.ek = model.Ek
        model.soma.eca = model.Eca
        print 'somatic active conductances enabled'
        
    if dend:
        for d in model.dends:
            d.insert('na'); d.gbar_na = model.gna_dend*dendNa
            d.insert('kv'); d.gbar_kv = model.gkv_dend
            d.insert('km'); d.gbar_km = model.gkm_dend
            d.insert('kca'); d.gbar_kca = model.gkca_dend
            d.insert('ca'); d.gbar_ca = model.gca_dend*dendCa
            d.insert('it'); d.gbar_it = model.git_dend*dendCa
            # d.insert('cad')
            d.ena = model.Ena
            d.ek = model.Ek
            d.eca = model.Eca
        print 'active dendrites enabled'

def add_somaStim(model, p=0.5, onset=20, dur=1, amp=0):
    model.stim = h.IClamp(model.soma(p))
    model.stim.delay = onset
    model.stim.dur = dur
    model.stim.amp = amp    # nA
    
def add_dendStim(model, p=0.5, dend=10, onset=20, dur=1, amp=0):
    model.stim = h.IClamp(model.dends[dend](p))
    model.stim.delay = onset
    model.stim.dur = dur
    model.stim.amp = amp    # nA


def synDist(model,locs):
    nsyn = len(locs)
    DSyn = np.zeros([nsyn, nsyn])
    fromSyn = 0
    for loc in locs:
        fromDend = loc[0]
        fromX = loc[1]
        fromSection = model.dends[fromDend]
        h.distance(0, fromX, sec=fromSection)
        toSyn = 0
        for toLoc in locs:
            toDend = toLoc[0]
            toX = toLoc[1]
            toSection = model.dends[toDend]
            x = h.distance(toX, sec=toSection)
            DSyn[toSyn, fromSyn] = x
            toSyn = toSyn + 1
        fromSyn = fromSyn + 1
    return DSyn


def synDistFromSoma(model,locs):
    nsyn = len(locs)
    DSyn = np.zeros([nsyn,1])
    fromSyn = 0
    h.distance(0, 0.5, sec=model.soma)
    toSyn = 0
    for toLoc in locs:
        toDend = toLoc[0]
        toX = toLoc[1]
        toSection = model.dends[toDend]
        if (toDend < 0):
            x = h.distance(toX, sec=model.soma)
        else:
            x = h.distance(toX, sec=toSection)
        DSyn[toSyn,0] = x
        toSyn = toSyn + 1
    return DSyn


def add_AMPAsyns(model, locs=[[0, 0.5]], gmax=0.5, tau1=0.1, tau2=2, NoSynDends=[]):
    model.AMPAlist = []
    model.ncAMPAlist = []
    gmax = gmax/1000.   # Set in nS and convert to muS
    for loc in locs:
        locInd = int(loc[0])
        if (locInd == -1):
            synloc = model.soma
        else:
            synloc = model.dends[int(loc[0])]
        AMPA = h.Exp2Syn(float(loc[1]), sec=synloc)
        AMPA.tau1 = tau1
        AMPA.tau2 = tau2
        if (int(loc[0]) in NoSynDends):
            gg = 0 # same input into a single branch
        else:
            gg = gmax
        NC = h.NetCon(h.nil, AMPA, 0, 0, gg) # NetCon(source, target, threshold, delay, weight)
        model.AMPAlist.append(AMPA)
        model.ncAMPAlist.append(NC)
    print 'AMPA synapses added'

def add_NMDAsyns(model, locs=[[0, 0.5]], gmax=0.5, tau1=3, tau2=40, NoSynDends=[]):
    model.NMDAlist = []
    model.ncNMDAlist = []
    gmax = gmax/1000.   # Set in nS and convert to muS
    for loc in locs:
        locInd = int(loc[0])
        if (locInd == -1):
            synloc = model.soma
        else:
            synloc = model.dends[int(loc[0])]
        NMDA = h.Exp2SynNMDA(float(loc[1]), sec=synloc) 
        NMDA.tau1 = tau1
        NMDA.tau2 = tau2
#        NMDA.e = 0
        if (int(loc[0]) in NoSynDends):
            gg = 0 # same input into a single branch
        else:
            gg = gmax
        NC = h.NetCon(h.nil, NMDA, 0, 0, gg)
        x = float(loc[1])
        model.NMDAlist.append(NMDA)
        model.ncNMDAlist.append(NC)   
    print 'dExp NMDA synapses generated'

def add_NMDAkin_syns(model, locs=[[0, 0.5]], gmax=0.5, tau1=3, tau2=40, NoSynDends=[]):
    model.NMDAlist = []
    model.ncNMDAlist = []
    gmax = gmax * 1000 * 5.4  # scaling gmax to match the 2exp synapse - see desens_test/test_NMDA.py
    for loc in locs:
        locInd = int(loc[0])
        if (locInd == -1):
            synloc = model.soma
        else:
            synloc = model.dends[int(loc[0])]
        NMDA = h.NMDA5d2nc(float(loc[1]), sec=synloc) 
        if (int(loc[0]) in NoSynDends):
            gg = 0 # same input into a single branch
        else:
            gg = gmax
        NMDA.gmax = gg
        NC = h.NetCon(h.nil, NMDA, 0, 0, 1)
        model.NMDAlist.append(NMDA)
        model.ncNMDAlist.append(NC)   
    print 'kinetic NMDA synapses generated'

        
def add_GABAsyns(model, locs=[[0, 0.5]], gmax=0.5, tau1=0.1, tau2=4,
                     rev=-75, NoSynDends=[]):
    model.GABAlist = []
    model.ncGABAlist = []
    gmax = gmax/1000.   # Set in nS and convert to muS
    for loc in locs:
        locInd = int(loc[0])
        if (locInd == -1):
            synloc = model.soma
        else:
            synloc = model.dends[int(loc[0])]
        GABA = h.Exp2Syn(float(loc[1]), sec=synloc) 
        GABA.tau1 = tau1
        GABA.tau2 = tau2
        GABA.e = rev
        if (int(loc[0]) in NoSynDends):
            gg = 0 # same input into a single branch
        else:
            gg = gmax
        NC = h.NetCon(h.nil, GABA, 0, 0, gg)
        model.GABAlist.append(GABA)
        model.ncGABAlist.append(NC)
    print 'inhibitory synapses generated'
        

def addSpines(model):
    for sec in model.dends:
        sec.cm = model.CM*1.5
        sec.g_pas = 1.0/(model.RM/1.5)
    print 'fake spines added'

def hotSpot(model):
    spot = np.ceil(7/2.)
    for section in model.dends:
        s = 0
        for seg in section:
            if s==spot:
                seg.gbar_na = model.gna_dend_hotSpot    
#                print 'hotspot added', s
            else:
                seg.gbar_na = 0
            s+=1


# ----------------------------------------------------------
# SIMULATION RUN
def simulate(model, t_stop=100, NMDA=False, recDend=False, i_recDend=11):
    trec, vrec = h.Vector(), h.Vector()
    gRec, iRec, vDendRec = [], [], []
    gNMDA_rec, iNMDA_rec = [], []
    trec.record(h._ref_t)
    vrec.record(model.soma(0.5)._ref_v)

    # if NMDA:        
    #     for n in np.arange(0, len(model.NMDAlist)):
    #         loc = model.NMDAlist[n].get_loc()
    #         h.pop_section()                        
    #         gNMDA_rec.append(h.Vector())
    #         iNMDA_rec.append(h.Vector())
    #         gNMDA_rec[n].record(model.NMDAlist[n]._ref_g)
    #         iNMDA_rec[n].record(model.NMDAlist[n]._ref_i)
    #     gRec.append(gNMDA_rec)
    #     iRec.append(iNMDA_rec)
    if recDend:
        n = 0
        for i_dend in i_recDend:
            vDendRec.append(h.Vector())
            vDendRec[n].record(model.dends[i_dend](0.5)._ref_v)
            n+=1
 
    h.celsius = model.CELSIUS
    h.finitialize(model.E_PAS)
    neuron.run(t_stop)
    return np.array(trec), np.array(vrec), np.array(vDendRec)
#, np.array(caDendRec), np.array(vSecRec)

