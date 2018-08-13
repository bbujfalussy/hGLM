import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.cm as cm

from neuron import h

def plot_syns(data, model, plotInh=False):
    #----------------------------------------------------------------------------
    # Some parameters
    nsyn = sum(data.Ensyn)
    if (plotInh): nsyn = int(nsyn + sum(data.Insyn))
    #----------------------------------------------------------------------------
    # Create neuron and generate locations
    cell_secs = []
    for sec in h.allsec(): 
        cell_secs.append(sec)
        sec.nseg = 100
    locs = np.array(data.Elocs)
    if (plotInh): locs = np.vstack((locs, np.array(data.Ilocs)))
    #----------------------------------------------------------------------------
    # Add dummy compartments at synapse locations
    synapses = []
    print ('number of synapses:', nsyn)
    for s in np.arange(0,nsyn):
        synapses.append(h.Section())
        synapses[s].L = 5# 1e-6
        synapses[s].diam = 1#1e-6
        idend = np.int(locs[s,0])
        if (idend > -1):
            synapses[s].connect(model.dends[idend], np.int(locs[s,1]), 0)
        else :
            synapses[s].connect(model.soma, 0.5, 0)
    #----------------------------------------------------------------------------
    # Get numbers for plotting cell morphology
    # cell: [xstart xend ystart yend diamstart diamend]
    # synapse: [x, y]
    h.define_shape()
    cell_coordinates = []
    for sec in cell_secs:
        sec.push()
        for stepCount in np.arange(1, h.n3d()):
            stepCount =  float(stepCount)
            cell_coordinates.append([h.x3d(stepCount-1), h.x3d(stepCount),
                                 h.y3d(stepCount-1), h.y3d(stepCount),
                                     h.diam3d(stepCount-1), h.diam3d(stepCount)])
        h.pop_section()


    set_dends = set(locs[:,0])
    # I don't know how to work with sets, so convert it to an array
    set_dends = np.array(list(set_dends))
    dend_coord = []
    for i_dend in np.arange(0,len(set_dends)):
        syn_coordinates = []
        d_dend = np.int(set_dends[i_dend])
        if (d_dend > -1):
            model.dends[d_dend].push()
        else :
            model.soma.push()
        for stepCount in np.arange(1, h.n3d()):
            stepCount =  float(stepCount)
            syn_coordinates.append([h.x3d(stepCount-1), h.x3d(stepCount),
                                    h.y3d(stepCount-1), h.y3d(stepCount), 
                                    h.diam3d(stepCount-1), h.diam3d(stepCount)])
        h.pop_section()
        dend_coord.append(syn_coordinates)

    #----------------------------------------------------------------------------
    # Make plot
    plt.figure()
    # Cell
    for pt in np.arange(len(cell_coordinates)):
        xstart, xend = cell_coordinates[pt][0], cell_coordinates[pt][1]
        ystart, yend = cell_coordinates[pt][2], cell_coordinates[pt][3]
        lx = xend-xstart
        ly = yend-ystart
        l = np.sqrt(np.dot([lx,ly],[lx,ly]))
        diamstart = cell_coordinates[pt][4] / 2
        diamend = cell_coordinates[pt][5] / 2

        if diamstart > 8: diamstart=diamstart/3.5
        if diamend > 8: diamend=diamend/3.5

        if l>0:
            plt.fill([xstart+2*diamstart*ly/(2*l), xend+2*diamend*ly/(2*l),
                     xend-2*diamend*ly/(2*l), xstart-2*diamstart*ly/(2*l)], 
                     [ystart-2*diamstart*lx/(2*l), yend-2*diamend*lx/(2*l),
                     yend+2*diamend*lx/(2*l), ystart+2*diamstart*lx/(2*l)],
                     'k', lw=0) 

        c1 = plt.Circle((xstart, ystart), diamstart, color='k', lw=0)
        c2 = plt.Circle((xend, yend), diamend, color='k', lw=0) 
        plt.gca().add_patch(c1)
        plt.gca().add_patch(c2)
        
    # Synapses
    nsyn_d = data.Ensyn
    if (plotInh): nsyn_d = nsyn_d + [sum(data.Insyn)] # all inhibitory synapses are plotted with the same color 

    
    idsyn = nsyn_d[:]
    for ii in np.arange(len(idsyn)):
        idsyn[ii] = sum(nsyn_d[0:(ii+1)])

    ecols13 = ['#6d4022', '#a05e32',
               '#512100', '#883600', '#b44800', '#e35b00', '#ff802b',
               '#715b00', '#ac8a00', '#d4aa00', '#f9c800',
               '#6d6649', '#a89f7d']     
    ecols4 = ['#6d4022', '#b44800', '#d4aa00', '#a89f7d']     
    ecols1 = ['#a05e32']     
    if (len(data.Ensyn)==13):
         ecols = ecols13
    elif (len(data.Ensyn)==4):
         ecols = ecols4
    elif (len(data.Ensyn)==1):
         ecols = ecols1
    #else : # a colormap with random colors for each cluster
    cmap = cm.get_cmap(name="jet")
    nden = len(data.Ensyn)
    ecols = [None] * nden
    np.random.seed(seed=16)
    for icol in np.arange(nden):
        icc = np.random.randint(256)
        ecols[icol] = colors.rgb2hex(cmap(icc))
    icols = ['#000000']
    # print (ecols)
    # print (idsyn)
    for syn in np.arange(nsyn):
        # select the branch for the synapse
        idend = np.int(locs[syn,0])
        iidend = np.int(np.argwhere(set_dends == idend))
        syn_coordinates = dend_coord[iidend]
        ptmax=np.shape(syn_coordinates)[0]
        # find the coordinates along the branch
        ptt = locs[syn,1] * ptmax
        pt = np.int(np.floor(ptt))
        ptd = ptt - pt
        xstart, xend = syn_coordinates[pt][0], syn_coordinates[pt][1]
        ystart, yend = syn_coordinates[pt][2], syn_coordinates[pt][3]
        lx = xend-xstart
        ly = yend-ystart
        x = xstart + lx * ptd
        y = ystart + lx * ptd
        nEsyn = sum(data.Ensyn)
        if (syn < nEsyn):
            iden = -1
            for ii in np.arange(len(idsyn)):
                 if syn >= idsyn[ii]: iden = ii
            iden = iden + 1
            # print (syn, iden)
            colsyn = ecols[iden]
            c = plt.Circle((x,y), 3, color=colsyn, lw=0)
        else :
            c = plt.Circle((x,y), 3, color=icols[0], lw=2)
        plt.gca().add_patch(c)   

    plt.axis('equal')
    plt.xlim(-200,200)
    plt.ylim(-160,240)
    plt.show(block=False)
            
