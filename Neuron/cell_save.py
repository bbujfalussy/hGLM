import pickle
import saveClass as sc
import libcell as lb
import numpy as np
import struct
import os

def save_syn(data, outdir='data', clust='0'):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    filename='clust'+str(clust)
    mat = np.asarray(data.Ensyn)
    mat = np.expand_dims(mat, axis=1)


    bfname = './'+outdir+'/Syn_'+filename+'.bin'
    print bfname
    # create a binary file
    binfile = file(bfname, 'wb')
    # and write out two integers with the row and column dimension
    header = struct.pack('2I', mat.shape[0], mat.shape[1])
    binfile.write(header)
    # then loop over columns and write each
    for i in range(mat.shape[1]):
        ddata = struct.pack('%id' % mat.shape[0], *mat[:,i])
        binfile.write(ddata)
    binfile.close()

def save_sim(data, out_binary=False, out_vdend=False, out_pickle=False, outdir='data', rate='0', omega='0', iAMP='0'):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    modelData = sc.emptyObject()
    lb.props(modelData)
    description = 'Simulation type:'+data.simType+'. Fixed patterns,  distributed.'

    if data.NMDA:
        data.Ntype='NMDA'
    else:
        data.Ntype='AMPA'
    if data.NMDAkinetic:
        data.Ntype='NMDAk'
        
    if data.ACTIVE:
        if data.ACTIVEdend:
            data.act='Adend'
        else :
            data.act='Asoma'
    else:
        data.act='passive'


    if (data.stimType == 'allDend'): 
        filename=data.Ntype+'_ApN'+str(data.ApN)+'_'+data.locBias+'_'+data.act+'_r'+str(rate)+'_o'+str(omega)+'_g'
    elif (data.stimType == 'mixedori'): 
        filename=data.Ntype+'_ApN'+str(data.ApN)+'_'+data.locBias+'_'+data.act+'_r'+str(rate)+'_o'+str(omega)+'_i'+str(iAMP)+'_g'
        filename = filename + '_b' + str(data.inburstrate)
    elif (data.stimType == 'orientations'): 
        filename=data.Ntype+'_ApN'+str(data.ApN)+'_'+data.locBias+'_'+data.act+'_r'+str(rate)+'_o'+str(omega)+'_i'+str(iAMP)+'_g'
    elif (data.stimType == 'minis'): 
        filename=data.Ntype+'_ApN'+str(data.ApN)+'_'+data.locBias+'_'+data.act+'_r'+str(rate)+'_o'+str(omega)+'_i'+str(iAMP)+'_g'
    else :
        filename=data.Ntype+'_NAR'+str(data.NAR)+'_'+data.locBias+'_'+data.act+'_r'+str(rate)+'_o'+str(omega)+'_i'+str(iAMP)+'_g'
    #filename='temp'

    if out_pickle:
        dataList = [data, modelData]
        fname = './'+outdir+'/'+filename+'.pkl'
        f = open(fname, 'wb')
        pickle.dump(dataList, f)
        f.close()


    if out_binary:
        #---------------------------------------------
        # WRITE the response in a binary file to read it with R
        mat = np.array(data.vdata)
        bfname = './'+outdir+'/vdata_'+filename+'.bin'
        print bfname
        # create a binary file
        binfile = file(bfname, 'wb')
        # and write out two integers with the row and column dimension
        header = struct.pack('2I', mat.shape[0], mat.shape[1])
        binfile.write(header)
        # then loop over columns and write each
        for i in range(mat.shape[1]):
            ddata = struct.pack('%id' % mat.shape[0], *mat[:,i])
            binfile.write(ddata)
        binfile.close()

        if out_vdend:
            # WRITE the dendritic response
            nRep = len(data.vDdata)
            mat = np.array(data.vDdata[0])
            for i in range(1, nRep):
                mat = np.hstack((mat, data.vDdata[i]))
            
            bfname = './'+outdir+'/vDdata_'+filename+'.bin'
            # create a binary file
            binfile = file(bfname, 'wb')
            # and write out two integers with the row and column dimension
            header = struct.pack('2I', mat.shape[0], mat.shape[1])
            binfile.write(header)
            # then loop over columns and write each
            for i in range(mat.shape[1]):
                ddata = struct.pack('%id' % mat.shape[0], *mat[:,i])
                binfile.write(ddata)
            binfile.close()
        

        #---------------------------------------------
        # WRITE the location of the synapses        
        if (data.SYN):
            if (data.GABA) :
                Ilocs = np.array(data.Ilocs) 
                Ilocs[:,1] = 1 + Ilocs[:,1] # code that these are inhibitory synapses
                Elocs = np.array(data.Elocs)
                Locs = np.row_stack((Elocs, Ilocs))
            else :
                Locs = np.array(data.Elocs)

            bfname = './'+outdir+'/synlocs_'+filename+'.bin'
            # create a binary file
            binfile = file(bfname, 'wb')
            # and write out two integers with the row and column dimension
            header = struct.pack('2I', Locs.shape[0], Locs.shape[1])
            binfile.write(header)
            # then loop over columns and write each
            for i in range(Locs.shape[1]):
                ddata = struct.pack('%id' % Locs.shape[0], *Locs[:,i])
                binfile.write(ddata)
            binfile.close()

        #---------------------------------------------
        # Write the input spike train
        if (len(data.stim)>0):
            stim = data.stim
            bfname = './'+outdir+'/stim_'+filename+'.bin'
            # create a binary file
            binfile = file(bfname, 'wb')
            # and write out two integers with the row and column dimension
            header = struct.pack('2I', stim.shape[0], stim.shape[1])
            binfile.write(header)
            # then loop over columns and write each
            for i in range(stim.shape[1]):
                ddata = struct.pack('%id' % stim.shape[0], *stim[:,i])
                binfile.write(ddata)
            binfile.close()
