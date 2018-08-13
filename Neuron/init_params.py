
if (data.simType == 'nsyn'):
    data.st_onset = 0.1 # in seconds
    data.st_duration = 0.4
    data.TSTOP = 0.35

    data.GABA = False

    data.Ensyn = [40] # 160
    data.Insyn = [40] # 160
    data.locBias = 'none'
    data.locDend = [95] # 1dend: 11; 4dend: 11, 24, 70, 95
    data.tInterval = 1

    # # for minis
    # data.Ensyn = [49, 26, 45, 40] # 160
    # data.Insyn = [16, 5, 2, 5, 4] # 32 - 20 %
    # data.locDend = [11, 24, 70, 95] # 1dend: 11; 4dend: 11, 24, 70, 95
    # data.GABA = True

if (data.simType == 'nsynPoiss'):
    data.GABA = False

    data.Ensyn = [40] # 160
    data.Erate = 1
    burstRate = 3.2
    maxN = 20

    data.locBias = 'midle'
    data.locDend = [95] # 1dend: 11; 4dend: 11, 24, 70, 95
    data.locDendRec = [38, 47, 69, 70] 

    data.iAmp = maxN
    data.ind_Rates = burstRate
    iomega = data.Erate

if (data.simType == 'iSteps'):
    data.st_onset = .25 # in seconds
    data.st_duration = .5
    data.TSTOP = 1

    data.ICLAMP = True
#    data.iclampLoc = ['dend', 0.5, 28]
    data.iclampLoc = ['soma', 0.5]
    data.iclampOnset = data.st_onset * 1000
    data.iclampDur = data.st_duration * 1000
    data.iclampAmp = 0.1
    data.iRange = np.arange(0.24,0.48,0.04)

    data.EbGroundRate = 4
    data.IbGroundRate = 10

else:
    data.ICLAMP = False

    

if (data.simType == 'rates'):
    # to reproduce Spencer's data I used this:
    data.st_onset = .25 # in seconds
    data.st_duration = .5
    data.TSTOP = 1

    data.EbGroundRate = 2
    data.IbGroundRate = 2
    data.rateE = np.array([2,4,6,8,10]) # excess rate. true rate: rate + bground! 
    data.rateI = np.array([2,4,6,8,10]) # 

    data.ICLAMP = True
    data.iclampLoc = ['soma', 0.5]
    data.iclampOnset = data.st_onset * 1000
    data.iclampDur = data.st_duration * 1000
    data.iclampAmp = 0.2
    data.locDendRec = [38, 47, 69, 70] 


if (data.simType == 'random'):
    data.st_onset = 0 # in seconds
    data.st_duration = data.TSTOP

    data.EbGroundRate = 0
    data.IbGroundRate = 0
    data.rateE = np.array([2,4,6,8,10]) # excess rate. true rate: rate + bground! 
    data.rateI = np.array([2,4,6,8,10]) # 


if (data.simType == 'allDend'):
    # not excess rates but true rates!
    #                         fgErate,        fgIrate,              bgErate,         bgIrate
    ### balanceRange = np.array([[20, 26, 30, 36],[30, 80, 120, 200], [5, 6, 7, 9], [20, 70, 80, 140]]) # this was the default!
    ## balanceRange = np.array([[16, 24, 32, 50],[60, 100, 200, 400], [4, 4, 4, 4], [15, 15, 15, 15]])
    balanceRange = np.array([[5, 10, 20, 26, 40, 55],[7, 15, 30, 80, 100, 300], [1, 3, 5, 6, 8, 10], [5, 10, 20, 55, 80, 150]])
    if (data.stimType == 'orientations'):
		balanceRange = np.array([[20, 26, 30, 36],[30, 80, 120, 140], [5, 6, 7, 9], [20, 70, 80, 100]])
    if (data.stimType == 'mixedori'):
        balanceRange = np.array([[20, 26, 30, 36],[30, 80, 120, 140], [5, 6, 7, 9], [20, 70, 80, 100]])
    if (data.actType == 'Poisson'):
        if (data.ACTIVEdend == True):
            data.EIrates = balanceRange[:,2]
        else:
            data.EIrates = balanceRange[:,1]            
    else :
        data.EIrates = balanceRange[:,data.ind_Rates]

    data.omega=[data.omegai, 10] # omegaUp, omegaDown (Hz)
    
    if (data.stimType == 'orientations'):
        data.ICLAMP = True
        if (data.ind_Rates < 2):
            clamp_dend = 24
        else :
            clamp_dend = 70
        # data.iclampLoc = ['dend', 0.5, clamp_dend]
        data.iclampLoc = ['soma', 0.5]
        data.iclampOnset = 0
        data.iclampDur = data.TSTOP * 1000
        data.iclampAmp = 0.01 * data.iAmp
    else :
        data.ICLAMP = False

    data.EbGroundRate = 0
    data.IbGroundRate = 0

    data.rateE = np.array([3, 5, 7, 9, 11, 13]) # excess rate. true rate: rate + bground! 
    data.rateI = np.array([3, 14, 24, 34, 44, 54]) #

    ##
# worked fine for 10 um E and 100 um I synapses + 10 I synapse on the soma
#    data.rateE = np.array([2,3,4,5,6]) # excess rate. true rate: rate + bground! 
#    data.rateI = np.array([0.1,5,16,26,38]) # 


if (data.stimType == 'balanceIteration'):
    # Specific parameters
    balanceRange = np.array([[12,18,24,36],[8, 30, 60, 100], [4,4,4,4], [8, 12, 14, 16]])
#    data.balanceRange = np.array([[11,13,15,17],[10, 30, 50, 70], [4,4,4,4], [ibg, ibg, ibg, ibg]])
    #                              fgErate,        fgIrate,      bgErate,   bgIrate
    # not excess rates but true rates!
    data.EIrates = balanceRange[:,data.ind_Rates]
    data.omega=[data.omegai, 10] # omegaUp, omegaDown (Hz)
    
    data.ICLAMP = True
#    data.iclampLoc = ['dend', 0.5, 28]
    data.iclampLoc = ['soma', 0.5]
    data.iclampOnset = 0
    data.iclampDur = data.TSTOP * 1000
    data.iclampAmp = 0.1

    if (data.simType == 'clustered'):        
        if (data.omegai < 0.1):
            data.iclampAmp = 0.175
            if (data.ind_Rates == 1): data.iclampAmp = 0.185
            if (data.ind_Rates == 2): data.iclampAmp = 0.2
            if (data.ind_Rates == 3): data.iclampAmp = 0.2

        if (data.omegai == 0.5): data.iclampAmp = 0.15

        if (data.omegai == 1):
            data.iclampAmp = 0.12 #* 0
            # if (data.ind_Rates == 1): data.iclampAmp = 0.07
            # if (data.ind_Rates == 2): data.iclampAmp = 0.07
            # if (data.ind_Rates == 3): data.iclampAmp = 0.07

        if (data.omegai == 2):
            data.iclampAmp = 0.07
            if (data.ind_Rates == 0): data.iclampAmp = 0.09

        if (data.omegai == 4):
            data.iclampAmp = 0.05
            if (data.ind_Rates == 1): data.iclampAmp = 0.03
            if (data.ind_Rates == 2): data.iclampAmp = 0.01
            if (data.ind_Rates == 3): data.iclampAmp = 0.01

        if (data.omegai == 10):
            data.iclampAmp = 0
            if (data.ind_Rates == 1): data.iclampAmp = -0.04
            if (data.ind_Rates == 2): data.iclampAmp = -0.05
            if (data.ind_Rates == 3): data.iclampAmp = -0.08

        # if (data.ApN < 0.8): data.iclampAmp = 0.01
        # if (data.ApN > 0.8): data.iclampAmp = -0.07

    
    if (data.simType == 'local'):
        data.locBias = 'midle'
        data.iclampAmp = 0.14
        # if (data.ind_Rates == 1): data.iclampAmp = 0.07
        # if (data.ind_Rates == 2): data.iclampAmp = 0.07
        if (data.ind_Rates == 3): data.iclampAmp = 0.17
        # if (data.ApN < 0.8): data.iclampAmp = 0.05
        # if (data.ApN > 0.8): data.iclampAmp = 0.01

    if (data.simType == 'distributed'):
        data.iclampAmp = 0.13
        if (data.ind_Rates == 1): data.iclampAmp = 0.09
        if (data.ind_Rates == 2): data.iclampAmp = 0.08
        if (data.ind_Rates == 3): data.iclampAmp = 0.06
        # if (data.ApN < 0.8): data.iclampAmp = 0.01
        # if (data.ApN > 0.8): data.iclampAmp = -0.07

    if (data.simType == 'somatic'):
        data.iclampAmp = 0.13 # 
        if (data.ind_Rates == 1): data.iclampAmp = 0.1
        if (data.ind_Rates == 2): data.iclampAmp = 0.09
        if (data.ind_Rates == 3): data.iclampAmp = 0.02
        # if (data.ApN < 0.8): data.iclampAmp = 0.03
        # if (data.ApN > 0.8): data.iclampAmp = -0.4

    if (data.simType == 'Poisson'):
        data.iclampAmp = 0.175
        if (data.ind_Rates == 1): data.iclampAmp = 0.185
        if (data.ind_Rates == 2): data.iclampAmp = 0.2
        if (data.ind_Rates == 3): data.iclampAmp = 0.2
        # if (data.ApN < 0.8): data.iclampAmp = 0.02
        # if (data.ApN > 0.8): data.iclampAmp = -0.05
    



if data.locBias == 'proximal': 
    data.locSeg = [0.0,0.1]
#    data.locBorders = [30,100]
    data.dendSpread = data.locSeg
if data.locBias == 'distal':
    data.locSeg = [0.7,0.9]
#    data.locBorders = [100,250]
    data.dendSpread = data.locSeg
if data.locBias == 'midle':
    data.locSeg = [0.4,0.6]
#    data.locBorders = [100,250]
    data.dendSpread = data.locSeg
if data.locBias == 'none':
    data.locSeg = [0.0,1.0]
 #   data.locBorders = [10,350]
    data.dendSpread = data.locSeg


if (data.simType == 'nsyn') : nseg_den = 50
else: nseg_den = 10
