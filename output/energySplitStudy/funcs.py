import sys

ROOT.gStyle.SetOptStat(0)

def getEtaPhiBins(tower_eta, tower_iEta, tower_iPhi):
    etaBin = -100
    phiBin = -100
    if (int(np.sign(tower_eta))==1):
        etaBin = 51 + tower_iEta
        phiBin = 72 - tower_iPhi
    elif (int(np.sign(tower_eta))==-1):
        etaBin = 20 - tower_iEta
        phiBin = 1 + (tower_iPhi+36)%72
    else:
        print('ERROR: found eta=0?')
        sys.exit(1)
    return (etaBin, phiBin)

def sumTowers(hist, gen_eta, gen_phi, index, numNeighbors):
    """
    numNeighbors = 1 means 9 (=3x3) towers with 8 around the center + 1 center
    numNeighbors = 2 means 25 (=5x5) towers with 24 around the center tower + 1 center
    """
    if int(numNeighbors)!=numNeighbors:
        print("ERROR: number of towers must be an integer!")
        sys.exit(1)
    
    genEtaBin = hist.GetXaxis().FindBin(abs(gen_eta[index]))
    genPhiBin = hist.GetYaxis().FindBin(gen_phi[index])
    
    sums = 0
    for ix in range(-numNeighbors, 1+numNeighbors):
        for iy in range(-numNeighbors, 1+numNeighbors):
            tower_eta = ix+genEtaBin
            tower_phi = iy+genPhiBin
            if tower_phi > 72:
                tower_phi = tower_phi - 72
            elif tower_phi < 1:
                tower_phi = tower_phi + 72
            sums += hist.GetBinContent(tower_eta, tower_phi)
            ##########test###########
            #x = hist.GetBinContent(ix+genEtaBin, iy+genPhiBin)
            #if x>5:
            #    print("sum+=", x, ix+genEtaBin, iy+genPhiBin)           
            #########test finish#########

    return sums

