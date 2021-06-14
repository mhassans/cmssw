from array import array
import ROOT
import math
import numpy as np
from funcs import getEtaPhiBins, sumTowers

fInPath = '../ntuple_photons_NoPU.root'
fIn = ROOT.TFile.Open(fInPath, "READ")
treeIn = fIn.Get("hgcalTriggerNtuplizer/HGCalTriggerNtuple")

fOut = ROOT.TFile("output.root", "RECREATE")
treeOut = ROOT.TTree("variables", "variables")

floatVars = ['GenEnergy', 'GenEta', 'GenPhi', 'IntegralEM', 'IntegralHad', 'EM1x1', 'Had1x1', 'EM3x3', 'Had3x3', 'EM5x5', 'Had5x5']

for etaSide in ['zPlus_', 'zMinus_']:
    for var in floatVars:
        exec(etaSide + var + " = array('f', [-100.])")
        exec("treeOut.Branch('{0}', {0}, '{0}/F')".format(etaSide + var))

ROOT.gStyle.SetOptStat(0)

etaBinStep = 0.0870
minBinEta = -35
maxBinEta = 35
minEta = minBinEta * etaBinStep
maxEta = maxBinEta * etaBinStep
nBinsEta = maxBinEta - minBinEta

phiBinStep = 2*math.pi/72
minBinPhi = -36
maxBinPhi = 36
minPhi = minBinPhi * phiBinStep
maxPhi = maxBinPhi * phiBinStep
nBinsPhi = maxBinPhi - minBinPhi

histEM = ROOT.TH2D("histEM","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
histHad = ROOT.TH2D("histHad","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)


treeInVars = ['tower_iPhi', 'tower_iEta', 'tower_etEm', 'tower_etHad', 'tower_n', 'tower_eta', 'gen_eta', 'gen_phi', 'gen_energy']
#for entryNum in range(0, treeIn.GetEntries()):
for entryNum in range(0, 100):
    treeIn.GetEntry(entryNum)
    for var in treeInVars:
        exec("{0} = getattr(treeIn, '{0}')".format(var))
    histEM.Reset()
    histHad.Reset()

    for towerID in range(tower_n):
        etaPhiBins = getEtaPhiBins(tower_eta[towerID], tower_iEta[towerID], tower_iPhi[towerID])
        histEM.SetBinContent(etaPhiBins[0], etaPhiBins[1], tower_etEm[towerID])
        histHad.SetBinContent(etaPhiBins[0], etaPhiBins[1], tower_etHad[towerID])

    ZplusIndex = 0 if (int(np.sign(gen_eta[0]))==1) else 1

    zPlus_GenEnergy[0] = gen_energy[ZplusIndex]/np.cosh(gen_eta[ZplusIndex])
    zPlus_GenEta[0] = gen_eta[ZplusIndex]
    zPlus_GenPhi[0] = gen_phi[ZplusIndex]
    zMinus_GenEnergy[0] = gen_energy[1-ZplusIndex]/np.cosh(gen_eta[1-ZplusIndex])
    zMinus_GenEta[0] = gen_eta[1-ZplusIndex]
    zMinus_GenPhi[0] = gen_phi[1-ZplusIndex]

    zPlus_IntegralEM[0] = histEM.Integral(int(nBinsEta/2), nBinsEta, 1, nBinsPhi)
    zPlus_IntegralHad[0] = histHad.Integral(int(nBinsEta/2), nBinsEta, 1, nBinsPhi)
    zMinus_IntegralEM[0] = histEM.Integral(1, int(nBinsEta/2), 1, nBinsPhi)
    zMinus_IntegralHad[0] = histHad.Integral(1, int(nBinsEta/2), 1, nBinsPhi)

    zPlus_EM1x1[0] = sumTowers(histEM, gen_eta[ZplusIndex], gen_phi[ZplusIndex], numNeighbors=0)
    zPlus_EM3x3[0] = sumTowers(histEM, gen_eta[ZplusIndex], gen_phi[ZplusIndex], numNeighbors=1)
    zPlus_EM5x5[0] = sumTowers(histEM, gen_eta[ZplusIndex], gen_phi[ZplusIndex], numNeighbors=2)
    zPlus_Had1x1[0] = sumTowers(histHad, gen_eta[ZplusIndex], gen_phi[ZplusIndex], numNeighbors=0)
    zPlus_Had3x3[0] = sumTowers(histHad, gen_eta[ZplusIndex], gen_phi[ZplusIndex], numNeighbors=1)
    zPlus_Had5x5[0] = sumTowers(histHad, gen_eta[ZplusIndex], gen_phi[ZplusIndex], numNeighbors=2)
    zMinus_EM1x1[0] = sumTowers(histEM, gen_eta[1-ZplusIndex], gen_phi[1-ZplusIndex], numNeighbors=0)
    zMinus_EM3x3[0] = sumTowers(histEM, gen_eta[1-ZplusIndex], gen_phi[1-ZplusIndex], numNeighbors=1)
    zMinus_EM5x5[0] = sumTowers(histEM, gen_eta[1-ZplusIndex], gen_phi[1-ZplusIndex], numNeighbors=2)
    zMinus_Had1x1[0] = sumTowers(histHad, gen_eta[1-ZplusIndex], gen_phi[1-ZplusIndex], numNeighbors=0)
    zMinus_Had3x3[0] = sumTowers(histHad, gen_eta[1-ZplusIndex], gen_phi[1-ZplusIndex], numNeighbors=1)
    zMinus_Had5x5[0] = sumTowers(histHad, gen_eta[1-ZplusIndex], gen_phi[1-ZplusIndex], numNeighbors=2)

    treeOut.Fill()

fOut.Write()
fOut.Close()
