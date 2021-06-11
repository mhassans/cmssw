from array import array

fInPath = '../ntuple_photons_NoPU.root'
fIn = ROOT.TFile.Open(fInPath, "READ")
treeIn = fIn.Get("hgcalTriggerNtuplizer/HGCalTriggerNtuple")

fOut = ROOT.TFile("output.root", "RECREATE")
treeOut = ROOT.TTree("variables", "variables")



ROOT.gStyle.SetOptStat(0)

for entryNum in range(0, treeOut.GetEntries()):
