import FWCore.ParameterSet.Config as cms 

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process('DIGI',Phase2C9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2000)
)

# Input source
filePfx = '/store/mc/Phase2HLTTDRWinter20DIGI/SingleElectron_PT2to200/GEN-SIM-DIGI-RAW/PU200_110X_mcRun4_realistic_v3_ext2-v2/40000/'
process.source = cms.Source("PoolSource",
       fileNames = cms.untracked.vstring(filePfx + '/00582F93-5A2A-5847-8162-D81EE503500F.root',
                                         filePfx + '/01C92E48-6B8E-6D4B-8B57-57540AFBC9D5.root',
                                         filePfx + '/01F69B83-AC48-1242-9C6D-339335FABC1E.root', 
                                         filePfx + '/022471A5-F18F-3D45-AAF0-88BAAA4B159C.root', 
                                         filePfx + '/025E005F-AF77-9F45-AC99-4A78ADBAFF40.root', 
                                         filePfx + '/0287BD18-6928-8844-9919-E2183E3C5BF4.root', 
                                         filePfx + '/02C1307B-78C5-BC4F-B854-6988CD65D4C8.root', 
                                         filePfx + '/03B69335-90DB-EA48-9D5C-F88786E854F0.root',
                                         filePfx + '/0457C237-9151-ED46-8FD5-F1B152AA824F.root', 
                                         filePfx + '/0457CFA7-3B54-324F-B61F-CA1F9384266B.root'),
       inputCommands=cms.untracked.vstring(
           'keep *',
           'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
           'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
           'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
           'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
           'drop l1tEMTFTrack2016s_simEmtfDigis__HLT',
           'drop FTLClusteredmNewDetSetVector_mtdClusters_FTLBarrel_RECO',
           'drop FTLClusteredmNewDetSetVector_mtdClusters_FTLEndcap_RECO',
           'drop MTDTrackingRecHitedmNewDetSetVector_mtdTrackingRecHits__RECO',
           'drop BTLDetIdBTLSampleFTLDataFrameTsSorted_mix_FTLBarrel_HLT',
           'drop ETLDetIdETLSampleFTLDataFrameTsSorted_mix_FTLEndcap_HLT',
           )
       )

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('SingleElectronPt10_cfi nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("ntuple.root")
    )

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

# load HGCAL TPG simulation
process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')

from L1Trigger.L1THGCal.customTriggerSums import custom_full_trigger_sums
process = custom_full_trigger_sums(process)

from L1Trigger.L1THGCal.customTriggerGeometry import custom_geometry_decentralized_V11
process = custom_geometry_decentralized_V11(process, implementation=2)

from L1Trigger.L1THGCal.customTowers import custom_towers_energySplit
process = custom_towers_energySplit(process)

process.hgcl1tpg_step = cms.Path(process.hgcalTriggerPrimitives)


# load ntuplizer
process.load('L1Trigger.L1THGCalUtilities.hgcalTriggerNtuples_cff')
process.ntuple_step = cms.Path(process.hgcalTriggerNtuples)

# Schedule definition
process.schedule = cms.Schedule(process.hgcl1tpg_step, process.ntuple_step)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

