import FWCore.ParameterSet.Config as cms
#Tracks without extra and hits

#AOD content
RecoTrackerAOD = cms.PSet(
    outputCommands = cms.untracked.vstring(
	'keep recoTracks_ctfWithMaterialTracksP5_*_*',
        'keep recoTracks_ctfWithMaterialTracksP5LHCNavigation_*_*',
        'keep recoTracks_rsWithMaterialTracksP5_*_*',
        'keep recoTracks_cosmictrackfinderP5_*_*',
        'keep recoTracks_beamhaloTracks_*_*',
        'keep recoTracks_splittedTracksP5_*_*',
        'keep recoTracks_ctfWithMaterialTracksP5Top_*_*',
        'keep recoTracks_rsWithMaterialTracksP5Top_*_*',
        'keep recoTracks_cosmictrackfinderP5Top_*_*',
        'keep recoTracks_ctfWithMaterialTracksP5Bottom_*_*',
        'keep recoTracks_rsWithMaterialTracksP5Bottom_*_*',
        'keep recoTracks_cosmictrackfinderP5Bottom_*_*',
        'keep recoTracks_regionalCosmicTracks_*_*',
        'keep *_dedxHitInfo_*_*',
        'keep *_dedxHarmonic2_*_*',
        'keep *_dedxHitInfoCTF_*_*',
        'keep *_dedxHarmonic2CTF_*_*',
        'keep *_dedxHitInfoCosmicTF_*_*',
        'keep *_dedxHarmonic2CosmicTF_*_*')
)

#RECO content
RecoTrackerRECO = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep recoTrackExtras_ctfWithMaterialTracksP5_*_*',
        'keep TrackingRecHitsOwned_ctfWithMaterialTracksP5_*_*',
        'keep recoTrackExtras_ctfWithMaterialTracksP5LHCNavigation_*_*',
        'keep TrackingRecHitsOwned_ctfWithMaterialTracksP5LHCNavigation_*_*',
        'keep recoTrackExtras_rsWithMaterialTracksP5_*_*',
        'keep TrackingRecHitsOwned_rsWithMaterialTracksP5_*_*',
        'keep recoTrackExtras_cosmictrackfinderP5_*_*',
        'keep TrackingRecHitsOwned_cosmictrackfinderP5_*_*',
        'keep recoTrackExtras_beamhaloTracks_*_*',
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*',
        'keep recoTrackExtras_splittedTracksP5_*_*',
        'keep TrackingRecHitsOwned_splittedTracksP5_*_*',
        'keep recoTrackExtras_ctfWithMaterialTracksP5Top_*_*',
        'keep TrackingRecHitsOwned_ctfWithMaterialTracksP5Top_*_*',
        'keep recoTrackExtras_rsWithMaterialTracksP5Top_*_*',
        'keep TrackingRecHitsOwned_rsWithMaterialTracksP5Top_*_*',
        'keep recoTrackExtras_cosmictrackfinderP5Top_*_*',
        'keep TrackingRecHitsOwned_cosmictrackfinderP5Top_*_*',
        'keep recoTrackExtras_ctfWithMaterialTracksP5Bottom_*_*',
        'keep TrackingRecHitsOwned_ctfWithMaterialTracksP5Bottom_*_*',
        'keep recoTrackExtras_rsWithMaterialTracksP5Bottom_*_*',
        'keep TrackingRecHitsOwned_rsWithMaterialTracksP5Bottom_*_*',
        'keep recoTrackExtras_cosmictrackfinderP5Bottom_*_*',
        'keep TrackingRecHitsOwned_cosmictrackfinderP5Bottom_*_*',
        'keep recoTrackExtras_regionalCosmicTracks_*_*',
        'keep TrackingRecHitsOwned_regionalCosmicTracks_*_*',
        'keep *_dedxTruncated40_*_*',
        'keep *_dedxTruncated40CTF_*_*',
        'keep *_dedxTruncated40CosmicTF_*_*',
        'keep recoTracks_cosmicDCTracks_*_*',
        'keep recoTrackExtras_cosmicDCTracks_*_*',
        'keep TrackingRecHitsOwned_cosmicDCTracks_*_*')
)
RecoTrackerRECO.outputCommands.extend(RecoTrackerAOD.outputCommands)

#Full Event content 
RecoTrackerFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring()
)
RecoTrackerFEVT.outputCommands.extend(RecoTrackerRECO.outputCommands)
