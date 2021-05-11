import FWCore.ParameterSet.Config as cms
import os
import sys
import commands

process = cms.Process("muonAnalyzer")

# MessageLogger & co.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False), SkipEvent = cms.untracked.vstring('ProductNotFound'))

# input files (up to 255 files accepted)

#fileName = "root://eoscms.cern.ch//eos/cms/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/00000/0760B9BE-600E-9C4D-9878-93D13947D728.root"
#filName = '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/00000/0760B9BE-600E-9C4D-9878-93D13947D728.root'

fileName = "/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/DYToLL_M-10To50_TuneCP5_14TeV-pythia8/FEVT/PU200_BSzpz35_BSzpz35_111X_mcRun4_realistic_T15_v1_ext1-v3/50000/0126A814-62FF-5D4D-872C-73999BCE3934.root"

process.source = cms.Source('PoolSource',
fileNames = cms.untracked.vstring( 
    fileName 
                                  ),
skipEvents =  cms.untracked.uint32(0),
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))
 
# PostLS1 geometry used
process.load('Configuration.Geometry.GeometryExtended2016Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

####Event Setup Producer, OMTF
process.load('L1Trigger.L1TMuonOverlap.fakeOmtfParams_cff')
process.load('L1Trigger.L1TMuonOverlap.simOmtfDigis_cfi')
process.simOmtfDigis.dumpResultToXML=cms.bool(True)

####Tracking partiles configuration https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideTrackingTruth
process.trackingParticle = cms.EDAnalyzer(
	"MuonAnalyzer",
	trackingTruth = cms.untracked.InputTag('mix', 'MergedTrackTruth'),
	dumpVertexes = cms.untracked.bool(True),
	dumpOnlyBremsstrahlung = cms.untracked.bool(True)
)


#process.MuonAnalizer = cms.EDAnalyzer("MuonAnalizer1")
process.MyPath = cms.Path(process.simOmtfDigis + process.trackingParticle)
#process.MyPath = cms.Path(process.simOmtfDigis + process.trackingParticle + process.MuonAnalizer)
#process.MyPath = cms.Path(process.simOmtfDigis  + process.MuonAnalizer)
process.schedule = cms.Schedule(process.MyPath)
