import FWCore.ParameterSet.Config as cms
import os
import sys
import commands
import subprocess
import files_list.py

process = cms.Process("muonAnalyzer")

# MessageLogger & co.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False), SkipEvent = cms.untracked.vstring('ProductNotFound'))

# input files (up to 255 files accepted)

#fileName = "root://eoscms.cern.ch//eos/cms/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/00000/0760B9BE-600E-9C4D-9878-93D13947D728.root"
fileName1 = '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/00000/0760B9BE-600E-9C4D-9878-93D13947D728.root'

fileName2 = "/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/DYToLL_M-10To50_TuneCP5_14TeV-pythia8/FEVT/PU200_BSzpz35_BSzpz35_111X_mcRun4_realistic_T15_v1_ext1-v3/50000/0126A814-62FF-5D4D-872C-73999BCE3934.root"

fileName3 = "/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/FEVT/PU200_BSzpz35_BSzpz35_111X_mcRun4_realistic_T15_v1-v3/270000/00611D76-B4B6-A046-B6B5-2DF8CD92FE5D.root"

fileName4 = "/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/FEVT/PU140_111X_mcRun4_realistic_T15_v1-v1/110000/A0BB7A6E-7CA3-854A-8F67-F0C185CC6424.root"

fileName5 = "/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v2/280000/003ACFBC-23B2-EA45-9A12-BECFF07760FC.root"


fileName6 = "/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/DoubleMuon_gun_FlatPt-1To100/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/100000/00554A1E-61A5-DE49-A5F3-B47F2F4039EF.root"

outRootFileName = "test1.root"

process.source = cms.Source('PoolSource',
fileNames = cms.untracked.vstring( 
    fileName4
                                  ),
skipEvents =  cms.untracked.uint32(0),
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1))
 
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
	rootFileName = cms.untracked.string(outRootFileName),
	dumpVertexes = cms.untracked.bool(True),
	dumpOnlyBremsstrahlung = cms.untracked.bool(True)
)


#process.MuonAnalizer = cms.EDAnalyzer("MuonAnalizer1")
process.MyPath = cms.Path(process.simOmtfDigis + process.trackingParticle)
#process.MyPath = cms.Path(process.simOmtfDigis + process.trackingParticle + process.MuonAnalizer)
#process.MyPath = cms.Path(process.simOmtfDigis  + process.MuonAnalizer)
process.schedule = cms.Schedule(process.MyPath)
