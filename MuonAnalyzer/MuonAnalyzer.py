import FWCore.ParameterSet.Config as cms
import os
import sys
import commands
import subprocess
from files_list import files_list

process = cms.Process("muonAnalyzer")

# MessageLogger & co.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False), SkipEvent = cms.untracked.vstring('ProductNotFound'))

# input files (up to 255 files accepted)

if len(sys.argv) == 3:
	fileNamesList, outRootFileName = files_list(sys.argv[2])
	n = 1000
elif len(sys.argv) == 4:
	fileNamesList, outRootFileName = files_list(sys.argv[2], False)
	n = int(sys.argv[3])
else:
	fileNamesList, outRootFileName = files_list("MinBias")
	n = 1000;

process.source = cms.Source('PoolSource',
fileNames = cms.untracked.vstring( 
    fileNamesList
                                  ),
skipEvents =  cms.untracked.uint32(0),
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(n))
 
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
