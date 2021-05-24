import os
import sys
import commands
import subprocess

def files_list(dataName, printFileNames = False):
	
	if dataName == "MinBias":
		dataDir = '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/00000/'
	elif dataName == "DYToLL":
		dataDir = "/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/FEVT/PU200_BSzpz35_BSzpz35_111X_mcRun4_realistic_T15_v1-v3/270000/"
	elif dataName == "ZprimeToMuMu":
		dataDir = "/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/FEVT/PU140_111X_mcRun4_realistic_T15_v1-v1/110000/"
	elif dataName == "TT":
		dataDir = "/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v2/280000/"
	elif dataName == "DoubleMuon":
		dataDir = "/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/DoubleMuon_gun_FlatPt-1To100/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/100000/"

	else:
		print "Incorrect dataset"

	rootFileName = dataName + '.root'

	commands = ['eos']
	lsDataDir='/eos/cms' + dataDir
	commands.append('ls ' + lsDataDir)
	dir=subprocess.Popen(commands, stdout=subprocess.PIPE)
	lsOutput=dir.communicate()[0]

	fileList = []
	for f in lsOutput.split():
		fileList.append(dataDir + f)
		if printFileNames:
			print dataDir + f
	return fileList, rootFileName

