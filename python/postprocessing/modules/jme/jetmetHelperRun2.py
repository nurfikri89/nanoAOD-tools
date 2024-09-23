#!/usr/bin/env python
import os
import sys
import subprocess

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.fatJetUncertainties import *

# JEC dict
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC#Recommended_for_MC
jecVersionsMC = {
    'UL2016_preVFP': 'Summer19UL16APV_V7_MC',#106X_mcRun2_asymptotic_preVFP_v11
    'UL2016': 'Summer19UL16_V7_MC', #106X_mcRun2_asymptotic_v17
    'UL2017': 'Summer19UL17_V5_MC', #106X_mc2017_realistic_v10
    'UL2018': 'Summer19UL18_V5_MC', #106X_upgrade2018_realistic_v15_L1v1
    '2022'     :'Summer22_22Sep2023_V2_MC',
    '2022_EE'  :'Summer22EE_22Sep2023_V2_MC',
    '2023'     :'Summer23Prompt23_V1_MC',
    '2023_BPix':'Summer23BPixPrompt23_V1_MC',
    '2024':     'Winter24Prompt24_V5M_MC',
}

jecVersionsFastSim = {
    '2016': 'Spring16_25nsFastSimV1_MC',
    '2017': 'Fall17_FastSimV1_MC',
    '2018': 'Autumn18_FastSimV1_MC',
}

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC#Recommended_for_Data
# NOTE: Always check that for a particular year, we must have one entry
# in "archiveTagsDATA".
archiveTagsDATA = {
    'UL2016_preVFP': 'Summer19UL16APV_V7_DATA',
    'UL2016':  'Summer19UL16_V7_DATA',
    'UL2017':  'Summer19UL17_V5_DATA',
    'UL2018':  'Summer19UL18_V5_DATA',
    '2022':    'Summer22_22Sep2023_V2_DATA',
    '2022_EE': 'Summer22EE_22Sep2023_V2_DATA',
    '2023Cv123': 'Summer23Prompt23_RunCv123_V1_DATA',
    '2023Cv4':   'Summer23Prompt23_RunCv4_V1_DATA',
    '2023_BPix': 'Summer23BPixPrompt23_RunD_V1_DATA',
    '2024':      'Winter24Prompt24_RunBCDEF_V5M_DATA',
}

jecVersionsDATA = {
    'UL2016_preVFPB': 'Summer19UL16APV_RunBCD_V7_DATA',
    'UL2016_preVFPC': 'Summer19UL16APV_RunBCD_V7_DATA',
    'UL2016_preVFPD': 'Summer19UL16APV_RunBCD_V7_DATA',
    'UL2016_preVFPE': 'Summer19UL16APV_RunEF_V7_DATA',
    'UL2016_preVFPF': 'Summer19UL16APV_RunEF_V7_DATA',
    'UL2016F': 'Summer19UL16_RunFGH_V7_DATA',
    'UL2016G': 'Summer19UL16_RunFGH_V7_DATA',
    'UL2016H': 'Summer19UL16_RunFGH_V7_DATA',
    'UL2017B': 'Summer19UL17_RunB_V5_DATA',
    'UL2017C': 'Summer19UL17_RunC_V5_DATA',
    'UL2017D': 'Summer19UL17_RunD_V5_DATA',
    'UL2017E': 'Summer19UL17_RunE_V5_DATA',
    'UL2017F': 'Summer19UL17_RunF_V5_DATA',
    'UL2018A': 'Summer19UL18_RunA_V5_DATA',
    'UL2018B': 'Summer19UL18_RunB_V5_DATA',
    'UL2018C': 'Summer19UL18_RunC_V5_DATA',
    'UL2018D': 'Summer19UL18_RunD_V5_DATA',
    '2022C'  : 'Summer22_22Sep2023_RunCD_V2_DATA',
    '2022D'  : 'Summer22_22Sep2023_RunCD_V2_DATA',
    '2022_EEE': 'Summer22EE_22Sep2023_RunE_V2_DATA',
    '2022_EEF': 'Summer22EE_22Sep2023_RunF_V2_DATA',
    '2022_EEG': 'Summer22EE_22Sep2023_RunG_V2_DATA',
    '2023Cv123' : 'Summer23Prompt23_RunCv123_V1_DATA',
    '2023Cv4'   : 'Summer23Prompt23_RunCv4_V1_DATA',
    '2023_BPixD': 'Summer23BPixPrompt23_RunD_V1_DATA',
    '2024B':   'Winter24Prompt24_RunBCD_V5M_DATA',
    '2024C':   'Winter24Prompt24_RunBCD_V5M_DATA',
    '2024D':   'Winter24Prompt24_RunBCD_V5M_DATA',
    '2024Ev1': 'Winter24Prompt24_RunE_V5M_DATA',
    '2024Ev2': 'Winter24Prompt24_RunE_V5M_DATA',
    '2024F':   'Winter24Prompt24_RunF_V5M_DATA',
}

# https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
jerVersionsMC = {
    'UL2016_preVFP': 'Summer20UL16APV_JRV3_MC',
    'UL2016'   : 'Summer20UL16_JRV3_MC',
    'UL2017'   : 'Summer19UL17_JRV2_MC',
    'UL2018'   : 'Summer19UL18_JRV2_MC',
    '2022'     : 'Summer22_22Sep2023_JRV1_MC',
    '2022_EE'  : 'Summer22EE_22Sep2023_JRV1_MC',
    '2023'     : 'Summer23Prompt23_RunCv1234_JRV1_MC',
    '2023_BPix': 'Summer23BPixPrompt23_RunD_JRV1_MC',
    '2024'     : 'Summer23BPixPrompt23_RunD_JRV1_MC',#TEMP Must not use it
}

# https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
jerAK8VersionsMC = {
    'UL2016_preVFP': jerVersionsMC['UL2016_preVFP'],
    'UL2016': jerVersionsMC['UL2016'],
    'UL2017': 'Summer19UL17_JRV3_MC', # Cannot use JRV2. Missing JER stuff.
    'UL2018': jerVersionsMC['UL2018'],
    '2022': jerVersionsMC['2022'],
    '2022_EE': jerVersionsMC['2022_EE'],
    '2023': jerVersionsMC['2023'],
    '2023_BPix': jerVersionsMC['2023_BPix'],
    '2024': jerVersionsMC['2023_BPix'],
}

# jet mass resolution: https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging
#nominal, up, down
jmrValues = {
    '2016': [1.0, 1.2, 0.8],
    '2017': [1.09, 1.14, 1.04],
    # Use 2017 values for 2018 until 2018 are released
    '2018': [1.09, 1.14, 1.04],
    'UL2016_preVFP': [1.00, 1.00, 1.00],  # placeholder
    'UL2016': [1.00, 1.00, 1.00],  # placeholder
    'UL2017': [1.00, 1.00, 1.00],  # placeholder
    'UL2018': [1.00, 1.00, 1.00],  # placeholder
    '2022': [1.00, 1.00, 1.00],  # placeholder
    '2022_EE': [1.00, 1.00, 1.00],  # placeholder
    '2023': [1.00, 1.00, 1.00],  # placeholder
    '2023_BPix': [1.00, 1.00, 1.00],  # placeholder
    '2024': [1.00, 1.00, 1.00],  # placeholder
}

# jet mass scale
# W-tagging PUPPI softdrop JMS values: https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging
# 2016 values
jmsValues = {
    '2016': [1.00, 0.9906, 1.0094],  # nominal, down, up
    '2017': [0.982, 0.978, 0.986],
    # Use 2017 values for 2018 until 2018 are released
    '2018': [0.982, 0.978, 0.986],
    'UL2016_preVFP': [1.000, 1.000, 1.000],  # placeholder
    'UL2016': [1.000, 1.000, 1.000],  # placeholder
    'UL2017': [1.000, 1.000, 1.000],  # placeholder
    'UL2018': [1.000, 1.000, 1.000],  # placeholder
    '2022': [1.000, 1.000, 1.000],  # placeholder
    '2022_EE': [1.000, 1.000, 1.000],  # placeholder
    '2023': [1.000, 1.000, 1.000],  # placeholder
    '2023_BPix': [1.000, 1.000, 1.000],  # placeholder
    '2024': [1.000, 1.000, 1.000],  # placeholder
}

def createJMECorrector(isMC=True,
                       dataYear=2016,
                       runPeriod="B",
                       jesUncert="Total",
                       jetType="AK4PFchs",
                       noGroom=False,
                       metBranchName="MET",
                       applySmearing=True,
                       isFastSim=False,
                       applyHEMfix=False,
                       saveMETUncs=['T1', 'T1Smear']):

    dataYear = str(dataYear)

    if isMC and not isFastSim:
        jecVersion_ = jecVersionsMC[dataYear]
    elif isMC and isFastSim:
        jecVersion_ = jecVersionsFastSim[dataYear]
    else:
        jecVersion_ = jecVersionsDATA[dataYear+runPeriod]

    jecUncertainties_ = [x for x in jesUncert.split(",")]
    jerVersion_ = jerVersionsMC[dataYear]
    jerAK8Version_ = jerAK8VersionsMC[dataYear]
    jmrValues_ = jmrValues[dataYear]
    jmsValues_ = jmsValues[dataYear]
    if dataYear in archiveTagsDATA:
        archiveTag_ = archiveTagsDATA[dataYear]
    elif dataYear+runPeriod in archiveTagsDATA:
        archiveTag_ = archiveTagsDATA[dataYear+runPeriod]
    met_ = metBranchName
    print(f'JEC : {jecVersion_} \t JER : {jerVersion_}')
    print(f'MET branch : {met_}')
    jmeCorrections = None
    # jme corrections
    if 'AK4' in jetType:
        if isMC:
            jmeCorrections = lambda: jetmetUncertaintiesProducer(
                era=dataYear,
                jecVersion=jecVersion_,
                jesUncertainties=jecUncertainties_,
                jerVersion=jerVersion_,
                jetType=jetType,
                metBranchName=met_,
                applySmearing=applySmearing,
                applyHEMfix=applyHEMfix,
                saveMETUncs=saveMETUncs)
        else:
            jmeCorrections = lambda: jetmetUncertaintiesProducer(
                era=dataYear,
                archive=archiveTag_,
                jecVersion=jecVersion_,
                jesUncertainties=jecUncertainties_,
                jetType=jetType,
                metBranchName=met_,
                isData=True)
    # no MET variations calculated
    elif 'AK8' in jetType:
        if isMC:
            jmeCorrections = lambda: fatJetUncertaintiesProducer(
                era=dataYear,
                jecVersion=jecVersion_,
                jesUncertainties=jecUncertainties_,
                jetType="AK8PFPuppi",
                jerVersion=jerAK8Version_,
                jmrVals=jmrValues_,
                jmsVals=jmsValues_,
                applySmearing=applySmearing,
                applyHEMfix=applyHEMfix)
        else:
            jmeCorrections = lambda: fatJetUncertaintiesProducer(
                era=dataYear,
                archive=archiveTag_,
                jecVersion=jecVersion_,
                jesUncertainties=jecUncertainties_,
                jetType="AK8PFPuppi",
                jmrVals=jmrValues_,
                jmsVals=jmsValues_,
                isData=True)

    return jmeCorrections