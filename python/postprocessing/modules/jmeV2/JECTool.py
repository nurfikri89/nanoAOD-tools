from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
import math
import correctionlib._core as core
ROOT.PyConfig.IgnoreCommandLineOptions = True

class JECTool(Module):
  def __init__(self,
    era,
    isData,
    jetType,
    jecVersion,
    uncNames=[]):

    self.jetType = jetType
    self.jecVersion = jecVersion
    self.uncNames = uncNames

    corrLibErasDict = {
      "UL2016_preVFP": "2016preVFP_UL",
      "UL2016": "2016postVFP_UL",
      "UL2017": "2017_UL",
      "UL2018": "2018_UL"
    }
    if era in corrLibErasDict:
      self.corrLibEra = corrLibErasDict[era]
    else:
      raise ValueError(f"ERROR: era = {era} unrecognised. Please check")

    jsonJetNamesDict = {
      "AK4PFchs":   "jet",
      "AK4PFPuppi": "jet",
      "AK8PFPuppi": "fatjet",
    }
    if self.jetType in jsonJetNamesDict:
      self.jsonJetName = jsonJetNamesDict[self.jetType]
    else:
      raise ValueError(f"ERROR: jetType = {self.jetType} not in json. The options are AK4PFchs, AK4PFPuppi and AK8PFPuppi")

    # Load the correctionlib json file
    self.jsonName  = f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/"
    self.jsonName += f"{self.corrLibEra}/{self.jsonJetName}_jerc.json.gz"
    self.cset = core.CorrectionSet.from_file(self.jsonName)

    self.jecLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
    if isData: self.jecLevels += ['L2L3Residual']

    # Check JEC version must be in json file
    self.jecVersion = jecVersion
    if not(any(self.jecVersion in key for key in self.cset)):
      raise ValueError(f"ERROR: jecVersion = {self.jecVersion} not in json")

    # Check JEC uncertainties must be in json file
    if not(isData):
      for name in self.uncNames:
        if not(any(f"{self.jecVersion}_{name}_{self.jetType}" in key for key in self.cset)):
          raise ValueError(f"ERROR: jesUncName = {name} not in json!")
    else:
      self.uncNames = [] #If its data, this should be empty.

  def beginJob(self):
    print(f"JECTool::Loading jet energy scale (JES) corrections, uncertainties from file {self.jsonName}")

    # JEC
    self.jecFactors = {}
    print(f"JECTool::Loading compound JEC: {self.jecVersion}_L1L2L3Res_{self.jetType}")
    self.jecFactors["L1L2L3Res"] = self.cset.compound[f"{self.jecVersion}_L1L2L3Res_{self.jetType}"]

    # JEC by levels
    for level in self.jecLevels:
      print(f"JECTool::Loading JEC: {self.jecVersion}_{level}_{self.jetType}")
      self.jecFactors[level] = self.cset[f"{self.jecVersion}_{level}_{self.jetType}"]

    # JEC uncertainty
    self.jecUncertainties = {}
    for name in self.uncNames:
      print(f"JECTool::Loading JEC uncertainty: {self.jecVersion}_{name}_{self.jetType}")
      self.jecUncertainties[name] = self.cset[f"{self.jecVersion}_{name}_{self.jetType}"]

  def endJob(self):
    pass

  ################################################
  #
  # Numpy-array based arguments
  #
  ################################################
  def getJECFactorArray(self, jets_area, jets_eta, jets_pt, rhos, level):
    return self.jecFactors[level].evalv(jets_area, jets_eta, jets_pt, rhos)

  def getJECFactorArrayL1L2L3Res(self, jets_area, jets_eta, jets_pt, rhos):
    return self.getJECFactorArray(jets_area, jets_eta, jets_pt, rhos, "L1L2L3Res")

  def getJECFactorArrayL1(self, jets_area, jets_eta, jets_pt, rhos):
    return self.getJECFactorArray(jets_area, jets_eta, jets_pt, rhos, "L1FastJet")

  def getJECFactorsArrayByLevel(self, jets_area, jets_eta, jets_pt, rhos):
    jec = {}
    jet_pt = jet_pt_raw
    for level in self.jecLevels:
      jec[level] = self.getJECFactorArray(jets_area, jets_eta, jets_pt, rhos, level)
      jet_pt *= jec[level]
    return jec

  def getJECUncertaintyArray(self, jets_eta, jets_pt, name):
    return self.jecUncertainties[name].evalv(jets_eta, jets_pt)

  def getJECUncertaintiesArrayAll(self, jets_eta, jets_pt):
    uncertainties = {}
    for name in self.uncNames:
      uncertainties[name] = self.getJECUncertaintyArray(jets_eta, jets_pt, name)
    return uncertainties

  ################################################
  #
  # Per-jet arguments
  #
  ################################################
  def getJECFactor(self, jet_area, jet_eta, jet_pt, rho, level):
    return self.jecFactors[level].evaluate(jet_area, jet_eta, jet_pt, rho)

  def getJECFactorL1L2L3Res(self, jet_area, jet_eta, jet_pt_raw, rho):
    return self.getJECFactor(jet_area, jet_eta, jet_pt_raw, rho, "L1L2L3Res")

  def getJECFactorL1(self, jet_area, jet_eta, jet_pt_raw, rho):
    return self.getJECFactor(jet_area, jet_eta, jet_pt_raw, rho, "L1FastJet")

  def getJECFactorsByLevel(self, jet_area, jet_eta, jet_pt_raw, rho):
    jec = {}
    jet_pt = jet_pt_raw
    for level in self.jecLevels:
      jec[level] = self.getJECFactor(jet_area, jet_eta, jet_pt, rho, level)
      jet_pt *= jec[level]
    return jec

  def getJECUncertainty(self, jet_eta, jet_pt, name):
    return self.jecUncertainties[name].evaluate(jet_eta, jet_pt)

  def getJECUncertaintiesAll(self, jet_eta, jet_pt):
    uncertainties = {}
    for name in self.uncNames:
      uncertainties[name] = self.getJECUncertainty(jet_eta, jet_pt, name)
    return uncertainties
