from PhysicsTools.NanoAODHelpers.postprocessing.tools import matchObjectCollection, matchObjectCollectionMultiple
from PhysicsTools.NanoAODHelpers.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODHelpers.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODHelpers.postprocessing.modules.jme.JECTool import JECTool
from PhysicsTools.NanoAODHelpers.postprocessing.modules.jme.JERTool import JERTool
import ROOT
import os
ROOT.PyConfig.IgnoreCommandLineOptions = True

import gzip
import correctionlib
import correctionlib._core as core

class JetEnergyCalibrator(Module):
  def __init__(self,
    era,
    isData,
    jetType,
    jetCollection,
    jecVersion,
    jecUncNames=["Total"],
    jerVersion="",
    getSmearingFactor=False,
    applySmearing=False):

    self.era = era
    self.isData = isData
    self.rhoBranchName = "fixedGridRhoFastjetAll"
    #
    # if set to true, Jet_pt_nom will have JER applied. not to be
    # switched on for data.
    self.applySmearing = applySmearing if not isData else False
    #
    #
    #
    self.jecUncNames = jecUncNames if not(isData) else []

    #
    self.isAK4Jet = False
    self.isLowPtJet = False
    self.isAK8Jet = False
    self.isSubJet = False
    if jetCollection == "Jet":
      self.jetBranchName = "Jet"
      self.genJetBranchName = "GenJet"
      self.jsonJetName = "jet"
      self.jetType = jetType
      self.isAK4Jet = True
    elif jetCollection == "CorrT1METJet":
      self.jetBranchName = "CorrT1METJet"
      self.genJetBranchName = "GenJet"
      self.jsonJetName = "jet"
      self.jetType = jetType
      self.isAK4Jet = True
      self.isLowPtJet = True
    elif jetCollection == "FatJet":
      self.jetBranchName = "FatJet"
      self.genJetBranchName = "GenJetAK8"
      self.jsonJetName = "fatJet"
      self.jetType = jetType
      self.isAK8Jet = True
    elif jetCollection == "SubJet":
      self.jetBranchName = "SubJet"
      self.genJetBranchName = "SubGenJet"
      self.jsonJetName = "jet"#UseAK4
      self.jetType = jetType #Should use AK4PFPuppi
      self.isSubJet = True
      if self.jetType != "AK4PFPuppi":
        raise ValueError(f"ERROR: Invalid jetType = {self.jetType} for SubJet collection. Should be AK4PFPuppi!")
    else:
        raise ValueError(f"ERROR: Invalid jetCollection = {jetCollection}!")
    self.lenVar = f"n{self.jetCollection}"

    #
    if era == "UL2016APV":
      self.corrLibEra = "2016preVFP_UL"
    elif era == "UL2016":
      self.corrLibEra = "2016postVFP_UL"
    elif era == "UL2017":
      self.corrLibEra = "2016postVFP_UL"
    elif era == "UL2018":
      self.corrLibEra = "2016postVFP_UL"


    # Load the correctionlib json file
    self.jsonName = f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/"
    self.jsonName += f"{self.corrlibEra}/{self.jsonJetName}_jerc.json.gz"
    self.cset = core.CorrectionSet.from_file(jsonName)

    #
    self.jecTool = None
    if not(self.isData) and getSmearingFactor:
      self.jecTool = JECTool(self.cset, self.isData, self.jecVersion, self.jetType)

    #
    self.jerTool = None
    self.getSmearingFactor = getEnergySmearing
    if not(self.isData) and getSmearingFactor:
      self.jerTool = JERTool(self.cset, self.jerVersion, self.jetType)

    def beginJob(self):
      self.jecTool.beginJob()
      if self.jerTool:
        self.jerTool.beginJob()

    def endJob(self):
      pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
      self.out = wrappedOutputTree
      self.out.branch(f"{self.jetBranchName}_pt_raw","F",lenVar=self.lenVar)
      self.out.branch(f"{self.jetBranchName}_pt_nom","F",lenVar=self.lenVar)
      self.out.branch(f"{self.jetBranchName}_mass_raw","F",lenVar=self.lenVar)
      self.out.branch(f"{self.jetBranchName}_mass_nom","F",lenVar=self.lenVar)
      self.out.branch(f"{self.jetBranchName}_corr_JEC","F",lenVar=self.lenVar)
      self.out.branch(f"{self.jetBranchName}_corr_JER","F",lenVar=self.lenVar)
      # self.out.branch(f"{self.jetBranchName}_corr_JECL1","F",lenVar=self.lenVar)
      # self.out.branch(f"{self.jetBranchName}_corr_JECL1L2","F",lenVar=self.lenVar)
      if self.isAK4Jet:
        self.out.branch(f"{self.jetBranchName}_pt_noMuRaw","F",lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_pt_noMuJECL1","F",lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_pt_noMuJECL1L2L3","F",lenVar=self.lenVar)

      if not(isData):
        for uncName in self.jecUncNames:
          self.out.branch(f"{self.jetBranchName}_pt_jes{uncName}Up","F",lenVar=self.lenVar)
          self.out.branch(f"{self.jetBranchName}_pt_jes{uncName}Down","F",lenVar=self.lenVar)
          self.out.branch(f"{self.jetBranchName}_mass_jes{uncName}Up","F",lenVar=self.lenVar)
          self.out.branch(f"{self.jetBranchName}_mass_jes{uncName}Down","F",lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_pt_jerUp","F",lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_pt_jerDown","F",lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_mass_jerUp","F",lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_mass_jerDown","F",lenVar=self.lenVar)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
      pass

    def analyze(self, event):
      """process event, return True (go to next module) or False (fail,
      go to next event)"""

      rho = getattr(event, self.rhoBranchName)
      jets = Collection(event, self.jetBranchName)
      nJet = getattr(event, f"n{self.jetBranchName}")

      # If we are correcting low pt jets,
      # fill some missing values
      if self.isLowPtJet:
        for iJet, jet in jets:
          jet.pt = jet.rawPt
          jet.mass = 0
          jet.rawFactor = 0

      if not self.isData:
        genJets = Collection(event, self.genJetBranchName)

      jets_pt_raw = []
      jets_pt_nom = []
      jets_mass_raw = []
      jets_mass_nom = []

      jets_corr_JEC = []
      jets_corr_JER = []

      # jets_corr_JECL1 = []
      # jets_corr_JECL1L2 = []

      jets_pt_jesUp = {}
      jets_pt_jesDown = {}
      jets_mass_jesUp = {}
      jets_mass_jesDown = {}

      for uncName in self.jecUncNames:
        jets_pt_jesUp[uncName] = []
        jets_pt_jesDown[uncName] = []
        jets_mass_jesUp[uncName] = []
        jets_mass_jesDown[uncName] = []

      jets_pt_jerUp = []
      jets_pt_jerDown = []
      jets_mass_jerUp = []
      jets_mass_jerDown = []

      jets_pt_noMuRaw = []
      jets_pt_noMuJECL1 = []
      jets_pt_noMuJECL1L2L3 = []

      # match reconstructed jets to generator level ones
      # (needed to evaluate JER scale factors and uncertainties)
      def resolution_matching(jet, genjet):
        '''Helper function to match to gen based on pt difference'''
        resolution = self.jerPtReso.eval(jet.eta, jet.pt, rho)
        return abs(jet.pt - genjet.pt) < 3 * resolution * jet.pt

      recoToGenMatch = None
      if not self.isData:
        recoToGenMatch = matchObjectCollection(jets,genJets,dRmax=0.2,presel=resolution_matching)

      for iJet, jet in enumerate(jets):
        jet_pt = jet.pt
        jet_mass = jet.mass
        jet_eta = jet.eta
        jet_area = jet.area
        jet_pt_orig = jet_pt
        jet_mass_orig = jet_mass
        rawFactor = jet.rawFactor if hasattr(jet, "rawFactor") else 0
        ###################################
        #
        # Get raw pt and raw mass
        #
        ###################################
        if hasattr(jet, "rawFactor"):
          jet_pt_raw = jet_pt * (1 - rawFactor)
          jet_mass_raw = jet_mass * (1 - rawFactor)
        else:
          # If factor not present factor will be saved as -1
          jet_pt_raw = -1.0 * jet_pt
          jet_mass_raw = -1.0 * jet_mass

        ###################################
        #
        # Get JEC factor and apply it
        #
        ###################################
        jecFactor = self.jecTool.getJECFactorL1L2L3Res(jet_area, jet_eta, jet_pt_raw, rho)
        jet_pt_jec = jet_pt_raw * jecFactor
        jet_mass_jec = jet_mass_raw * jecFactor

        jecFactorsByLevel = self.jecTool.getJECFactorsByLevel(jet_area, jet_eta, jet_pt_raw, rho)
        ###################################
        #
        # Get smearing factor for JER
        #
        ###################################
        smearFactor = smearFactorUp = smearFactorDown = 1.

        # Retrieve the smearing factor for this jet
        if getSmearingFactor:
          genJet = recoToGenMatch[jet]
          genjet_pt = genJet.pt if genJet else -1.
          self.jerTool.setSeed(event, self.jetBranchName)
          smearFactor, smearFactorUp, smearFactorDown = \
            self.jerTool.getSmearingFactor(jet_pt_jec, jet_eta, rho, genjet_pt)

        jet_pt_nom = jet_pt_jec * smearFactor if self.applySmearing else jet_pt_jec
        jet_mass_nom = jet_mass_jec * smearFactor  if self.applySmearing else jet_mass_jec

        if jet_mass_nom < 0.0:
          jet_mass_nom *= -1.0

        jets_pt_raw.append(jet_pt_raw)
        jets_pt_nom.append(jet_pt_nom)
        jets_mass_raw.append(jet_mass_raw)
        jets_mass_nom.append(jet_mass_nom)
        jets_corr_JEC.append(jecFactor)
        jets_corr_JER.append(smearFactor)

        if not self.isData:
          ###################################
          #
          # evaluate JES uncertainties
          #
          ###################################
          jet_pt_jesUp = {}
          jet_pt_jesDown = {}
          jet_mass_jesUp = {}
          jet_mass_jesDown = {}
          for uncName in self.jecUncNames:
            delta = self.jecTool.getJECUncertainty.evaluate(jet_eta, jet_pt_nom, uncName)
            jet_pt_jesUp[uncName] = jet_pt_nom * (1. + delta)
            jet_pt_jesDown[uncName] = jet_pt_nom * (1. - delta)
            jet_mass_jesUp[uncName] = jet_mass_nom * (1. + delta)
            jet_mass_jesDown[uncName] = jet_mass_nom * (1. - delta)

            jets_pt_jesUp[uncName].append(jet_pt_jesUp[uncName])
            jets_pt_jesDown[uncName].append(jet_pt_jesDown[uncName])
            jets_mass_jesUp[uncName].append(jet_mass_jesUp[uncName])
            jets_mass_jesDown[uncName].append(jet_mass_jesDown[uncName])
          ###################################
          #
          # evaluate JER uncertainties
          #
          ###################################
          if getSmearingFactor:
            jet_pt_jerUp = jet_pt_nom * smearFactorUp if self.applySmearing else jet_pt_nom
            jet_pt_jerDown = jet_pt_nom * smearFactorDown if self.applySmearing else jet_pt_nom
            jet_mass_jerUp = jet_mass_nom * smearFactorUp if self.applySmearing else jet_mass_nom
            jet_mass_jerDown = jet_mass_nom * smearFactorDown if self.applySmearing else jet_mass_nom


        #
        #
        #
        if self.isAK4Jet and not(self.isLowPtJet):
          jet_pt_raw_noMu = jet_pt_orig * (1 - jet.rawFactor) * (1 - jet.muonSubtrFactor)
          muon_pt = jet_pt_orig * (1 - jet.rawFactor) * jet.muonSubtrFactor
          newjet = ROOT.TLorentzVector()
          newjet.SetPtEtaPhiM(jet_pt_raw_noMu, jet.eta, jet.phi, jet.mass)

          jet_pt_noMuRaw = newjet.Pt()
          jecFactorL1 = self.jecTool.getJECFactorL1(jet_area, jet_eta, jet_pt_noMuRaw, rho) # To CHECK: Should it be the original jet raw pt?
          jet_pt_noMuJECL1 = jet_pt_noMuRaw * jecFactorL1

        #
        #
        #
        ###################################
        #
        # Jet Loop end
        #
        ###################################

      ###################################
      #
      # Fill Branches
      #
      ###################################
      self.out.fillBranch(f"{self.jetBranchName}_pt_raw",jets_pt_raw)
      self.out.fillBranch(f"{self.jetBranchName}_pt_nom",jets_pt_nom)
      self.out.fillBranch(f"{self.jetBranchName}_mass_raw",jets_mass_raw)
      self.out.fillBranch(f"{self.jetBranchName}_mass_nom",jets_mass_nom)
      self.out.fillBranch(f"{self.jetBranchName}_corr_JEC",jets_corr_JEC)
      self.out.fillBranch(f"{self.jetBranchName}_corr_JER",jets_corr_JER)
      # self.out.fillBranch(f"{self.jetBranchName}_corr_JECL1",jets_corr_JECL1)
      # self.out.fillBranch(f"{self.jetBranchName}_corr_JECL1L2",jets_corr_JECL1L2)
      for uncName in self.jecUncNames:
        self.out.fillBranch(f"{self.jetBranchName}_pt_jes{uncName}Up",jets_pt_jesUp[uncName])
        self.out.fillBranch(f"{self.jetBranchName}_pt_jes{uncName}Down",jets_pt_jesDown[uncName])
        self.out.fillBranch(f"{self.jetBranchName}_mass_jes{uncName}Up",jets_mass_jesUp[uncName])
        self.out.fillBranch(f"{self.jetBranchName}_mass_jes{uncName}Down",jets_mass_jesDown[uncName])
      return True
