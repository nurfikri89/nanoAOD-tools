from PhysicsTools.NanoAODTools.postprocessing.tools import matchObjectCollection, matchObjectCollectionMultiple
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
import ROOT
import math
import correctionlib._core as core
ROOT.PyConfig.IgnoreCommandLineOptions = True

class JERTool(Module):
  def __init__(self,
    era,
    jetType,
    jerVersion):

    self.jerVersion = jerVersion
    self.jetType = jetType

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

    self.rnd = ROOT.TRandom3(12345)

    #
    # Check JER version must be in json file
    #
    self.jerVersion = jerVersion
    hasJERSF = any(f"{self.jerVersion}_ScaleFactor_{self.jetType}" in key for key in self.cset)
    hasJetPtReso = any(f"{self.jerVersion}_PtResolution_{self.jetType}" in key for key in self.cset)
    if not(hasJERSF and hasJetPtReso):
      raise ValueError(f"ERROR: jerVersion = {self.jerVersion} not in json!")

  def beginJob(self):
    print(f"Loading JER: {self.jerVersion}_ScaleFactor_{self.jetType}")
    self.cset_jerSF = self.cset[f"{self.jerVersion}_ScaleFactor_{self.jetType}"]
    print(f"Loading Pt Reso for JER: {self.jerVersion}_PtResolution_{self.jetType}")
    self.cset_jerPtReso = self.cset[f"{self.jerVersion}_PtResolution_{self.jetType}"]

  def endJob(self):
    pass

  def setSeed(self, event, jetBranchName):
    """Set seed deterministically."""
    # (cf. https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h)
    runnum = int(event.run) << 20
    luminum = int(event.luminosityBlock) << 10
    evtnum = event.event
    jet0eta = int(getattr(event,f"{jetBranchName}_eta")[0] / 0.01 if getattr(event,f"n{jetBranchName}") > 0 else 0)
    seed = 1 + runnum + evtnum + luminum + jet0eta
    self.rnd.SetSeed(seed)

  def getSmearingFactorForJet(self, jetIn, genJetIn, rho):
    return self.getSmearingFactor(jetIn.pt, jetIn.eta, jetIn.phi, jetIn.mass, rho, genJetIn.pt if genJetIn else -1.)

  def getSmearingFactor(self, jet_pt, jet_eta, jet_phi, jet_mass, rho, genjet_pt=-1.):
    hasGenJet = genjet_pt > 0.
    jet_p4 = ROOT.TLorentzVector()
    jet_p4.SetPtEtaPhiM(jet_pt, jet_eta, jet_phi, jet_mass)
    # --------------------------------------------------------------------------------------------
    # CV: Smear jet pT to account for measured difference in JER between data and simulation.
    #     The function computes the nominal smeared jet pT simultaneously with the JER up and down shifts,
    #     in order to use the same random number to smear all three (for consistency reasons).
    #
    #     The implementation of this function follows PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
    #
    # --------------------------------------------------------------------------------------------
    variations = ["nom", "up", "down"]
    jerSFAndVariations= {}
    smearFactorAndVariations = {}

    for var in variations:
      jerSFAndVariations[var] = self.cset_jerSF.evaluate(jet_eta, var)
      smearFactorAndVariations[var] = 1.

    if hasGenJet:
      for var in variations:
        #
        # Case 1: we have a "good" generator level jet matched to
        # the reconstructed jet
        #
        dPt = jet_pt - genjet_pt
        smearFactorAndVariations[var] = 1.+(jerSFAndVariations[var] - 1.) * dPt / jet_pt
    else:
      jet_pt_resolution = self.cset_jerPtReso.evaluate(jet_eta, jet_pt, rho)

      randNumber = self.rnd.Gaus(0, jet_pt_resolution)
      for var in variations:
        if jerSFAndVariations[var] > 1.:
          #
          # Case 2: we don't have a generator level jet. Smear jet
          # pT using a random Gaussian variation
          #
          smearFactor = 1. + randNumber * math.sqrt(jerSFAndVariations[var]**2 - 1.)
        else:
          #
          # Case 3: we cannot smear this jet, as we don't have a
          # generator level jet and the resolution in data is better
          # than the resolution in the simulation, so we would need
          # to randomly "unsmear" the jet, which is impossible
          #
          smearFactor = 1.
        smearFactorAndVariations[var] = smearFactor

    for var in variations:
      # check that smeared jet energy remains positive,
      # as the direction of the jet would change ("flip")
      # otherwise - and this is not what we want
      if smearFactorAndVariations[var] * jet_p4.E() < 1.e-2:
        smearFactorAndVariations[var] = 1.e-2 / jet_p4.E()

    return (smearFactorAndVariations["nom"], smearFactorAndVariations["nom"], smearFactorAndVariations["down"])
