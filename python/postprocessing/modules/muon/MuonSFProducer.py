from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
import os
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

import gzip
from correctionlib._core import CorrectionSet

class MuonSFProducer(Module):
    def __init__(self, era):
        self.era=era
        self.corrlibEra = ""
        if self.era == "UL2016APV":
          self.corrlibEra += "2016preVFP_UL"
        elif self.era == "UL2016":
          self.corrlibEra += "2016postVFP_UL"
        elif self.era == "UL2017":
          self.corrlibEra += "2017_UL"
        elif self.era == "UL2018":
          self.corrlibEra += "2018_UL"

        def GetCorrLib(jsonName):
            evaluator = None
            if jsonName.endswith(".json.gz"):
              with gzip.open(jsonName,'rt') as file:
                data = file.read().strip()
              evaluator = CorrectionSet.from_string(data)
            else:

              evaluator = CorrectionSet.from_file(jsonName)
            return evaluator

        #
        # muon corrections
        #
        muo_jsonName = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/"+self.corrlibEra+"/muon_Z.json.gz"
        print("muo_jsonName:"+muo_jsonName)
        self.corrlib_muon = GetCorrLib(muo_jsonName)

        #
        # SFs for reconstruction, ID and isolation efficiencies.
        #
        self.effSFTypes = [
            "NUM_GlobalMuons_DEN_genTracks",
            "NUM_TrackerMuons_DEN_genTracks",
            "NUM_SoftID_DEN_TrackerMuons",
            "NUM_SoftID_DEN_genTracks",
            "NUM_LooseID_DEN_TrackerMuons",
            "NUM_LooseID_DEN_genTracks",
            "NUM_MediumID_DEN_TrackerMuons",
            "NUM_MediumID_DEN_genTracks",
            "NUM_MediumPromptID_DEN_TrackerMuons",
            "NUM_MediumPromptID_DEN_genTracks",
            "NUM_HighPtID_DEN_TrackerMuons",
            "NUM_HighPtID_DEN_genTracks",
            "NUM_TightID_DEN_TrackerMuons",
            "NUM_TightID_DEN_genTracks",
            "NUM_TrkHighPtID_DEN_TrackerMuons",
            "NUM_TrkHighPtID_DEN_genTracks",
            "NUM_LooseRelIso_DEN_LooseID",
            "NUM_LooseRelIso_DEN_MediumID",
            "NUM_LooseRelIso_DEN_MediumPromptID",
            "NUM_LooseRelIso_DEN_TightIDandIPCut",
            "NUM_LooseRelTkIso_DEN_HighPtIDandIPCut",
            "NUM_LooseRelTkIso_DEN_TrkHighPtIDandIPCut",
            "NUM_TightRelIso_DEN_MediumID",
            "NUM_TightRelIso_DEN_MediumPromptID",
            "NUM_TightRelIso_DEN_TightIDandIPCut",
            "NUM_TightRelTkIso_DEN_HighPtIDandIPCut",
            "NUM_TightRelTkIso_DEN_TrkHighPtIDandIPCut",
        ]
        #
        # SFs for triggers
        #
        if self.era == "UL2016APV":
            self.effSFTypes += [
                "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight",
                "NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose",
            ]
        elif self.era == "UL2016":
            self.effSFTypes += [
                "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight",
                "NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose",
            ]
        elif self.era == "UL2017":
            self.effSFTypes += [
                "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight",
                "NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose",
            ]
        elif self.era == "UL2018":
            self.effSFTypes += [
                "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
                "NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose",
            ]
        #
        # Use this if you want to store ALL scale factors in the json file.
        #
        #### self.effSFTypes = []
        #### for name in self.corrlib_muon:
        ####     self.effSFTypes += [name]

        self.corrlibProvider = {}

        for sfType in self.effSFTypes:
            self.corrlibProvider["SF_Eff_"+sfType] = self.corrlib_muon[sfType]

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for sfType in self.effSFTypes:
            self.out.branch("Muon_SF_Eff_"+sfType,                  "F", lenVar="nMuon")
            self.out.branch("Muon_SF_Eff_"+sfType+"_systTotalUp",   "F", lenVar="nMuon")
            self.out.branch("Muon_SF_Eff_"+sfType+"_systTotalDown", "F", lenVar="nMuon")
 
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = Collection(event, "Muon")

        #
        #
        #
        for sfType in self.effSFTypes:
            minPt = 15.
            if "Mu50" in sfType: minPt = 52.
            if "IsoMu24" in sfType: minPt = 26.
            if "IsoMu27" in sfType: minPt = 29.
            setattr(self,"SF_Eff_"+sfType, 
                [self.corrlibProvider["SF_Eff_"+sfType].evaluate(self.corrlibEra, abs(mu.eta), mu.pt, "sf") if (mu.pt >= minPt and abs(mu.eta) <= 2.4) else 1.0 for mu in muons]
            )
            setattr(self,"SF_Eff_"+sfType+"_systTotalUp", 
                [self.corrlibProvider["SF_Eff_"+sfType].evaluate(self.corrlibEra, abs(mu.eta), mu.pt, "systup") if (mu.pt >= minPt and abs(mu.eta) <= 2.4) else 1.0 for mu in muons]
            )
            setattr(self,"SF_Eff_"+sfType+"_systTotalDown", 
                [self.corrlibProvider["SF_Eff_"+sfType].evaluate(self.corrlibEra, abs(mu.eta), mu.pt, "systdown") if (mu.pt >= minPt and abs(mu.eta) <= 2.4) else 1.0 for mu in muons]
            )

        for sfType in self.effSFTypes:
            self.out.fillBranch("Muon_SF_Eff_"+sfType,                  getattr(self,"SF_Eff_"+sfType))
            self.out.fillBranch("Muon_SF_Eff_"+sfType+"_systTotalUp",   getattr(self,"SF_Eff_"+sfType+"_systTotalUp"))
            self.out.fillBranch("Muon_SF_Eff_"+sfType+"_systTotalDown", getattr(self,"SF_Eff_"+sfType+"_systTotalDown"))

        return True

MuonProducer_UL2016APV = lambda: MuonSFProducer("UL2016APV")
MuonProducer_UL2016    = lambda: MuonSFProducer("UL2016")
MuonProducer_UL2017    = lambda: MuonSFProducer("UL2017")
MuonProducer_UL2018    = lambda: MuonSFProducer("UL2018")

