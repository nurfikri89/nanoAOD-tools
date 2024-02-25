from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import ROOT
import os
from itertools import chain
import gzip

from correctionlib._core import CorrectionSet

ROOT.PyConfig.IgnoreCommandLineOptions = True

class BTagSFProducerCorrLib(Module):

    def __init__(
        self, era, algo='deepJet', selectedWPs=['L','M','T']
    ):
        self.era = era
        self.algo = algo
        self.selectedWPs = selectedWPs

        self.corrlibEra = ""
        if self.era == "UL2016APV": 
          self.corrlibEra += "2016preVFP_UL"
        elif self.era == "UL2016":
          self.corrlibEra += "2016postVFP_UL"
        elif self.era == "UL2017":
          self.corrlibEra += "2017_UL"
        elif self.era == "UL2018":
          self.corrlibEra += "2018_UL"

        self.max_abs_eta = 2.5

        self.jet_flav_names = {
            5: "bjet",
            4: "cjet",
            0: "ljet",
        }
        
        #
        # NOTE
        # - tagger name must follow exactly as in json file.
        #
        supported_btagSF = {
            'deepJet': { 
                'measurement_types': {
                    5: "comb",  # b
                    4: "comb",  # c
                    0: "incl"   # light
                },
                'supported_wp': ["L", "M", "T"]
            },
        }
        self.measurement_types = supported_btagSF[self.algo]['measurement_types']
        self.supported_wp = supported_btagSF[self.algo]['supported_wp']

        # 
        # Define systematic uncertainties 
        #
        self.systs = []
        self.systs.append("up")
        self.systs.append("down")
        self.central_and_systs = ["central"]
        self.central_and_systs.extend(self.systs)

        self.branchNames_central_and_systs = {}
        for wp in self.selectedWPs:
            branchNames = {}
            central_and_systs = self.central_and_systs
            baseBranchName = 'Jet_btagSF_{}_{}'.format(self.algo, wp)
            for central_or_syst in central_and_systs:
                if central_or_syst == "central":
                    branchNames[central_or_syst] = baseBranchName
                else:
                    branchNames[central_or_syst] = baseBranchName+'_'+central_or_syst
            self.branchNames_central_and_systs[wp] = branchNames

    def beginJob(self):
        btv_jsonName = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/"+self.corrlibEra+"/btagging.json.gz"
        print("btv_jsonName:"+btv_jsonName)

        def GetCorrLib(jsonName):
            evaluator = None
            if jsonName.endswith(".json.gz"):
              with gzip.open(jsonName,'rt') as file:
                data = file.read().strip()
              evaluator = CorrectionSet.from_string(data)
            else:
              evaluator = CorrectionSet.from_file(jsonName)
            return evaluator

        self.corrlib_btv = GetCorrLib(btv_jsonName)

        self.corrlib_sf = {}
        for key in self.measurement_types:
            self.corrlib_sf[self.algo+"_"+self.jet_flav_names[key]] = self.corrlib_btv[self.algo+"_"+self.measurement_types[key]]

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for central_or_syst in list(self.branchNames_central_and_systs.values()):
            for branch in list(central_or_syst.values()):
                self.out.branch(branch, "F", lenVar="nJet")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def getSFs(self, jet_data, wp, syst):
        scale_factors = []
        for idx, (pt, eta, flav) in enumerate(jet_data):
            if abs(eta) < self.max_abs_eta:
                sf = self.corrlib_sf[self.algo+"_"+self.jet_flav_names[flav]].evaluate(syst, wp, flav, abs(eta), pt)
            else:
                sf = 1.0
            # check if SF is OK. 
            if sf < 0.01:
                sf = 1.0
            scale_factors.append(sf)

        return scale_factors

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        jets = Collection(event, "Jet")

        preloaded_jets = [
            (jet.pt_nom if hasattr(jet,"pt_nom") else jet.pt , jet.eta, jet.hadronFlavour) for jet in jets
        ] # use pt_nom from nanoAOD-tools JEC corrections. If not available, use pt in NanoAOD.

        for wp in self.selectedWPs:
            for central_or_syst in self.central_and_systs:
                scale_factors = list(self.getSFs(preloaded_jets, wp, central_or_syst))
                self.out.fillBranch(self.branchNames_central_and_systs[wp][central_or_syst], scale_factors)
        return True


algo="deepJet"
selectedWPs=["L","M","T"]
BTagSF_DeepJet_UL2016APV = lambda: BTagSFProducerCorrLib("UL2016APV",algo=algo,selectedWPs=selectedWPs)
BTagSF_DeepJet_UL2016 = lambda: BTagSFProducerCorrLib("UL2016",algo=algo,selectedWPs=selectedWPs)
BTagSF_DeepJet_UL2017 = lambda: BTagSFProducerCorrLib("UL2017",algo=algo,selectedWPs=selectedWPs)
BTagSF_DeepJet_UL2018 = lambda: BTagSFProducerCorrLib("UL2018",algo=algo,selectedWPs=selectedWPs)

