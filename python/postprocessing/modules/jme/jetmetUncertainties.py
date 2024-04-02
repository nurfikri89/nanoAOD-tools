from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetReCalibrator import JetReCalibrator
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetSmearer import JetSmearer
from PhysicsTools.NanoAODTools.postprocessing.modules.jmeV2.JECTool import JECTool
from PhysicsTools.NanoAODTools.postprocessing.modules.jmeV2.JERTool import JERTool

from PhysicsTools.NanoAODTools.postprocessing.tools import matchObjectCollection, matchObjectCollectionMultiple
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
import ROOT
import math
import os
import re
import tarfile
import tempfile
import shutil
import numpy as np
import itertools
ROOT.PyConfig.IgnoreCommandLineOptions = True

class jetmetUncertaintiesProducer(Module):
    def __init__(self,
        era,
        jecVersion,
        jesUncertainties=["Total"],
        archive=None,
        jetType="AK4PFchs",
        metBranchName="MET",
        jerVersion="",
        isData=False,
        applySmearing=False,
        applyHEMfix=False,
        saveMETUncs=['T1', 'T1Smear']
     ):
        self.era = era
        self.isData = isData
        # if set to true, Jet_pt_nom will have JER applied. not to be
        # switched on for data.
        self.applySmearing = applySmearing if not isData else False
        self.splitJERIDs = [""]  # "empty" ID for the overall JER
        self.metBranchName = metBranchName
        self.rhoBranchName = "fixedGridRhoFastjetAll"
        if "2022" in era or "2023" in era:
            self.rhoBranchName = "Rho_fixedGridRhoFastjetAll" # TEMP. Should be re-checked in the future

        # --------------------------------------------------------------------
        #
        # Setup JEC
        #
        # --------------------------------------------------------------------
        if "AK4" in jetType:
            self.jetBranchName = "Jet"
            self.genJetBranchName = "GenJet"
        else:
            raise ValueError(f"ERROR: Invalid jet type = '{jetType}'!")
        self.lenVar = f"n{self.jetBranchName}"

        # --------------------------------------------------------------------
        # CV: globalTag and jetType not yet used in the jet smearer, as there
        # is no consistent set of txt files for JES uncertainties and JER scale
        # factors and uncertainties yet
        # --------------------------------------------------------------------
        self.jesUncertainties = jesUncertainties

        # Calculate and save uncertainties on T1Smear MET if this flag is set
        # to True. Otherwise calculate and save uncertainties on T1 MET
        self.saveMETUncs = saveMETUncs

        #
        #
        self.useCorrLib = False
        # --------------------------------------------------------------------
        #
        # Setup JEC
        #
        # --------------------------------------------------------------------
        self.jecTool = None
        self.jetReCalibrator = None
        self.jetReCalibratorL1 = None
        #####################################################################################
        # 1) correction lib
        #####################################################################################
        if self.useCorrLib:
            self.jecTool = JECTool(self.era, self.isData, jetType, jecVersion, self.jesUncertainties if not isData else [])
        #####################################################################################
        # 2) txt-file based approach
        # load libraries for accessing JEC factors and uncertainties from txt files
        #####################################################################################
        else:
            for library in ["libCondFormatsJetMETObjects", "libPhysicsToolsNanoAODTools"]:
                if library not in ROOT.gSystem.GetLibraries():
                    print("Load Library '{}'".format(library.replace("lib", "")))
                    ROOT.gSystem.Load(library)


            # read jet energy scale (JES) uncertaintiesJetSmearer
            # (downloaded from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC )
            self.jesInputArchivePath = f"{os.environ['CMSSW_BASE']}/src/PhysicsTools/NanoAODTools/data/jme/"
            # Text files are now tarred so must extract first into temporary
            # directory (gets deleted during python memory management at
            # script exit)
            fileExt = "tgz"
            if "2022" in era or "2023" in era: fileExt = "tar.gz" # TEMP. Should be re-checked in the future
            self.jesArchive = tarfile.open(f"{self.jesInputArchivePath}{jecVersion}.{fileExt}", "r:gz") if not archive else tarfile.open(f"{self.jesInputArchivePath}{archive}.{fileExt}", "r:gz")
            self.jesInputFilePath = tempfile.mkdtemp()
            self.jesArchive.extractall(self.jesInputFilePath)

            # Define the jet recalibrator
            print(f"jesInputFilePath={self.jesInputFilePath}")
            self.jetReCalibrator = JetReCalibrator(jecVersion,jetType,True,self.jesInputFilePath,calculateSeparateCorrections=False,calculateType1METCorrection=False)
            # Define the recalibrator for level 1 corrections only
            self.jetReCalibratorL1 = JetReCalibrator(jecVersion,jetType,False,self.jesInputFilePath,calculateSeparateCorrections=True,calculateType1METCorrection=False,upToLevel=1)

            # to fully re-calculate type-1 MET the JEC that are currently
            # applied are also needed.
            if len(jesUncertainties) == 1 and jesUncertainties[0] == "Total":
                self.jesUncertaintyInputFileName = f"{jecVersion}_Uncertainty_{jetType}.txt"
            elif jesUncertainties[0] == "Merged" and not self.isData:
                self.jesUncertaintyInputFileName = f"Regrouped_{jecVersion}_UncertaintySources_{jetType}.txt"
            else:
                self.jesUncertaintyInputFileName = f"{jecVersion}_UncertaintySources_{jetType}.txt"

            # read all uncertainty source names from the loaded file
            if jesUncertainties[0] in ["All", "Merged"]:
                with open(f'{self.jesInputFilePath}/{self.jesUncertaintyInputFileName}') as f:
                    lines = f.read().split("\n")
                    sources = [x for x in lines if x.startswith("[") and x.endswith("]")]
                    sources = [x[1:-1] for x in sources]
                    self.jesUncertainties = sources
            if applyHEMfix:
                self.jesUncertainties.append("HEMIssue")

        # --------------------------------------------------------------------
        #
        # Setup JER
        #
        # --------------------------------------------------------------------
        self.jerTool = None
        self.jetSmearer = None
        if self.applySmearing:
            #####################################################################################
            # 1) correction lib
            #####################################################################################
            if self.useCorrLib:
                self.jerTool = JERTool(self.era, jetType, jerVersion)
            #####################################################################################
            # 2) txt-file based approach
            # load libraries for accessing JER scale factors and uncertainties from txt files
            #####################################################################################
            else:
                # smear jet pT to account for measured difference in JER between data
                # and simulation.
                if jerVersion != "":
                    self.jerInputFileName = f"{jerVersion}_PtResolution_{jetType}.txt"
                    self.jerUncertaintyInputFileName = f"{jerVersion}_SF_{jetType}.txt"
                else:
                    raise ValueError(f"!")
                self.jetSmearer = JetSmearer(era, jetType, self.jerInputFileName, self.jerUncertaintyInputFileName)

        # --------------------------------------------------------------------
        # define energy threshold below which jets are considered as "unclustered energy"
        # cf. JetMETCorrections/Type1MET/python/correctionTermsPfMetType1Type2_cff.py
        # --------------------------------------------------------------------
        self.unclEnThreshold = 15.

    def beginJob(self):
        if not self.useCorrLib:
            print(f"Loading jet energy scale (JES) uncertainties from file {os.path.join(self.jesInputFilePath,self.jesUncertaintyInputFileName)}")
            self.jesUncertainty = {}
            # implementation didn't seem to work for factorized JEC,
            # try again another way
            for jesUncertainty in self.jesUncertainties:
                jesUncertainty_label = jesUncertainty
                if jesUncertainty == "Total" and (len(self.jesUncertainties) == 1 or (len(self.jesUncertainties) == 2 and "HEMIssue" in self.jesUncertainties)):
                    jesUncertainty_label = ''
                if jesUncertainty != "HEMIssue":
                    pars = ROOT.JetCorrectorParameters(os.path.join(self.jesInputFilePath,self.jesUncertaintyInputFileName),jesUncertainty_label)
                    self.jesUncertainty[jesUncertainty] = ROOT.JetCorrectionUncertainty(pars)
        else:
            self.jecTool.beginJob()

        if self.applySmearing:
            if self.useCorrLib:
                self.jerTool.beginJob()
            else:
                self.jetSmearer.beginJob()

    def endJob(self):
        if self.applySmearing:
            if not self.useCorrLib:
                self.jetSmearer.endJob()
        if not self.useCorrLib:
            shutil.rmtree(self.jesInputFilePath)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(f"{self.jetBranchName}_pt_raw"  ,"F",lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_pt_nom"  ,"F",lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_mass_raw","F",lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_mass_nom","F",lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_corr_JEC","F",lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_corr_JER","F",lenVar=self.lenVar)
        self.out.branch(f"{self.metBranchName}_T1_pt" ,  "F")
        self.out.branch(f"{self.metBranchName}_T1_phi",  "F")

        if not self.isData:
            self.out.branch(f"{self.metBranchName}_T1Smear_pt",  "F")
            self.out.branch(f"{self.metBranchName}_T1Smear_phi", "F")

            for shift in itertools.chain(["Up", "Down"]):
                self.out.branch(f"{self.jetBranchName}_pt_jer{shift}", "F", lenVar=self.lenVar)
                self.out.branch(f"{self.jetBranchName}_mass_jer{shift}", "F", lenVar=self.lenVar)
                if 'T1' in self.saveMETUncs:
                    self.out.branch(f"{self.metBranchName}_T1_pt_jer{shift}", "F")
                    self.out.branch(f"{self.metBranchName}_T1_phi_jer{shift}", "F")
                if 'T1Smear' in self.saveMETUncs:
                    self.out.branch(f"{self.metBranchName}_T1Smear_pt_jer{shift}", "F")
                    self.out.branch(f"{self.metBranchName}_T1Smear_phi_jer{shift}", "F")

                for jesUncertainty in itertools.chain(self.jesUncertainties):
                    self.out.branch(f"{self.jetBranchName}_pt_jes{jesUncertainty}{shift}","F", lenVar=self.lenVar)
                    self.out.branch(f"{self.jetBranchName}_mass_jes{jesUncertainty}{shift}","F", lenVar=self.lenVar)
                    if 'T1' in self.saveMETUncs:
                        self.out.branch(f"{self.metBranchName}_T1_pt_jes{jesUncertainty}{shift}", "F")
                        self.out.branch(f"{self.metBranchName}_T1_phi_jes{jesUncertainty}{shift}", "F")
                    if 'T1Smear' in self.saveMETUncs:
                        self.out.branch(f"{self.metBranchName}_T1Smear_pt_jes{jesUncertainty}{shift}", "F")
                        self.out.branch(f"{self.metBranchName}_T1Smear_phi_jes{jesUncertainty}{shift}", "F")

                if 'T1' in self.saveMETUncs:
                    self.out.branch(f"{self.metBranchName}_T1_pt_unclustEn{shift}", "F")
                    self.out.branch(f"{self.metBranchName}_T1_phi_unclustEn{shift}", "F")
                if 'T1Smear' in self.saveMETUncs:
                    self.out.branch(f"{self.metBranchName}_T1Smear_pt_unclustEn{shift}", "F")
                    self.out.branch(f"{self.metBranchName}_T1Smear_phi_unclustEn{shift}", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail,
        go to next event)"""
        jets = Collection(event, self.jetBranchName)
        nJet = event.nJet
        lowPtJets = Collection(event, "CorrT1METJet")
        # to subtract out of the jets for proper type-1 MET corrections
        if not self.isData:
            genJets = Collection(event, self.genJetBranchName)

        # rho needed to retrieve JEC
        rho = getattr(event, self.rhoBranchName)

        # prepare the low pt jets (they don't have a rawFactor)
        for jet in itertools.chain(lowPtJets):
            jet.pt = jet.rawPt
            jet.rawFactor = 0
            jet.mass = 0
            # the following dummy values should be removed once the values
            # are kept in nanoAOD
            jet.neEmEF = 0
            jet.chEmEF = 0

        # Set the seed for the random number generator that is used
        # for the JER smearing
        if self.applySmearing:
            if self.useCorrLib:
                self.jerTool.setSeed(event, self.jetBranchName)
            else:
                self.jetSmearer.setSeed(event)

        ########################################################
        #
        # Setup for JEC and JER corrected for jet pt and mass
        #
        ########################################################
        jets_pt_raw = []
        jets_pt_jer = []
        jets_pt_nom = []

        jets_mass_raw = []
        jets_mass_nom = []

        jets_corr_JEC = []
        jets_corr_JER = []

        ###########################################################
        #
        # Setup JEC and JER uncertainties for jet pt and mass
        #
        ###########################################################
        jets_pt_jesUp = {}
        jets_pt_jesDown = {}

        jets_mass_jesUp = {}
        jets_mass_jesDown = {}

        for jesUncertainty in itertools.chain(self.jesUncertainties):
            jets_pt_jesUp[jesUncertainty] = []
            jets_pt_jesDown[jesUncertainty] = []
            jets_mass_jesUp[jesUncertainty] = []
            jets_mass_jesDown[jesUncertainty] = []

        jets_pt_jerUp = []
        jets_pt_jerDown = []
        jets_mass_jerUp = []
        jets_mass_jerDown = []

        ###################################
        #
        # MET
        #
        ###################################
        met = Object(event, self.metBranchName)
        if "Puppi" in self.metBranchName:
            rawmet = Object(event, "RawPuppiMET")
        else:
            rawmet = Object(event, "RawMET")

        (t1met_px, t1met_py) = (met.pt * math.cos(met.phi),met.pt * math.sin(met.phi))
        (met_px, met_py) = (rawmet.pt * math.cos(rawmet.phi),rawmet.pt * math.sin(rawmet.phi))
        (met_T1_px, met_T1_py) = (met_px, met_py)
        (met_T1Smear_px, met_T1Smear_py) = (met_px, met_py)

        if 'T1' in self.saveMETUncs:
            (met_T1_px_jerUp, met_T1_py_jerUp) = (met_px, met_py)
            (met_T1_px_jerDown, met_T1_py_jerDown) = (met_px, met_py)
        if 'T1Smear' in self.saveMETUncs:
            (met_T1Smear_px_jerUp, met_T1Smear_py_jerUp) = (met_px, met_py)
            (met_T1Smear_px_jerDown, met_T1Smear_py_jerDown) = (met_px, met_py)

        if 'T1' in self.saveMETUncs:
            (met_T1_px_jesUp, met_T1_py_jesUp) = ({}, {})
            (met_T1_px_jesDown, met_T1_py_jesDown) = ({}, {})
            for jesUncertainty in itertools.chain(self.jesUncertainties):
                met_T1_px_jesUp[jesUncertainty] = met_px
                met_T1_py_jesUp[jesUncertainty] = met_py
                met_T1_px_jesDown[jesUncertainty] = met_px
                met_T1_py_jesDown[jesUncertainty] = met_py
        if 'T1Smear' in self.saveMETUncs:
            (met_T1Smear_px_jesUp, met_T1Smear_py_jesUp)     = ({}, {})
            (met_T1Smear_px_jesDown, met_T1Smear_py_jesDown) = ({}, {})
            for jesUncertainty in itertools.chain(self.jesUncertainties):
                met_T1Smear_px_jesUp[jesUncertainty] = met_px
                met_T1Smear_py_jesUp[jesUncertainty] = met_py
                met_T1Smear_px_jesDown[jesUncertainty] = met_px
                met_T1Smear_py_jesDown[jesUncertainty] = met_py

        ###################################
        #
        # Things to do with JER
        #
        ###################################
        # match reconstructed jets to generator level ones
        # (needed to evaluate JER scale factors and uncertainties)
        def resolution_matching(jet, genjet):
            '''Helper function to match to gen based on pt difference'''
            if self.useCorrLib:
                resolution = self.jerTool.cset_jerPtReso.evaluate(jet.eta, jet.pt, rho)
            else:
                params = ROOT.PyJetParametersWrapper()
                params.setJetEta(jet.eta)
                params.setJetPt(jet.pt)
                params.setRho(rho)
                resolution = self.jetSmearer.jer.getResolution(params)

            return abs(jet.pt - genjet.pt) < 3 * resolution * jet.pt

        pairs = None
        if self.applySmearing:
            pairs = matchObjectCollection(jets,genJets,dRmax=0.2,presel=resolution_matching)
            lowPtPairs = matchObjectCollection(lowPtJets,genJets,dRmax=0.2,presel=resolution_matching)
            pairs.update(lowPtPairs)

        ###############################################
        #
        # Loop over "Jet" and "LowPtJets" collections
        #
        ###############################################
        for iJet, jet in enumerate(itertools.chain(jets, lowPtJets)):
            jet_pt, jet_mass = jet.pt, jet.mass
            jet_pt_nano = jet.pt
            jet_mass_nano = jet.mass
            rawFactor = jet.rawFactor
            jet_area = jet.area
            jet_eta = jet.eta
            jet_phi = jet.phi

            if hasattr(jet, "rawFactor"):
                jet_pt_raw = jet_pt_nano * (1 - jet.rawFactor)
                jet_mass_raw = jet_mass_nano * (1 - jet.rawFactor)
            else:
                jet_pt_raw = -1.0 * jet_pt_nano  # If factor not present factor will be saved as -1
                jet_mass_raw = -1.0 * jet_mass_nano  # If factor not present factor will be saved as -1

            #
            # Get the JEC-corrected jet_pt
            #
            ################################################################################
            # corrlib-tools
            ################################################################################
            if self.useCorrLib:
                jec = self.jecTool.getJECFactorL1L2L3Res(jet_area, jet_eta, jet_pt_raw, rho)
                jet_pt_corr, jet_mass_corr = (jet_pt_raw * jec), (jet_mass_raw * jec)
                jet.pt = jet_pt_corr
                jet.mass = jet_mass_corr

                jecL1 = self.jecTool.getJECFactorL1(jet_area, jet_eta, jet_pt_raw, rho)
                jet_pt_L1corr, jet_mass_L1corr = (jet_pt_raw * jecL1), (jet_mass_raw * jecL1)
            ################################################################################
            else:
                jet_pt_corr, jet_mass_corr = self.jetReCalibrator.correct(jet, rho)
                jet.pt = jet_pt_corr
                jet.mass = jet_mass_corr
                # Get the JEC factors
                jec = jet_pt_corr / jet_pt_raw

                jet_pt_L1corr, jet_mass_L1corr = self.jetReCalibratorL1.correct(jet, rho)
                jecL1 = jet_pt_L1corr / jet_pt_raw


            genJet = None
            if self.applySmearing:
                genJet = pairs[jet]

            # get the jet for type-1 MET
            jet_p4_raw_noMu = ROOT.TLorentzVector()
            jet_p4_raw_noMu.SetPtEtaPhiM(jet_pt_nano * (1 - jet.rawFactor) * (1 - jet.muonSubtrFactor), jet.eta, jet.phi, jet.mass)
            muon_pt = jet_pt_nano * (1 - jet.rawFactor) * jet.muonSubtrFactor

            # get the proper jet pts for type-1 MET
            jet_pt_raw_noMu = jet_p4_raw_noMu.Pt()
            jet_pt_noMuL1L2L3 = jet_pt_raw_noMu * jec
            jet_pt_noMuL1 = jet_pt_raw_noMu * jecL1

            # setting jet back to central values
            jet.pt = jet_pt_corr
            jet.mass = jet_mass_corr
            jet.rawFactor = rawFactor

            # evaluate JER scale factors and uncertainties
            # cf. https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution and
            # https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution
            if self.applySmearing:
                if self.useCorrLib:
                    ################################################################################
                    # corrlib-tools
                    ################################################################################
                    (jet_pt_jerNomVal, jet_pt_jerUpVal, jet_pt_jerDownVal) = self.jerTool.getSmearingFactorForJet(jet, genJet, rho)
                    ################################################################################
                else:
                    (jet_pt_jerNomVal, jet_pt_jerUpVal, jet_pt_jerDownVal) = self.jetSmearer.getSmearValsPt(jet, genJet, rho)
            else:
                # if you want to do something with JER in data, please add it here.
                (jet_pt_jerNomVal, jet_pt_jerUpVal, jet_pt_jerDownVal) = (1., 1., 1.)

            #
            # Set nominal pt to be jec-applied+smeared pt if asked to apply smearing. If not, just the
            # jec-applied pt.
            #
            jet_pt_nom = jet_pt_jerNomVal * jet_pt_corr if self.applySmearing else jet_pt_corr
            jet_pt_L1L2L3 = jet_pt_noMuL1L2L3 + muon_pt
            jet_pt_L1 = jet_pt_noMuL1 + muon_pt

            jet_mass_nom = jet_pt_jerNomVal * jet_mass_corr if self.applySmearing else jet_mass_corr
            if jet_mass_nom < 0.0:
                jet_mass_nom *= -1.0

            # don't store the low pt jets in the Jet_pt_nom branch
            if iJet < nJet:
                jets_pt_raw.append(jet_pt_raw)
                jets_pt_nom.append(jet_pt_nom)
                jets_mass_raw.append(jet_mass_raw)
                jets_mass_nom.append(jet_mass_nom)
                jets_corr_JEC.append(jec)
                # can be used to undo JER
                jets_corr_JER.append(jet_pt_jerNomVal)

            if not self.isData:
                ################################################
                #
                # evaluate JES uncertainties
                #
                ################################################
                jet_pt_jesUp = {}
                jet_pt_jesDown = {}

                jet_pt_jesUpT1 = {}
                jet_pt_jesDownT1 = {}

                jet_mass_jesUp = {}
                jet_mass_jesDown = {}

                for jesUncertainty in itertools.chain(self.jesUncertainties):
                    # cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorUncertainties
                    # cf. https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
                    if jesUncertainty == "HEMIssue":
                        delta = 1.
                        if iJet < nJet and jet_pt_nom > 15 and jet.jetId & 2 and jet.phi > -1.57 and jet.phi < -0.87:
                            if jet.eta > -2.5 and jet.eta < -1.3:
                                delta = 0.8
                            elif jet.eta <= -2.5 and jet.eta > -3:
                                delta = 0.65
                        jet_pt_jesUp[jesUncertainty]     = jet_pt_nom
                        jet_pt_jesDown[jesUncertainty]   = delta * jet_pt_nom
                        jet_mass_jesUp[jesUncertainty]   = jet_mass_nom
                        jet_mass_jesDown[jesUncertainty] = delta * jet_mass_nom

                        jet_pt_jesUpT1[jesUncertainty]   = jet_pt_L1L2L3
                        jet_pt_jesDownT1[jesUncertainty] = delta * jet_pt_L1L2L3

                    else:
                        if self.useCorrLib:
                            delta = self.jecTool.getJECUncertainty(jet_eta, jet_pt_nom, jesUncertainty)
                        else:
                            self.jesUncertainty[jesUncertainty].setJetPt(jet_pt_nom)
                            self.jesUncertainty[jesUncertainty].setJetEta(jet_eta)
                            delta = self.jesUncertainty[jesUncertainty].getUncertainty(True)

                        jet_pt_jesUp[jesUncertainty]     = jet_pt_nom * (1. + delta)
                        jet_pt_jesDown[jesUncertainty]   = jet_pt_nom * (1. - delta)
                        jet_mass_jesUp[jesUncertainty]   = jet_mass_nom * (1. + delta)
                        jet_mass_jesDown[jesUncertainty] = jet_mass_nom * (1. - delta)

                        # redo JES variations for T1 MET
                        if self.useCorrLib:
                            delta = self.jecTool.getJECUncertainty(jet_eta, jet_pt_L1L2L3, jesUncertainty)
                        else:
                            self.jesUncertainty[jesUncertainty].setJetPt(jet_pt_L1L2L3)
                            self.jesUncertainty[jesUncertainty].setJetEta(jet.eta)
                            delta = self.jesUncertainty[jesUncertainty].getUncertainty(True)

                        jet_pt_jesUpT1[jesUncertainty] = jet_pt_L1L2L3 * (1. + delta)
                        jet_pt_jesDownT1[jesUncertainty] = jet_pt_L1L2L3 * (1. - delta)

                    if iJet < nJet:
                        jets_pt_jesUp[jesUncertainty].append(jet_pt_jesUp[jesUncertainty])
                        jets_pt_jesDown[jesUncertainty].append(jet_pt_jesDown[jesUncertainty])
                        jets_mass_jesUp[jesUncertainty].append(jet_mass_jesUp[jesUncertainty])
                        jets_mass_jesDown[jesUncertainty].append(jet_mass_jesDown[jesUncertainty])

                ################################################
                #
                # evaluate JER uncertainties
                #
                ################################################
                jet_pt_jerUp     = jet_pt_jerUpVal   * jet_pt_corr
                jet_pt_jerDown   = jet_pt_jerDownVal * jet_pt_corr
                jet_mass_jerUp   = jet_pt_jerUpVal   * jet_mass_corr
                jet_mass_jerDown = jet_pt_jerDownVal * jet_mass_corr

                # don't store the low pt jets in the Jet_pt_nom branch
                if iJet < nJet:
                    jets_pt_jerUp.append(jet_pt_jerUp)
                    jets_pt_jerDown.append(jet_pt_jerDown)
                    jets_mass_jerUp.append(jet_mass_jerUp)
                    jets_mass_jerDown.append(jet_mass_jerDown)

            # progate JER and JES corrections and uncertainties to MET.
            # Only propagate JECs to MET if the corrected pt without the muon
            # is above the threshold
            if jet_pt_noMuL1L2L3 > self.unclEnThreshold and (jet.neEmEF + jet.chEmEF) < 0.9:
                # do not re-correct for jets that aren't included in METv2 recipe
                if not (self.metBranchName == 'METFixEE2017'and 2.65 < abs(jet.eta) < 3.14 and jet.pt *(1 - jet.rawFactor) < 50):
                    jet_cosPhi = math.cos(jet.phi)
                    jet_sinPhi = math.sin(jet.phi)
                    met_T1_px = met_T1_px - (jet_pt_L1L2L3 - jet_pt_L1) * jet_cosPhi
                    met_T1_py = met_T1_py - (jet_pt_L1L2L3 - jet_pt_L1) * jet_sinPhi
                    if not self.isData:
                        met_T1Smear_px = met_T1Smear_px - (jet_pt_L1L2L3 * jet_pt_jerNomVal - jet_pt_L1) * jet_cosPhi
                        met_T1Smear_py = met_T1Smear_py - (jet_pt_L1L2L3 * jet_pt_jerNomVal - jet_pt_L1) * jet_sinPhi
                        # Variations of T1 MET
                        if 'T1' in self.saveMETUncs:
                            # For uncertainties on T1 MET, the up/down
                            # variations are just the centrally smeared
                            # MET values
                            jerUpVal, jerDownVal = jet_pt_jerNomVal, jet_pt_jerNomVal
                            met_T1_px_jerUp = met_T1_px_jerUp - (jet_pt_L1L2L3 * jerUpVal - jet_pt_L1) * jet_cosPhi
                            met_T1_py_jerUp = met_T1_py_jerUp - (jet_pt_L1L2L3 * jerUpVal - jet_pt_L1) * jet_sinPhi
                            met_T1_px_jerDown = met_T1_px_jerDown - (jet_pt_L1L2L3 * jerDownVal - jet_pt_L1) * jet_cosPhi
                            met_T1_py_jerDown = met_T1_py_jerDown - (jet_pt_L1L2L3 * jerDownVal - jet_pt_L1) * jet_sinPhi
                            # Calculate JES uncertainties on unsmeared MET
                            for jesUncertainty in self.jesUncertainties:
                                met_T1_px_jesUp[jesUncertainty] = met_T1_px_jesUp[jesUncertainty] - (jet_pt_jesUpT1[jesUncertainty] - jet_pt_L1) * jet_cosPhi
                                met_T1_py_jesUp[jesUncertainty] = met_T1_py_jesUp[jesUncertainty] - (jet_pt_jesUpT1[jesUncertainty] - jet_pt_L1) * jet_sinPhi
                                met_T1_px_jesDown[jesUncertainty] = met_T1_px_jesDown[jesUncertainty] - (jet_pt_jesDownT1[jesUncertainty] - jet_pt_L1) * jet_cosPhi
                                met_T1_py_jesDown[jesUncertainty] = met_T1_py_jesDown[jesUncertainty] - (jet_pt_jesDownT1[jesUncertainty] - jet_pt_L1) * jet_sinPhi
                        # Variations of T1Smear MET
                        if 'T1Smear' in self.saveMETUncs:
                            jerUpVal, jerDownVal = jet_pt_jerNomVal, jet_pt_jerNomVal
                            met_T1Smear_px_jerUp = met_T1Smear_px_jerUp - (jet_pt_L1L2L3 * jerUpVal - jet_pt_L1) * jet_cosPhi
                            met_T1Smear_py_jerUp = met_T1Smear_py_jerUp - (jet_pt_L1L2L3 * jerUpVal - jet_pt_L1) * jet_sinPhi
                            met_T1Smear_px_jerDown = met_T1Smear_px_jerDown - (jet_pt_L1L2L3 * jerDownVal - jet_pt_L1) * jet_cosPhi
                            met_T1Smear_py_jerDown = met_T1Smear_py_jerDown - (jet_pt_L1L2L3 * jerDownVal - jet_pt_L1) * jet_sinPhi

                            # Calculate JES uncertainties on smeared MET
                            for jesUncertainty in self.jesUncertainties:
                                jesUp_correction_forT1SmearMET = (jet_pt_L1L2L3 * jet_pt_jerNomVal - jet_pt_L1) + (jet_pt_jesUpT1[jesUncertainty] - jet_pt_L1L2L3)
                                jesDown_correction_forT1SmearMET = (jet_pt_L1L2L3 * jet_pt_jerNomVal -jet_pt_L1) + (jet_pt_jesDownT1[jesUncertainty] - jet_pt_L1L2L3)
                                met_T1Smear_px_jesUp[jesUncertainty] = met_T1Smear_px_jesUp[jesUncertainty] - jesUp_correction_forT1SmearMET * jet_cosPhi
                                met_T1Smear_py_jesUp[jesUncertainty] = met_T1Smear_py_jesUp[jesUncertainty] - jesUp_correction_forT1SmearMET * jet_sinPhi
                                met_T1Smear_px_jesDown[jesUncertainty] = met_T1Smear_px_jesDown[jesUncertainty] - jesDown_correction_forT1SmearMET * jet_cosPhi
                                met_T1Smear_py_jesDown[jesUncertainty] = met_T1Smear_py_jesDown[jesUncertainty] - jesDown_correction_forT1SmearMET * jet_sinPhi

        if not self.isData:
            (met_T1_px_unclEnUp, met_T1_py_unclEnUp) = (met_T1_px, met_T1_py)
            (met_T1_px_unclEnDown, met_T1_py_unclEnDown) = (met_T1_px, met_T1_py)
            (met_T1Smear_px_unclEnUp, met_T1Smear_py_unclEnUp) = (met_T1Smear_px, met_T1Smear_py)
            (met_T1Smear_px_unclEnDown, met_T1Smear_py_unclEnDown) = (met_T1Smear_px, met_T1Smear_py)
            if "Puppi" in self.metBranchName:
                met_deltaPx_unclEn = 0.
                met_deltaPy_unclEn = 0.
            else:
                met_deltaPx_unclEn = getattr(event, f"{self.metBranchName}_MetUnclustEnUpDeltaX")
                met_deltaPy_unclEn = getattr(event, f"{self.metBranchName}_MetUnclustEnUpDeltaY")
            met_T1_px_unclEnUp = met_T1_px_unclEnUp + met_deltaPx_unclEn
            met_T1_py_unclEnUp = met_T1_py_unclEnUp + met_deltaPy_unclEn
            met_T1_px_unclEnDown = met_T1_px_unclEnDown - met_deltaPx_unclEn
            met_T1_py_unclEnDown = met_T1_py_unclEnDown - met_deltaPy_unclEn
            met_T1Smear_px_unclEnUp = met_T1Smear_px_unclEnUp + met_deltaPx_unclEn
            met_T1Smear_py_unclEnUp = met_T1Smear_py_unclEnUp + met_deltaPy_unclEn
            met_T1Smear_px_unclEnDown = met_T1Smear_px_unclEnDown - met_deltaPx_unclEn
            met_T1Smear_py_unclEnDown = met_T1Smear_py_unclEnDown - met_deltaPy_unclEn

        self.out.fillBranch(f"{self.jetBranchName}_pt_raw", jets_pt_raw)
        self.out.fillBranch(f"{self.jetBranchName}_pt_nom", jets_pt_nom)
        self.out.fillBranch(f"{self.jetBranchName}_corr_JEC", jets_corr_JEC)
        self.out.fillBranch(f"{self.jetBranchName}_corr_JER", jets_corr_JER)
        if not self.isData:
            self.out.fillBranch(f"{self.jetBranchName}_pt_jerUp", jets_pt_jerUp)
            self.out.fillBranch(f"{self.jetBranchName}_pt_jerDown", jets_pt_jerDown)

        self.out.fillBranch(f"{self.metBranchName}_T1_pt", math.sqrt(met_T1_px**2 + met_T1_py**2))
        self.out.fillBranch(f"{self.metBranchName}_T1_phi", math.atan2(met_T1_py, met_T1_px))
        self.out.fillBranch(f"{self.jetBranchName}_mass_raw", jets_mass_raw)
        self.out.fillBranch(f"{self.jetBranchName}_mass_nom", jets_mass_nom)

        if not self.isData:
            self.out.fillBranch(f"{self.jetBranchName}_mass_jerUp",  jets_mass_jerUp)
            self.out.fillBranch(f"{self.jetBranchName}_mass_jerDown",jets_mass_jerDown)

        if not self.isData:
            self.out.fillBranch(f"{self.metBranchName}_T1Smear_pt" , math.sqrt(met_T1Smear_px**2 + met_T1Smear_py**2))
            self.out.fillBranch(f"{self.metBranchName}_T1Smear_phi", math.atan2(met_T1Smear_py, met_T1Smear_px))

            if 'T1' in self.saveMETUncs:
                self.out.fillBranch(f"{self.metBranchName}_T1_pt_jerUp",math.sqrt(met_T1_px_jerUp**2 + met_T1_py_jerUp**2))
                self.out.fillBranch(f"{self.metBranchName}_T1_phi_jerUp",math.atan2(met_T1_py_jerUp, met_T1_px_jerUp))
                self.out.fillBranch(f"{self.metBranchName}_T1_pt_jerDown",math.sqrt(met_T1_px_jerDown**2 + met_T1_py_jerDown**2))
                self.out.fillBranch(f"{self.metBranchName}_T1_phi_jerDown", math.atan2(met_T1_py_jerDown, met_T1_px_jerDown))

            if 'T1Smear' in self.saveMETUncs:
                self.out.fillBranch(f"{self.metBranchName}_T1Smear_pt_jerUp",math.sqrt(met_T1Smear_px_jerUp**2 + met_T1Smear_py_jerUp**2))
                self.out.fillBranch(f"{self.metBranchName}_T1Smear_phi_jerUp",math.atan2(met_T1Smear_py_jerUp, met_T1Smear_px_jerUp))
                self.out.fillBranch(f"{self.metBranchName}_T1Smear_pt_jerDown",math.sqrt(met_T1Smear_px_jerDown**2 + met_T1Smear_py_jerDown**2))
                self.out.fillBranch(f"{self.metBranchName}_T1Smear_phi_jerDown",math.atan2(met_T1Smear_py_jerDown, met_T1Smear_px_jerDown))

            for jesUncertainty in itertools.chain(self.jesUncertainties):
                self.out.fillBranch(f"{self.jetBranchName}_pt_jes{jesUncertainty}Up",jets_pt_jesUp[jesUncertainty])
                self.out.fillBranch(f"{self.jetBranchName}_pt_jes{jesUncertainty}Down",jets_pt_jesDown[jesUncertainty])

                if 'T1' in self.saveMETUncs:
                    self.out.fillBranch(f"{self.metBranchName}_T1_pt_jes{jesUncertainty}Up",math.sqrt(met_T1_px_jesUp[jesUncertainty]**2 + met_T1_py_jesUp[jesUncertainty]**2))
                    self.out.fillBranch(f"{self.metBranchName}_T1_phi_jes{jesUncertainty}Up",math.atan2(met_T1_py_jesUp[jesUncertainty], met_T1_px_jesUp[jesUncertainty]))
                    self.out.fillBranch(f"{self.metBranchName}_T1_pt_jes{jesUncertainty}Down",math.sqrt(met_T1_px_jesDown[jesUncertainty]**2 + met_T1_py_jesDown[jesUncertainty]**2))
                    self.out.fillBranch(f"{self.metBranchName}_T1_phi_jes{jesUncertainty}Down",math.atan2(met_T1_py_jesDown[jesUncertainty],met_T1_px_jesDown[jesUncertainty]))

                if 'T1Smear' in self.saveMETUncs:
                    self.out.fillBranch(f"{self.metBranchName}_T1Smear_pt_jes{jesUncertainty}Up",math.sqrt(met_T1Smear_px_jesUp[jesUncertainty]**2 +met_T1Smear_py_jesUp[jesUncertainty]**2))
                    self.out.fillBranch(f"{self.metBranchName}_T1Smear_phi_jes{jesUncertainty}Up",math.atan2(met_T1Smear_py_jesUp[jesUncertainty], met_T1Smear_px_jesUp[jesUncertainty]))
                    self.out.fillBranch(f"{self.metBranchName}_T1Smear_pt_jes{jesUncertainty}Down", math.sqrt(met_T1Smear_px_jesDown[jesUncertainty]**2 + met_T1Smear_py_jesDown[jesUncertainty]**2))
                    self.out.fillBranch(f"{self.metBranchName}_T1Smear_phi_jes{jesUncertainty}Down", math.atan2(met_T1Smear_py_jesDown[jesUncertainty], met_T1Smear_px_jesDown[jesUncertainty]))

                self.out.fillBranch(f"{self.jetBranchName}_mass_jes{jesUncertainty}Up",  jets_mass_jesUp[jesUncertainty])
                self.out.fillBranch(f"{self.jetBranchName}_mass_jes{jesUncertainty}Down",jets_mass_jesDown[jesUncertainty])

            if 'T1' in self.saveMETUncs:
                self.out.fillBranch(f"{self.metBranchName}_T1_pt_unclustEnUp",math.sqrt(met_T1_px_unclEnUp**2 + met_T1_py_unclEnUp**2))
                self.out.fillBranch(f"{self.metBranchName}_T1_phi_unclustEnUp",math.atan2(met_T1_py_unclEnUp, met_T1_px_unclEnUp))
                self.out.fillBranch(f"{self.metBranchName}_T1_pt_unclustEnDown",math.sqrt(met_T1_px_unclEnDown**2 + met_T1_py_unclEnDown**2))
                self.out.fillBranch(f"{self.metBranchName}_T1_phi_unclustEnDown",math.atan2(met_T1_py_unclEnDown, met_T1_px_unclEnDown))
            if 'T1Smear' in self.saveMETUncs:
                self.out.fillBranch(f"{self.metBranchName}_T1Smear_pt_unclustEnUp",math.sqrt(met_T1Smear_px_unclEnUp**2 + met_T1Smear_py_unclEnUp**2))
                self.out.fillBranch(f"{self.metBranchName}_T1Smear_phi_unclustEnUp",math.atan2(met_T1Smear_py_unclEnUp, met_T1Smear_px_unclEnUp))
                self.out.fillBranch(f"{self.metBranchName}_T1Smear_pt_unclustEnDown",math.sqrt(met_T1Smear_px_unclEnDown**2 + met_T1Smear_py_unclEnDown**2))
                self.out.fillBranch(f"{self.metBranchName}_T1Smear_phi_unclustEnDown",math.atan2(met_T1Smear_py_unclEnDown, met_T1Smear_px_unclEnDown))

        return True

