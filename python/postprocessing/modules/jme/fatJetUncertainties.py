from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetReCalibrator import JetReCalibrator
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetSmearer import JetSmearer
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
ROOT.PyConfig.IgnoreCommandLineOptions = True

class fatJetUncertaintiesProducer(Module):
    def __init__(
            self,
            era,
            jecVersion,
            jesUncertainties=["Total"],
            archive=None,
            jetType="AK8PFPuppi",
            jerVersion="",
            jmrVals=[],
            jmsVals=[],
            isData=False,
            applySmearing=False,
            applyHEMfix=False
    ):
        self.era = era
        self.isData = isData
        self.applySmearing = applySmearing if not isData else False  # don't smear for data
        # ---------------------------------------------------------------------
        # CV: jecVersion and jetType not yet used in the jet smearer, as there
        # is no consistent set of txt files for JES uncertainties and JER scale
        # factors and uncertainties yet
        # ---------------------------------------------------------------------
        if "AK8" in jetType:
            self.jetBranchName = "FatJet"
            self.subJetBranchName = "SubJet"
            self.genJetBranchName = "GenJetAK8"
            self.genSubJetBranchName = "SubGenJetAK8"
        else:
            raise ValueError(f"ERROR: Invalid jet type = '{jetType}'!")
        self.rhoBranchName = "fixedGridRhoFastjetAll"
        if "2022" in era or "2023" in era:
            self.rhoBranchName = "Rho_fixedGridRhoFastjetAll" # TEMP. Should be re-checked in the future
        self.lenVar = "n" + self.jetBranchName
        self.lenVarSubjets = "n" + self.subJetBranchName

        #############################################################################
        #
        # AK8 Puppi JEC & JER
        #
        #############################################################################
        self.jesUncertainties = jesUncertainties
        # read jet energy scale (JES) uncertainties
        # (downloaded from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC )
        self.jesInputArchivePath = f"{os.environ['CMSSW_BASE']}/src/PhysicsTools/NanoAODTools/data/jme/"
        # Text files are now tarred so must extract first into temporary
        # directory (gets deleted during python memory management at
        # script exit)
        fileExt = "tgz"
        if "2022" in era or "2023" in era: fileExt = "tar.gz" # TEMP. Should be re-checked in the future
        if not archive:
            print(f"Open tarfile {self.jesInputArchivePath}{jecVersion}.{fileExt}")
        else:
            print(f"Open tarfile {self.jesInputArchivePath}{archive}.{fileExt}")
        self.jesArchive = tarfile.open(f"{self.jesInputArchivePath}{jecVersion}.{fileExt}", "r:gz") if not archive else tarfile.open(f"{self.jesInputArchivePath}{archive}.{fileExt}", "r:gz")
        self.jesInputFilePath = tempfile.mkdtemp()
        self.jesArchive.extractall(self.jesInputFilePath)

        if len(jesUncertainties) == 1 and jesUncertainties[0] == "Total":
            self.jesUncertaintyInputFileName =  f"{jecVersion}_Uncertainty_{jetType}.txt"
        elif jesUncertainties[0] == "Merged" and not self.isData:
            self.jesUncertaintyInputFileName = f"Regrouped_{jecVersion}_UncertaintySources_{jetType}.txt"
        else:
            self.jesUncertaintyInputFileName =  f"{jecVersion}_UncertaintySources_{jetType}.txt"

        # read all uncertainty source names from the loaded file
        if jesUncertainties[0] in ["All", "Merged"]:
            with open(f'{self.jesInputFilePath}/{self.jesUncertaintyInputFileName}') as f:
                lines = f.read().split("\n")
                sources = [x for x in lines if x.startswith("[") and x.endswith("]")]
                sources = [x[1:-1] for x in sources]
                self.jesUncertainties = sources
        if applyHEMfix:
            self.jesUncertainties.append("HEMIssue")

        self.jetReCalibrator = JetReCalibrator(jecVersion, jetType, True, self.jesInputFilePath,
            calculateSeparateCorrections=False, calculateType1METCorrection=False
        )

        self.jetSmearer = None
        if self.applySmearing:
            # smear jet pT to account for measured difference in JER between data
            # and simulation.
            if jerVersion != "":
                self.jerInputFileName = f"{jerVersion}_PtResolution_{jetType}.txt"
                self.jerUncertaintyInputFileName = f"{jerVersion}_SF_{jetType}.txt"
            # jet mass resolution: https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging
            self.jmrVals = jmrVals
            self.jetSmearer = JetSmearer(jecVersion, jetType, self.jerInputFileName, self.jerUncertaintyInputFileName, self.jmrVals)

        #############################################################################
        #
        # AK4 Puppi JEC & JER. To be applied on the soft-drop subjets.
        #
        #############################################################################
        if len(jesUncertainties) == 1 and jesUncertainties[0] == "Total":
            self.jesUncertaintyInputAK4PuppiFileName =  f"{jecVersion}_Uncertainty_AK4PFPuppi.txt"
        elif jesUncertainties[0] == "Merged" and not self.isData:
            self.jesUncertaintyInputAK4PuppiFileName = f"Regrouped_{jecVersion}_UncertaintySources_AK4PFPuppi.txt"
        else:
            self.jesUncertaintyInputAK4PuppiFileName =  f"{jecVersion}_UncertaintySources_AK4PFPuppi.txt"

        # read all uncertainty source names from the loaded file
        if jesUncertainties[0] in ["All", "Merged"]:
            with open(f'{self.jesInputFilePath}/{self.jesUncertaintyInputAK4PuppiFileName}') as f:
                lines = f.read().split("\n")
                sources = [x for x in lines if x.startswith("[") and x.endswith("]")]
                sources = [x[1:-1] for x in sources]
                self.jesUncertainties = sources

        self.jetReCalibratorAK4Puppi = JetReCalibrator(jecVersion, "AK4PFPuppi", True, self.jesInputFilePath,
            calculateSeparateCorrections=False, calculateType1METCorrection=False
        )

        #############################################################################
        #
        # load libraries for accessing JES scale factors and uncertainties
        # from txt files
        #
        #############################################################################
        for library in ["libCondFormatsJetMETObjects", "libPhysicsToolsNanoAODTools"]:
            if library not in ROOT.gSystem.GetLibraries():
                print("Load Library '{}'".format(library.replace("lib", "")))
                ROOT.gSystem.Load(library)

    def beginJob(self):

        print(f"Loading jet energy scale (JES) uncertainties from file '{os.path.join(self.jesInputFilePath,self.jesUncertaintyInputFileName)}'")
        # self.jesUncertainty = ROOT.JetCorrectionUncertainty(os.path.join(self.jesInputFilePath, self.jesUncertaintyInputFileName))

        #
        # AK8 Puppi jets JEC uncertainties
        #
        self.jesUncertainty = {}
        # implementation didn't seem to work for factorized JEC,try again
        # another way
        for jesUncertainty in self.jesUncertainties:
            jesUncertainty_label = jesUncertainty
            if jesUncertainty == 'Total' and (len(self.jesUncertainties) == 1 or len(self.jesUncertainties) == 2 and 'HEMIssue' in self.jesUncertainties):
                jesUncertainty_label = ''
            if jesUncertainty != "HEMIssue":
                pars = ROOT.JetCorrectorParameters(os.path.join(self.jesInputFilePath, self.jesUncertaintyInputFileName), jesUncertainty_label)
                self.jesUncertainty[jesUncertainty] = ROOT.JetCorrectionUncertainty(pars)

        #
        # AK4 Puppi jets JEC uncertainties for subjets
        #
        self.jesUncertaintyAK4Puppi = {}
        for jesUncertainty in self.jesUncertainties:
            jesUncertainty_label = jesUncertainty
            if jesUncertainty == 'Total' and (len(self.jesUncertainties) == 1 or len(self.jesUncertainties) == 2 and 'HEMIssue' in self.jesUncertainties):
                jesUncertainty_label = ''
            if jesUncertainty != "HEMIssue":
                pars = ROOT.JetCorrectorParameters(os.path.join(self.jesInputFilePath, self.jesUncertaintyInputAK4PuppiFileName), jesUncertainty_label)
                self.jesUncertaintyAK4Puppi[jesUncertainty] = ROOT.JetCorrectionUncertainty(pars)

        if self.applySmearing:
            self.jetSmearer.beginJob()

    def endJob(self):
        if self.applySmearing:
            self.jetSmearer.endJob()
        shutil.rmtree(self.jesInputFilePath)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch(f"{self.subJetBranchName}_pt_raw","F", lenVar=self.lenVarSubjets)
        self.out.branch(f"{self.subJetBranchName}_pt_AK4Puppicorr","F", lenVar=self.lenVarSubjets)
        self.out.branch(f"{self.subJetBranchName}_mass_raw","F", lenVar=self.lenVarSubjets)
        self.out.branch(f"{self.subJetBranchName}_mass_AK4Puppicorr","F", lenVar=self.lenVarSubjets)

        self.out.branch(f"{self.jetBranchName}_pt_raw","F", lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_pt_nom","F", lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_mass_raw","F", lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_mass_nom","F", lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_corr_JEC","F", lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_corr_JER","F", lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_msoftdrop_raw", "F", lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_msoftdrop_AK8corr", "F",lenVar=self.lenVar)
        self.out.branch(f"{self.jetBranchName}_msoftdrop_AK4Puppicorr", "F",lenVar=self.lenVar)
        if not self.isData:
            for shift in ["Up", "Down"]:
                for jesUncertainty in self.jesUncertainties:
                    self.out.branch(f"{self.jetBranchName}_pt_jes{jesUncertainty}{shift}", "F", lenVar=self.lenVar)
                    self.out.branch(f"{self.jetBranchName}_pt_jes{jesUncertainty}{shift}", "F", lenVar=self.lenVar)
                    self.out.branch(f"{self.jetBranchName}_mass_jes{jesUncertainty}{shift}", "F", lenVar=self.lenVar)
                    self.out.branch(f"{self.jetBranchName}_msoftdrop_AK8corr_jes{jesUncertainty}{shift}","F", lenVar=self.lenVar)
                self.out.branch(f"{self.jetBranchName}_pt_jer{shift}", "F", lenVar=self.lenVar)
                self.out.branch(f"{self.jetBranchName}_mass_jer{shift}", "F", lenVar=self.lenVar)
                self.out.branch(f"{self.jetBranchName}_msoftdrop_AK8corr_jer{shift}","F", lenVar=self.lenVar)
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        jets = Collection(event, self.jetBranchName)
        if not self.isData:
            genJets = Collection(event, self.genJetBranchName)

        subJets = Collection(event, self.subJetBranchName)
        if not self.isData:
            genSubJets = Collection(event, self.genSubJetBranchName)
            genSubJetMatcher = matchObjectCollectionMultiple(genJets, genSubJets, dRmax=0.8)

        if self.applySmearing:
            self.jetSmearer.setSeed(event)

        subjets_pt_raw = []
        subjets_pt_AK4Puppicorr = []
        subjets_mass_raw = []
        subjets_mass_AK4Puppicorr = []

        jets_pt_raw = []
        jets_pt_nom = []
        jets_mass_raw = []
        jets_mass_nom = []

        jets_corr_JEC = []
        jets_corr_JER = []

        jets_pt_jesUp = {}
        jets_pt_jesDown = {}
        jets_mass_jesUp = {}
        jets_mass_jesDown = {}
        for jesUncertainty in self.jesUncertainties:
            jets_pt_jesUp[jesUncertainty] = []
            jets_pt_jesDown[jesUncertainty] = []
            jets_mass_jesUp[jesUncertainty] = []
            jets_mass_jesDown[jesUncertainty] = []
        jets_pt_jerUp = []
        jets_pt_jerDown = []
        jets_mass_jerUp = []
        jets_mass_jerDown = []

        jets_msoftdrop_raw = []
        jets_msoftdrop_AK8corr = []
        jets_msoftdrop_AK8corr_jesUp = {}
        jets_msoftdrop_AK8corr_jesDown = {}
        for jesUncertainty in self.jesUncertainties:
            jets_msoftdrop_AK8corr_jesUp[jesUncertainty] = []
            jets_msoftdrop_AK8corr_jesDown[jesUncertainty] = []
        jets_msoftdrop_AK8corr_jerUp = []
        jets_msoftdrop_AK8corr_jerDown = []

        jets_msoftdrop_AK4Puppicorr = []

        rho = getattr(event, self.rhoBranchName)

        # match reconstructed jets to generator level ones
        # (needed to evaluate JER scale factors and uncertainties)
        if not self.isData:
            pairs = matchObjectCollection(jets, genJets)

        #
        # Loop over subjet collection. Get raw and re-apply AK4 Puppi correction
        #
        for subjet in subJets:
            subjet_pt = subjet.pt
            subjet_mass = subjet.mass
            if hasattr(subjet, "rawFactor"):
                subjet.p4_raw = subjet.p4() * (1 - subjet.rawFactor)
                subjet_pt_raw = subjet.p4_raw.Pt()
                subjet_mass_raw = subjet.p4_raw.M()
            else:
                subjet_pt_raw = -1.0 * subjet_pt  # If factor not present factor will be saved as -1
                subjet_mass_raw = -1.0 * subjet_mass  # If factor not present factor will be saved as -1
                subjet.p4_raw = ROOT.TLorentzVector()

            subjet.area = 1 #Dummy value. AK4PFPuppi corrections dont have L1 corrections so no need area (yet)
            subjet_pt_AK4Puppicorr, subjet_mass_AK4Puppicorr = self.jetReCalibratorAK4Puppi.correct(subjet, rho)
            subjet.p4_AK4Puppicorr = ROOT.TLorentzVector()
            subjet.p4_AK4Puppicorr.SetPtEtaPhiM(subjet_pt_AK4Puppicorr,subjet.eta,subjet.phi,subjet_mass_AK4Puppicorr)
            subjets_pt_raw.append(subjet_pt_raw)
            subjets_pt_AK4Puppicorr.append(subjet_pt_AK4Puppicorr)
            subjets_mass_raw.append(subjet_mass_raw)
            subjets_mass_AK4Puppicorr.append(subjet_mass_AK4Puppicorr)

        self.out.fillBranch(f"{self.subJetBranchName}_pt_raw",   subjets_pt_raw)
        self.out.fillBranch(f"{self.subJetBranchName}_pt_AK4Puppicorr",   subjets_pt_AK4Puppicorr)
        self.out.fillBranch(f"{self.subJetBranchName}_mass_raw", subjets_mass_raw)
        self.out.fillBranch(f"{self.subJetBranchName}_mass_AK4Puppicorr", subjets_mass_AK4Puppicorr)

        #
        # Loop over fatjet collection.
        #
        for jet in jets:
            # jet pt and mass corrections
            jet_pt = jet.pt
            jet_mass = jet.mass

            if hasattr(jet, "rawFactor"):
                jet_pt_raw = jet_pt * (1 - jet.rawFactor)
                jet_mass_raw = jet_mass * (1 - jet.rawFactor)
            else:
                jet_pt_raw = -1.0 * jet_pt  # If factor not present factor will be saved as -1
                jet_mass_raw = -1.0 * jet_mass  # If factor not present factor will be saved as -1

            (jet_pt, jet_mass) = self.jetReCalibrator.correct(jet, rho)
            jet.pt = jet_pt
            jet.mass = jet_mass
            jets_pt_raw.append(jet_pt_raw)
            jets_mass_raw.append(jet_mass_raw)

            jecAK8 = jet_pt / jet_pt_raw
            jets_corr_JEC.append(jecAK8)

            if not self.isData:
                genJet = pairs[jet]

            # evaluate JER scale factors and uncertainties
            # (cf. https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution and https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution )
            if self.applySmearing:
                (jet_pt_jerNomVal, jet_pt_jerUpVal, jet_pt_jerDownVal) = self.jetSmearer.getSmearValsPt(jet, genJet, rho)
            else:
                # set values to 1 for data so that jet_pt_nom is not smeared
                (jet_pt_jerNomVal, jet_pt_jerUpVal, jet_pt_jerDownVal) = (1., 1., 1.)
            jets_corr_JER.append(jet_pt_jerNomVal)

            jet_pt_nom = jet_pt_jerNomVal * jet_pt if self.applySmearing else jet_pt
            if jet_pt_nom < 0.0:
                jet_pt_nom *= -1.0
            jets_pt_nom.append(jet_pt_nom)

            jet_mass_nom = jet_pt_jerNomVal * jet_mass if self.applySmearing else jet_mass
            if jet_mass_nom < 0.0:
                jet_mass_nom *= -1.0
            jets_mass_nom.append(jet_mass_nom)

            if not self.isData:
                genGroomedSubJets = genSubJetMatcher[genJet] if genJet is not None else None
                genGroomedJet = genGroomedSubJets[0].p4() + genGroomedSubJets[1].p4() if genGroomedSubJets is not None and len(genGroomedSubJets) >= 2 else None
            else:
                genGroomedSubJets = None
                genGroomedJet = None

            #
            # Get raw subjets
            #
            subjet0 = None
            subjet1 = None
            if jet.subJetIdx1 >= 0 and jet.subJetIdx2 >= 0:
                subjet0 = subJets[jet.subJetIdx1]
                subjet1 = subJets[jet.subJetIdx2]
                groomedP4_raw = subjet0.p4_raw + subjet1.p4_raw
            else:
                groomedP4_raw = None
            jet_msoftdrop_raw = groomedP4_raw.M() if groomedP4_raw is not None else -1.
            jets_msoftdrop_raw.append(jet_msoftdrop_raw)

            #
            # Get msoftdrop after applying AK8 corrections on the subjets.
            #
            jet_msoftdrop_AK8corr = -1.
            if subjet0 and subjet1:
                jet_msoftdrop_AK8corr = jet_pt_jerNomVal * jecAK8 * jet_msoftdrop_raw  if self.applySmearing else jecAK8 * jet_msoftdrop_raw
            # store the soft-drop mass corrected with AK8 JEC (& JER)
            jets_msoftdrop_AK8corr.append(jet_msoftdrop_AK8corr)

            #
            # Get msoftdrop usng the AK4 Puppi JEC-corrected subjets.
            #
            jet_msoftdrop_AK4Puppicor = -1.
            if subjet0 and subjet1:
                jet_msoftdrop_AK4Puppicor = (subjet0.p4_AK4Puppicorr + subjet1.p4_AK4Puppicorr).M()
            jets_msoftdrop_AK4Puppicorr.append(jet_msoftdrop_AK4Puppicor)

            if not self.isData:
                #
                # JES
                #
                jet_pt_jesUp = {}
                jet_pt_jesDown = {}
                jet_mass_jesUp = {}
                jet_mass_jesDown = {}
                jet_msoftdrop_AK8corr_jesUp = {}
                jet_msoftdrop_AK8corr_jesDown = {}
                for jesUncertainty in self.jesUncertainties:
                    # (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorUncertainties)
                    # cf. https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
                    if jesUncertainty == "HEMIssue":
                        delta = 1.
                        if jet_pt_nom > 15 and jet.jetId & 2 and jet.phi > -1.57 and jet.phi < -0.87:
                            if jet.eta > -2.5 and jet.eta < -1.3:
                                delta = 0.8
                            elif jet.eta <= -2.5 and jet.eta > -3:
                                delta = 0.65
                        jet_pt_jesUp[jesUncertainty] = jet_pt_nom
                        jet_pt_jesDown[jesUncertainty] = delta * jet_pt_nom
                        jet_mass_jesUp[jesUncertainty] = jet_mass_nom
                        jet_mass_jesDown[jesUncertainty] = delta * jet_mass_nom
                        jet_msoftdrop_AK8corr_jesUp[jesUncertainty] = -1.
                        jet_msoftdrop_AK8corr_jesDown[jesUncertainty] = -1.
                        if subjet0 and subjet1:
                            jet_msoftdrop_AK8corr_jesUp[jesUncertainty] = jet_msoftdrop_AK8corr
                            jet_msoftdrop_AK8corr_jesDown[jesUncertainty] = delta * jet_msoftdrop_AK8corr
                    else:
                        self.jesUncertainty[jesUncertainty].setJetPt(jet_pt_nom)
                        self.jesUncertainty[jesUncertainty].setJetEta(jet.eta)
                        delta = self.jesUncertainty[jesUncertainty].getUncertainty(True)
                        jet_pt_jesUp[jesUncertainty] = jet_pt_nom * (1. + delta)
                        jet_pt_jesDown[jesUncertainty] = jet_pt_nom * (1. - delta)
                        jet_mass_jesUp[jesUncertainty] = jet_mass_nom * (1. + delta)
                        jet_mass_jesDown[jesUncertainty] = jet_mass_nom * (1. - delta)
                        jet_msoftdrop_AK8corr_jesUp[jesUncertainty] = -1.
                        jet_msoftdrop_AK8corr_jesDown[jesUncertainty] = -1.
                        if subjet0 and subjet1:
                            jet_msoftdrop_AK8corr_jesUp[jesUncertainty] = jet_msoftdrop_AK8corr * (1. + delta)
                            jet_msoftdrop_AK8corr_jesDown[jesUncertainty] = jet_msoftdrop_AK8corr * (1. - delta)

                    jets_pt_jesUp[jesUncertainty].append(jet_pt_jesUp[jesUncertainty])
                    jets_pt_jesDown[jesUncertainty].append(jet_pt_jesDown[jesUncertainty])
                    jets_mass_jesUp[jesUncertainty].append(jet_mass_jesUp[jesUncertainty])
                    jets_mass_jesDown[jesUncertainty].append(jet_mass_jesDown[jesUncertainty])
                    jets_msoftdrop_AK8corr_jesUp[jesUncertainty].append(jet_msoftdrop_AK8corr_jesUp[jesUncertainty])
                    jets_msoftdrop_AK8corr_jesDown[jesUncertainty].append(jet_msoftdrop_AK8corr_jesDown[jesUncertainty])
                #
                # JER
                #
                jet_pt_jerUp = jet_pt_jerUpVal * jet_pt
                jet_pt_jerDown = jet_pt_jerDownVal * jet_pt
                jet_mass_jerUp = jet_pt_jerUpVal *  jet_mass
                jet_mass_jerDown = jet_pt_jerDownVal * jet_mass

                jets_pt_jerUp.append(jet_pt_jerUp)
                jets_pt_jerDown.append(jet_pt_jerDown)
                jets_mass_jerUp.append(jet_mass_jerUp)
                jets_mass_jerDown.append(jet_mass_jerDown)

                jet_msoftdrop_AK8corr_jerUp = -1.
                jet_msoftdrop_AK8corr_jerDown = 1.
                if subjet0 and subjet1:
                    jet_msoftdrop_AK8corr_jerUp   = jet_pt_jerUpVal   * jecAK8 * jet_msoftdrop_raw  if self.applySmearing else jecAK8 * jet_msoftdrop_raw
                    jet_msoftdrop_AK8corr_jerDown = jet_pt_jerDownVal * jecAK8 * jet_msoftdrop_raw  if self.applySmearing else jecAK8 * jet_msoftdrop_raw
                jets_msoftdrop_AK8corr_jerUp.append(jet_msoftdrop_AK8corr_jerUp)
                jets_msoftdrop_AK8corr_jerDown.append(jet_msoftdrop_AK8corr_jerDown)

        self.out.fillBranch(f"{self.jetBranchName}_pt_raw",   jets_pt_raw)
        self.out.fillBranch(f"{self.jetBranchName}_pt_nom",   jets_pt_nom)
        self.out.fillBranch(f"{self.jetBranchName}_mass_raw", jets_mass_raw)
        self.out.fillBranch(f"{self.jetBranchName}_mass_nom", jets_mass_nom)
        self.out.fillBranch(f"{self.jetBranchName}_corr_JEC", jets_corr_JEC)
        if not self.isData:
            self.out.fillBranch(f"{self.jetBranchName}_corr_JER",     jets_corr_JER)
            for jesUncertainty in self.jesUncertainties:
                self.out.fillBranch(f"{self.jetBranchName}_pt_jes{jesUncertainty}Up",     jets_pt_jesUp[jesUncertainty])
                self.out.fillBranch(f"{self.jetBranchName}_pt_jes{jesUncertainty}Down",   jets_pt_jesDown[jesUncertainty])
                self.out.fillBranch(f"{self.jetBranchName}_mass_jes{jesUncertainty}Up",   jets_mass_jesUp[jesUncertainty])
                self.out.fillBranch(f"{self.jetBranchName}_mass_jes{jesUncertainty}Down", jets_mass_jesDown[jesUncertainty])
            self.out.fillBranch(f"{self.jetBranchName}_pt_jerUp",     jets_pt_jerUp)
            self.out.fillBranch(f"{self.jetBranchName}_pt_jerDown",   jets_pt_jerDown)
            self.out.fillBranch(f"{self.jetBranchName}_mass_jerUp",   jets_mass_jerUp)
            self.out.fillBranch(f"{self.jetBranchName}_mass_jerDown", jets_mass_jerDown)
        #
        #
        #
        self.out.fillBranch(f"{self.jetBranchName}_msoftdrop_raw", jets_msoftdrop_raw)
        self.out.fillBranch(f"{self.jetBranchName}_msoftdrop_AK8corr", jets_msoftdrop_AK8corr)
        if not self.isData:
            for jesUncertainty in self.jesUncertainties:
                self.out.fillBranch(f"{self.jetBranchName}_msoftdrop_AK8corr_jes{jesUncertainty}Up",   jets_msoftdrop_AK8corr_jesUp[jesUncertainty])
                self.out.fillBranch(f"{self.jetBranchName}_msoftdrop_AK8corr_jes{jesUncertainty}Down", jets_msoftdrop_AK8corr_jesDown[jesUncertainty])
            self.out.fillBranch(f"{self.jetBranchName}_msoftdrop_AK8corr_jerUp",   jets_msoftdrop_AK8corr_jerUp)
            self.out.fillBranch(f"{self.jetBranchName}_msoftdrop_AK8corr_jerDown", jets_msoftdrop_AK8corr_jerDown)
        #
        #
        #
        self.out.fillBranch(f"{self.jetBranchName}_msoftdrop_AK4Puppicorr", jets_msoftdrop_AK4Puppicorr)





        return True

