import FWCore.ParameterSet.Config as cms
import pprint
isMC = True

process = cms.Process("demo")

##---------  Load standard Reco modules ------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')



##----- this config frament brings you the generator information ----
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("Configuration.StandardSequences.Generator_cff")


##----- Detector geometry : some of these needed for b-tag -------
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")


##----- B-tags --------------
process.load("RecoBTag.Configuration.RecoBTag_cff")


##----- Global tag: conditions database ------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

##----- Counter module ------------
#process.load("ElectroWeakAnalysis.VPlusJets.AllPassFilter_cfi")

## import skeleton process
#from PhysicsTools.PatAlgos.patTemplate_cfg import *


############################################
if not isMC:
    process.GlobalTag.globaltag = 'GR_R_53_V10::All'
else:
    process.GlobalTag.globaltag = 'START53_V7E::All'

OutputFileName = "WmunuJetAnalysisntuple.root"
numEventsToRun = 100
############################################
########################################################################################
########################################################################################
#
###---------  W-->munu Collection ------------
#process.load("ElectroWeakAnalysis.VPlusJets.WmunuCollectionsPAT_cfi")
#
###---------  Jet Collection ----------------
#process.load("ElectroWeakAnalysis.VPlusJets.JetCollectionsPAT_cfi")
#
###---------  Vertex and track Collections -----------
#process.load("ElectroWeakAnalysis.VPlusJets.TrackCollections_cfi")
##


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(numEventsToRun)
)

process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

#process.options = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound')
#)
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
#       '/store/user/lnujj/PatTuples_8TeV_53X-v1/jdamgov/SingleMu/SQWaT_PAT_53X_2012B-13Jul2012-v1_part1v3/3e4086321697e2c39c90dad08848274b/pat_53x_test_v03_data_9_0_BNS.root'
#       '/store/user/lnujj/PatTuples_8TeV_53X-v1/jdamgov/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/SQWaT_PAT_53X_Summer12_v1/829f288d768dd564418efaaf3a8ab9aa/pat_53x_test_v03_995_1_wBa.root'
        '/store/user/lnujj/PatTuples_8TeV_53X-v1/zixu/GluGluToHToWWToLAndTauNuQQ_M-800_8TeV-powheg-pythia6/SQWaT_PAT_53X_ggH800/829f288d768dd564418efaaf3a8ab9aa/pat_53x_test_v03_ggH800_9_1_Zyr.root',                                                                            
) )




##-------- Muon events of interest --------
process.HLTMu =cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
     HLTPaths = cms.vstring('HLT_IsoMu24_*','HLT_IsoMu30_*'),
     eventSetupPathsKey = cms.string(''),
     andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
     throw = cms.bool(False) # throw exception on unknown path names
)

process.primaryVertex = cms.EDFilter("VertexSelector",
                             src = cms.InputTag("offlinePrimaryVertices"),                                
                             cut = cms.string("!isFake && ndof >= 4 && abs(z) <= 24 && position.Rho <= 2"), # tracksSize() > 3 for the older cut
                             filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
                             )
process.primaryVertex.src = cms.InputTag("goodOfflinePrimaryVertices");
process.primaryVertex.cut = cms.string(" ");

##-------- Save V+jets trees --------
process.substructureTester = cms.EDAnalyzer("SubstructureTester", 
#    jetType = cms.string("PF"),
#  #  srcPFCor = cms.InputTag("selectedPatJetsPFlow"),
#    srcPFCor = cms.InputTag("ak5PFJetsLooseId"),
#    srcPhoton = cms.InputTag("photons"),
#    IsoValPhoton = cms.VInputTag(cms.InputTag('phoPFIso:chIsoForGsfEle'),
#                                 cms.InputTag('phoPFIso:phIsoForGsfEle'),
#                                 cms.InputTag('phoPFIso:nhIsoForGsfEle'),
#                                                           ),
#    srcPFCorVBFTag = cms.InputTag("ak5PFJetsLooseIdVBFTag"), 
#    srcVectorBoson = cms.InputTag("bestWmunu"),
#    VBosonType     = cms.string('W'),
#    LeptonType     = cms.string('muon'),                          
#    TreeName    = cms.string('WJet'),
#    srcPrimaryVertex = cms.InputTag("goodOfflinePrimaryVertices"),                               
#    runningOverMC = cms.bool(isMC),			
#    runningOverAOD = cms.bool(False),			
##    srcMet = cms.InputTag("patType1CorrectedPFMet"),
#    srcMet = cms.InputTag("patMetShiftCorrected"), # type1 + shift corrections
#    srcMetMVA = cms.InputTag("pfMEtMVA"),
#    srcGen  = cms.InputTag("ak5GenJets"),
#    srcMuons  = cms.InputTag("selectedPatMuonsPFlow"),
#    srcBeamSpot  = cms.InputTag("offlineBeamSpot"),
#    srcCaloMet  = cms.InputTag("patMETs"),
#    srcgenMet  = cms.InputTag("genMetTrue"),
#    srcGenParticles  = cms.InputTag("genParticles"),
#    srcTcMet    = cms.InputTag("patMETsAK5TC"),
    srcJetsforRho = cms.string("kt6PFJetsPFlow"),                               
#    srcJetsforRho_lepIso = cms.string("kt6PFJetsForIsolation"),       
#    srcJetsforRhoCHS = cms.string("kt6PFJetsChsPFlow"),
#    srcJetsforRho_lepIsoCHS = cms.string("kt6PFJetsChsForIsolationPFlow"),
#    srcFlavorByValue = cms.InputTag("ak5tagJet"),
#    bTagger=cms.string("simpleSecondaryVertexHighEffBJetTags"),
#
#    applyJECToGroomedJets=cms.bool(True),
##    doGroomedAK5 = cms.bool(True),
#    doGroomedAK5 = cms.bool(False),                                   
##    doGroomedAK7 = cms.bool(True),
#    doGroomedAK7 = cms.bool(False),                                   
#    doGroomedAK8 = cms.bool(False),
#    doGroomedCA8 = cms.bool(True),
#    doGroomedCA12 = cms.bool(False),
#                                   
#    GroomedJet_doQJets = cms.bool( True ),
#    GroomedJet_QJetsPreclustering = cms.int32(50), 
#    GroomedJet_QJetsN = cms.int32(100)                                    
)

#if isMC:
#    process.VplusJets.JEC_GlobalTag_forGroomedJet = cms.string("START53_V15")
#else:
#    process.VplusJets.JEC_GlobalTag_forGroomedJet = cms.string("FT_53_V10_AN3")


process.TFileService = cms.Service(
    "TFileService", 
    fileName = cms.string( OutputFileName ),
    closeFileFast = cms.untracked.bool(False)
)

process.myseq = cms.Sequence(
    process.HLTMu
    )

if isMC:
    process.myseq.remove ( process.HLTMu)

#options.wantSummary = True
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

##---- if do not want to require >= 2 jets then disable that filter ---
##process.myseq.remove ( process.RequireTwoJets)  

#process.outpath.remove(process.out)
#process.p = cms.Path( process.myseq  * process.substructureTester)
process.p = cms.Path( process.substructureTester )
#try git
