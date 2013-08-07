#ifndef JETSUBSTRUCTURETOOLS_H
#define JETSUBSTRUCTURETOOLS_H


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h" 

//#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEventAuxInfo.h"
//#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"

#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"

//#include "VHbbAnalysis/VHbbDataFormats/interface/Nsubjettiness.h"
#include "JetSubstructure/SubstructureTools/interface/NjettinessPlugin.hh"
#include "JetSubstructure/SubstructureTools/interface/Nsubjettiness.hh"
#include "JetSubstructure/SubstructureTools/interface/GeneralizedEnergyCorrelator.hh"
#include "JetSubstructure/SubstructureTools/src/QjetsPlugin.h"

#include <iostream>

class JetSubstructureTools {
    
    // member functions
 public:
    
    // constructor
    JetSubstructureTools( fastjet::JetDefinition jetDef, std::vector< fastjet::PseudoJet > constits, float area );    
    ~JetSubstructureTools(){}
    
    // Grooming
    fastjet::PseudoJet getJet();
    fastjet::PseudoJet getJet_wGhostes();
    fastjet::PseudoJet getJet_basic();
    fastjet::PseudoJet getPrunedJet( float zcut = 0.1, float rcut = 0.5 );
    fastjet::PseudoJet getFilteredJet( float rfilt = 0.3, int nfilt = 3 );
    fastjet::PseudoJet getTrimmedJet( float rfilt = 0.2, float ptfrac = 0.03 );
    std::vector<fastjet::PseudoJet> getPrunedSubjets( int nToKeep = 2, float zcut = 0.1, float rcut = 0.5 );

    // N-subjettiness
    float getTau( int N, float kappa = 1 );

    // Qjets
    double getQjetVolatility(int seed, int QJetsN = 50, int QJetsPreclustering = 30 );

    // Jet charge
    double getJetCharge( float kappa );
    
    // Generalized Energy correlations
    double getJetECF( float beta );

    
  private: 
    float FindRMS( std::vector< float > qjetmasses );
    float FindMean( std::vector< float > qjetmasses );
    float getPdgIdCharge( float fid );
    
    double mJetRadius;
    int nJets_;
    
    fastjet::ClusterSequenceArea *thisClustering_;
    fastjet::ClusterSequenceArea *thisClustering_wGhosts_;    
    fastjet::ClusterSequence *thisClustering_basic_;
    std::vector<fastjet::PseudoJet> out_jets_;    
    std::vector<fastjet::PseudoJet> out_jets_wGhosts_;    
    std::vector<fastjet::PseudoJet> out_jets_basic_;        
    
    fastjet::PseudoJet jet0_full_;
    fastjet::PseudoJet jet0_full_wGhosts_;    
    fastjet::PseudoJet jet0_basic_;   
    
    std::vector<int> neutrals;
    std::vector<int> positives;
    std::vector<int> negatives;        
    
 
  public:
    
    std::vector<fastjet::PseudoJet> FJconstituents_;
    std::vector<float> c_pdgIds_; 
        
};

#endif






