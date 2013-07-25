// system include files
#include <memory>
#include <string>
#include <iostream>
#include <map>
#include <fstream>
    
// user include files
#include "JetSubstructure/SubstructureTools/interface/JetSubstructureTools.h" 
#include "JetSubstructure/SubstructureTools/interface/PseudoJetUserInfo.h"

///////////////////////////////////
////// some helpers
float JetSubstructureTools::FindRMS( std::vector< float > qjetmasses ){
    
    float total = 0.;
    float ctr = 0.;
    for (unsigned int i = 0; i < qjetmasses.size(); i++){
        total = total + qjetmasses[i];
        ctr++;
    }
    float mean =  total/ctr;
    
    float totalsquared = 0.;
    for (unsigned int i = 0; i < qjetmasses.size(); i++){
        totalsquared += (qjetmasses[i] - mean)*(qjetmasses[i] - mean) ;
    }
    float RMS = sqrt( totalsquared/ctr );
    return RMS;
}

float JetSubstructureTools::FindMean( std::vector< float > qjetmasses ){
    float total = 0.;
    float ctr = 0.;
    for (unsigned int i = 0; i < qjetmasses.size(); i++){
        total = total + qjetmasses[i];
        ctr++;
    }
    return total/ctr;
}
///////////////////////////////////



// -------------------------------------------
// Constructor
// -------------------------------------------
JetSubstructureTools::JetSubstructureTools( fastjet::JetDefinition jetDef, std::vector< fastjet::PseudoJet > constits, float area )
{
    
    
//    mJetRadius = radii;
//    
//    // check that they are all the same size
//    if ((c_px.size() == c_py.size())&&(c_py.size() == c_pz.size())&&(c_pz.size() == c_e.size())&&(c_e.size() == c_pdgId.size())){
//        for (unsigned int i = 0; i < c_px.size(); i++){
//            
//            FJconstituents_.push_back( fastjet::PseudoJet( c_px[i], c_py[i], c_pz[i], c_e[i] ) );
//            c_pdgIds_.push_back( c_pdgId[i] );
//            
//        }
//    }
//    else throw cms::Exception("JetSubstructureTools") << "Constituent size mismatch..." << std::endl;
//    
//    // -------------------------------------------
//    // recluster on the fly....
//    fastjet::JetDefinition jetDef(fastjet::cambridge_algorithm, mJetRadius);
    
//    FJconstituents_ = constits;    
//    
//    for (unsigned int i = 0; i < FJconstituents_.size(); i++){
//        
//        std::cout << "has info = " << FJconstituents_[i].has_user_info() << std::endl;
//        if ( FJconstituents_[i].has_user_info() ){
//         
//            std::cout << "charge = " << FJconstituents_[i].user_info<PseudoJetUserInfo>().charge() << std::endl;
//            
//        }
//        
//    }    
        
    FJconstituents_ = constits;
    
    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
//    fjActiveArea.set_fj2_placement(true);
    fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
    fastjet::AreaDefinition fjAreaDefinition_wGhosts( fastjet::active_area_explicit_ghosts, fjActiveArea );
    
    thisClustering_ = new fastjet::ClusterSequenceArea(FJconstituents_, jetDef, fjAreaDefinition);
    thisClustering_wGhosts_ = new fastjet::ClusterSequenceArea(FJconstituents_, jetDef, fjAreaDefinition_wGhosts);
    thisClustering_basic_ = new fastjet::ClusterSequence(FJconstituents_, jetDef);

    out_jets_ = sorted_by_pt(thisClustering_->inclusive_jets(5.0));
    out_jets_wGhosts_ = sorted_by_pt(thisClustering_wGhosts_->inclusive_jets(5.0));    
    out_jets_basic_ = sorted_by_pt(thisClustering_basic_->inclusive_jets(5.0));
    
    jet0_full_ = out_jets_[0];
    jet0_full_wGhosts_ = out_jets_wGhosts_[0];
    jet0_basic_ = out_jets_basic_[0];
        
    // define charges of pdgIds for jet charge
    neutrals.push_back( 22 ); neutrals.push_back( 130 ); neutrals.push_back( 310 ); neutrals.push_back( 311 ); neutrals.push_back( 111 ); 
    neutrals.push_back( 1 ); neutrals.push_back( 2 ); neutrals.push_back( 3 ); neutrals.push_back( 4 ); neutrals.push_back( 5 ); 
    neutrals.push_back( -1 ); neutrals.push_back( -2 ); neutrals.push_back( -3 ); neutrals.push_back( -4 ); neutrals.push_back( -5 ); 
    neutrals.push_back( 2112 );
    
    positives.push_back( 321 ); positives.push_back( 211 ); ; positives.push_back( -11 ); positives.push_back( -13); positives.push_back( 2212);
    negatives.push_back( -321 ); negatives.push_back( -211 ); negatives.push_back( 11 ); negatives.push_back( 13 );    
    
}

// -------------------------------------------
// N-subjettiness
// -------------------------------------------
float JetSubstructureTools::getTau( int N, float kappa ){
    
    fastjet::Nsubjettiness nSubNKT(N, Njettiness::onepass_kt_axes, kappa, mJetRadius, mJetRadius);
    float tauN = nSubNKT(jet0_full_);
    return tauN;
    
}

// -------------------------------------------
// Groomers
// -------------------------------------------
fastjet::PseudoJet JetSubstructureTools::getJet(){
    
    return jet0_full_;
    
}
fastjet::PseudoJet JetSubstructureTools::getPrunedJet( float zcut, float rcut ){
    
    fastjet::Pruner pruner(fastjet::cambridge_algorithm, zcut, rcut);
    fastjet::PseudoJet prunedJet = pruner( jet0_full_wGhosts_ );
    return prunedJet;
    
}
fastjet::PseudoJet JetSubstructureTools::getFilteredJet( float rfilt, int nfilt ){
    
    fastjet::Filter filter( fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, rfilt), fastjet::SelectorNHardest(nfilt)));
    fastjet::PseudoJet filteredJet = filter( jet0_full_wGhosts_ );
    return filteredJet;
    
}
fastjet::PseudoJet JetSubstructureTools::getTrimmedJet( float rfilt, float ptfrac ){
    
    fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03)));
    fastjet::PseudoJet trimmedJet = trimmer( jet0_full_wGhosts_ );
    return trimmedJet;
    
}
    
std::vector< fastjet::PseudoJet > JetSubstructureTools::getPrunedSubjets( int nToKeep, float zcut, float rcut ){
    
    fastjet::Pruner pruner(fastjet::cambridge_algorithm, zcut, rcut);
    fastjet::PseudoJet prunedJet = pruner( jet0_basic_ );
    std::vector<fastjet::PseudoJet> subjets = jet0_basic_.associated_cluster_sequence()->exclusive_subjets(jet0_basic_,nToKeep);    
    return subjets;
    
}

// -------------------------------------------
// QJets
// -------------------------------------------
double JetSubstructureTools::getQjetVolatility(int seed, int QJetsN, int QJetsPreclustering){
    
    
    std::vector< float > qjetmasses;
    
    double zcut(0.1), dcut_fctr(0.5), exp_min(0.), exp_max(0.), rigidity(0.1);          
    
    vector<fastjet::PseudoJet> constits;
    unsigned int nqjetconstits = jet0_full_.constituents().size();
    if (nqjetconstits < (unsigned int) QJetsPreclustering) constits = jet0_full_.constituents();
    else constits = jet0_full_.associated_cluster_sequence()->exclusive_subjets_up_to(jet0_full_,QJetsPreclustering);
    for(unsigned int ii = 0 ; ii < (unsigned int) QJetsN ; ii++){
        QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity);
        qjet_plugin.SetRandSeed(seed+ii); // new feature in Qjets to set the random seed
        fastjet::JetDefinition qjet_def(&qjet_plugin);
        fastjet::ClusterSequence qjet_seq(constits, qjet_def);
        vector<fastjet::PseudoJet> inclusive_jets2 = sorted_by_pt(qjet_seq.inclusive_jets(5.0));
        
        if (inclusive_jets2.size()>0) { qjetmasses.push_back( inclusive_jets2[0].m() ); }
        
    }
    
        // find RMS of a vector
    float qjetsRMS = FindRMS( qjetmasses );
        // find mean of a vector
    float qjetsMean = FindMean( qjetmasses );
    float qjetsVolatility = qjetsRMS/qjetsMean;
    return qjetsVolatility;
}

// -------------------------------------------
// Jet Charge
// -------------------------------------------
float JetSubstructureTools::getPdgIdCharge( float fid ){
    
    float qq = -99.;
    int id = (int) fid;
    if (std::find(neutrals.begin(), neutrals.end(), id) != neutrals.end()){
        qq = 0.;
    }
    else if (std::find(positives.begin(), positives.end(), id) != positives.end()){
        qq = 1.;
    }
    else if (std::find(negatives.begin(), negatives.end(), id) != negatives.end()){
        qq = -1.;
    }
    else{
        throw cms::Exception("GroomedJetFiller") << " unknown PDG id " << id << std::endl;
    }
    return qq;
}

double JetSubstructureTools::getJetCharge( float kappa ){
    
    float PTjet = jet0_basic_.pt(); 
    
    float val = 0.;
    for (unsigned int i = 0; i < FJconstituents_.size(); i++){
        float qq = FJconstituents_[i].user_info<PseudoJetUserInfo>().charge();
        val += qq*pow(FJconstituents_.at(i).pt(),kappa);
    }
    val /= PTjet;
    return val;
    
}

// -------------------------------------------
// Generalized Energy Correlations
// -------------------------------------------

double JetSubstructureTools::getJetECF( float beta ){

    // Generalized energy correlator
    fastjet::JetDefinition jet_def_forECF(fastjet::antikt_algorithm, 2.0);
    fastjet::ClusterSequence clust_seq_forECF(FJconstituents_, jet_def_forECF);
    vector<fastjet::PseudoJet> incluisve_jets_forECF = clust_seq_forECF.inclusive_jets(0);
    fastjet::GeneralizedEnergyCorrelatorRatio C2beta(2,beta,fastjet::pT_R); // beta = 1.7
    double jetGeneralizedECF = C2beta(incluisve_jets_forECF[0]);
    return jetGeneralizedECF;
    
}

