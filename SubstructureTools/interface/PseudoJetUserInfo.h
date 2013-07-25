#ifndef PseudoJetUserInfo_h
#define PseudoJetUserInfo_h

// system include files
#include <memory>
#include <string>
#include <iostream>
#include <map>
#include <fstream>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"

class PseudoJetUserInfo : public fastjet::PseudoJet::UserInfoBase{
public:
    // default ctor
    //  - pdg_id        the PDG id of the particle
    //  - charge        the charge of the particle
    PseudoJetUserInfo(const int & pdg_id, const int & charge) :
    _pdg_id(pdg_id), _charge(charge){}
    
    /// access to the PDG id
    int pdg_id() const { return _pdg_id; }
    
    /// access to the vertex number
    int charge() const { return _charge; }
    
protected:
    int _pdg_id;         // the associated pdg id
    int _charge;  // the associated vertex number
};

#endif
