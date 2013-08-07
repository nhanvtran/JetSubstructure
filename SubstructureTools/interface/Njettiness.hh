//  Njettiness Package
//  Version 0.4.1 (January 22, 2012)
//  Questions/Comments?  jthaler@jthaler.net

// Copyright (c) 2011-12, Jesse Thaler, Ken Van Tilburg, and Christopher K.
// Vermilion
//
//----------------------------------------------------------------------
// This file is part of the N-jettiness package ("N-jettiness").
//
//  N-jettiness is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  SpartyJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with SpartyJet; if not, write to the Free Software
//  Foundation, Inc.:
//      59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//----------------------------------------------------------------------


#ifndef __NJETTINESS_HH__
#define __NJETTINESS_HH__


#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <cmath>
#include <vector>
#include <list>

///////
//
// Helper classes
//
///////


// Parameters that define Njettiness
class NsubParameters {
private:
   double _beta;  // angular weighting exponent
   double _R0;    // characteristic jet radius (for normalization)
   double _Rcutoff;  // Cutoff scale for cone jet finding (default is large number such that no boundaries are used)
   
public:
   NsubParameters(const double mybeta, const double myR0, const double myRcutoff=10000.0) :
   _beta(mybeta), _R0(myR0), _Rcutoff(myRcutoff) {}
   double beta() const {return _beta;}
   double R0() const {return _R0;}
   double Rcutoff() const {return _Rcutoff;}
};

// Parameters that change minimization procedure.
// Set automatically when you choose NsubAxesMode, but can be adjusted manually as well
class KmeansParameters {
private:
   int _n_iterations;  // Number of iterations to run  (0 for no minimization, 1 for one-pass, >>1 for global minimum)
   double _precision;  // Desired precision in axes alignment
   int _halt;          // maximum number of steps per iteration
   double _noise_range;// noise range for random initialization
   
public:
   KmeansParameters() : _n_iterations(0), _precision(0.0), _halt(0), _noise_range(0.0) {}
   KmeansParameters(const int my_n_iterations, double my_precision, int my_halt, double my_noise_range) :
   _n_iterations(my_n_iterations),  _precision(my_precision), _halt(my_halt), _noise_range(my_noise_range) {}
   int n_iterations() const { return _n_iterations;}
   double precision() const {return _precision;}
   int halt() const {return _halt;}
   double noise_range() const {return _noise_range;}
};

// helper class for minimization
class LightLikeAxis {
private:
   double _rap, _phi, _weight, _mom;
   
public:
   LightLikeAxis() : _rap(0.0), _phi(0.0), _weight(0.0), _mom(0.0) {}
   LightLikeAxis(double my_rap, double my_phi, double my_weight, double my_mom) :
   _rap(my_rap), _phi(my_phi), _weight(my_weight), _mom(my_mom) {}
   double rap() const {return _rap;}
   double phi() const {return _phi;}
   double weight() const {return _weight;}
   double mom() const {return _mom;}
   void set_rap(double my_set_rap) {_rap = my_set_rap;}
   void set_phi(double my_set_phi) {_phi = my_set_phi;}
   void set_weight(double my_set_weight) {_weight = my_set_weight;}
   void set_mom(double my_set_mom) {_mom = my_set_mom;}
   void reset(double my_rap, double my_phi, double my_weight, double my_mom) {_rap=my_rap; _phi=my_phi; _weight=my_weight; _mom=my_mom;}
};

///////
//
// Functions for minimization.
// TODO:  Wrap these in N-subjettiness class
//
///////

inline double sq(double x);// {return x*x;}

// Calculates distance between two points in rapidity-azimuth plane
double DistanceSq(double rap1, double phi1, double rap2, double phi2);/* {
   double distRap = rap1-rap2;
   double distPhi = std::fabs(phi1-phi2);
   if (distPhi > M_PI) {distPhi = 2.0*M_PI - distPhi;}
   return sq(distRap) + sq(distPhi);
}*/

inline double Distance(double rap1, double phi1, double rap2, double phi2) ;

// Given starting axes, update to find better axes
template <int N>
std::vector<LightLikeAxis> UpdateAxesFast(const std::vector <LightLikeAxis> & old_axes, 
                                  const std::vector <fastjet::PseudoJet> & inputJets,
                                  NsubParameters paraNsub, double precision);

// Given starting axes, update to find better axes
// (This is just a wrapper for the templated version above.)
std::vector<LightLikeAxis> UpdateAxes(const std::vector <LightLikeAxis> & old_axes, 
                                      const std::vector <fastjet::PseudoJet> & inputJets, NsubParameters paraNsub, double precision) ;


// Go from internal LightLikeAxis to PseudoJet
// TODO:  Make part of LightLikeAxis class.
std::vector<fastjet::PseudoJet> ConvertToPseudoJet(const std::vector <LightLikeAxis>& axes) ;

// N-subjettiness pieces
std::vector<double> ConstituentTauValue(const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& axes, const NsubParameters& paraNsub);/* {// Returns the sub-tau values, i.e. a std::vector of the contributions to tau_N of each Voronoi region (or region within R_0)
   double beta = paraNsub.beta();
   double R0 = paraNsub.R0();
   double Rcutoff = paraNsub.Rcutoff();
   
   std::vector<double> tauNum(axes.size(), 0.0), tau(axes.size());
   double tauDen = 0.0;
   int j_min = 0;
   for (unsigned i = 0; i < particles.size(); i++) {
      // find minimum distance; start with 0'th axis for reference
      double minR = std::sqrt(particles[i].squared_distance(axes[0]));
      for (unsigned j = 1; j < axes.size(); j++) {
         double tempR = std::sqrt(particles[i].squared_distance(axes[j])); // delta R distance
         if (tempR < minR) {minR = tempR; j_min = j;}
      }
      if (minR > Rcutoff) {minR = Rcutoff;}
      tauNum[j_min] += particles[i].perp() * std::pow(minR,beta);
      tauDen += particles[i].perp() * std::pow(R0,beta);
   }
   for (unsigned j = 0; j < axes.size(); j++) {
      tau[j] = tauNum[j]/tauDen;
   }
   return tau;
} */

// N-subjettiness values
double TauValue(const std::vector <fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes,const NsubParameters& paraNsub) ;

// Get exclusive kT subjets
std::vector<fastjet::PseudoJet> GetKTAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputJets) ;

// Get exclusive CA subjets
std::vector<fastjet::PseudoJet> GetCAAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputJets) ;

// Get inclusive anti kT hardest subjets
std::vector<fastjet::PseudoJet> GetAntiKTAxes(int n_jets, double R0, const std::vector <fastjet::PseudoJet> & inputJets) ;


// Get minimization axes
std::vector<fastjet::PseudoJet> GetMinimumAxes(const std::vector <fastjet::PseudoJet> & seedAxes, const std::vector <fastjet::PseudoJet> & inputJets, KmeansParameters para, 
                                          NsubParameters paraNsub) ;

///////
//
// Main Njettiness Class
//
///////

class Njettiness {
public:
   enum AxesMode {
      kt_axes,  // exclusive kt axes
      ca_axes,  // exclusive ca axes
      antikt_0p2_axes,  // inclusive hardest axes with antikt-0.2
      min_axes, // axes that minimize N-subjettiness (100 passes by default)
      manual_axes, // set your own axes with setAxes()
      onepass_kt_axes, // one-pass minimization from kt starting point
      onepass_ca_axes, // one-pass minimization from ca starting point
      onepass_antikt_0p2_axes,  // one-pass minimization from antikt-0.2 starting point 
      onepass_manual_axes  // one-pass minimization from manual starting point
   };

private:
   AxesMode _axes;    //Axes mode choice
   NsubParameters _paraNsub;  //Parameters for Njettiness
   KmeansParameters _paraKmeans;  //Parameters for Minimization Procedure (set by NsubAxesMode automatically, but can change manually if desired)

   std::vector<fastjet::PseudoJet> _currentAxes;
   
   void establishAxes(unsigned n_jets, const std::vector <fastjet::PseudoJet> & inputs);
   
public:
   Njettiness(AxesMode axes, NsubParameters paraNsub);
   
   void setParaKmeans(KmeansParameters newPara) {_paraKmeans = newPara;}
   void setParaNsub(NsubParameters newPara) {_paraNsub = newPara;}
   
   // setAxes for Manual mode
   void setAxes(std::vector<fastjet::PseudoJet> myAxes) {
      assert((_axes == manual_axes) || (_axes == onepass_manual_axes));
      _currentAxes = myAxes;
   }
   
   // The value of N-subjettiness
   double getTau(unsigned n_jets, const std::vector<fastjet::PseudoJet> & inputJets) {
      if (inputJets.size() <= n_jets) {
         _currentAxes = inputJets;
         _currentAxes.resize(n_jets,fastjet::PseudoJet(0.0,0.0,0.0,0.0));
         return 0.0;
      }
      establishAxes(n_jets, inputJets);  // sets current Axes
      return TauValue(inputJets,_currentAxes,_paraNsub);
   }
   
   // get axes used by getTau.
   std::vector<fastjet::PseudoJet> currentAxes() {
      return _currentAxes;
   }
   
   // partition inputs by Voronoi (each vector stores indices corresponding to inputJets)
   std::vector<std::list<int> > getPartition(const std::vector<fastjet::PseudoJet> & inputJets);

   // partition inputs by Voronoi
   std::vector<fastjet::PseudoJet> getJets(const std::vector<fastjet::PseudoJet> & inputJets);

};

#endif

