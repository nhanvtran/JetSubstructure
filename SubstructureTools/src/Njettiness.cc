#include"JetSubstructure/SubstructureTools/interface/Njettiness.hh"


// Parameters that define Njettiness
///////
//
// Functions for minimization.
// TODO:  Wrap these in N-subjettiness class
//
///////

inline double sq(double x) {return x*x;}

// Calculates distance between two points in rapidity-azimuth plane
double DistanceSq(double rap1, double phi1, double rap2, double phi2) {
   double distRap = rap1-rap2;
   double distPhi = std::fabs(phi1-phi2);
   if (distPhi > M_PI) {distPhi = 2.0*M_PI - distPhi;}
   return sq(distRap) + sq(distPhi);
}

inline double Distance(double rap1, double phi1, double rap2, double phi2) {
   return std::sqrt(DistanceSq(rap1, phi1, rap2, phi2));
}


// Given starting axes, update to find better axes
template <int N>
std::vector<LightLikeAxis> UpdateAxesFast(const std::vector <LightLikeAxis> & old_axes, 
                                  const std::vector <fastjet::PseudoJet> & inputJets,
                                  NsubParameters paraNsub, double precision) {
   assert(old_axes.size() == N);
   
   // some storage, declared static to save allocation/re-allocation costs
   static LightLikeAxis new_axes[N];
   static fastjet::PseudoJet new_jets[N];
   for (int n = 0; n < N; ++n) {
      new_axes[n].reset(0.0,0.0,0.0,0.0);
#ifdef FASTJET2
      new_jets[n].reset(0.0,0.0,0.0,0.0);
#else
      // use cheaper reset if available
      new_jets[n].reset_momentum(0.0,0.0,0.0,0.0);
#endif
   }

   
   double beta = paraNsub.beta();
   double Rcutoff = paraNsub.Rcutoff();
   
   /////////////// Assignment Step //////////////////////////////////////////////////////////
   std::vector<int> assignment_index(inputJets.size()); 
   int k_assign = -1;
   
   for (unsigned i = 0; i < inputJets.size(); i++){
      double smallestDist = 1000000.0;
      for (int k = 0; k < N; k++) {
         double thisDist = DistanceSq(inputJets[i].rap(),inputJets[i].phi(),old_axes[k].rap(),old_axes[k].phi());
         if (thisDist < smallestDist) {
            smallestDist = thisDist;
            k_assign = k;
         }
      }
      if (smallestDist > sq(Rcutoff)) {k_assign = -1;}
      assignment_index[i] = k_assign;
   }
   
   //////////////// Update Step /////////////////////////////////////////////////////////////
   double distPhi, old_dist;
   for (unsigned i = 0; i < inputJets.size(); i++) {
      int old_jet_i = assignment_index[i];
      if (old_jet_i == -1) {continue;}

      const fastjet::PseudoJet& inputJet_i = inputJets[i];
      LightLikeAxis& new_axis_i = new_axes[old_jet_i];
      double inputPhi_i = inputJet_i.phi();
      double inputRap_i = inputJet_i.rap();
            
      // optimize pow() call
      // add noise (the precision term) to make sure we don't divide by zero
      if (beta == 1.0) {
         double DR = std::sqrt(sq(precision) + DistanceSq(inputRap_i, inputPhi_i, old_axes[old_jet_i].rap(), old_axes[old_jet_i].phi()));
         old_dist = 1.0/DR;
      } else if (beta == 2.0) {
         old_dist = 1.0;
      } else if (beta == 0.0) {
         double DRSq = sq(precision) + DistanceSq(inputRap_i, inputPhi_i, old_axes[old_jet_i].rap(), old_axes[old_jet_i].phi());
         old_dist = 1.0/DRSq;
      } else {
         old_dist = sq(precision) + DistanceSq(inputRap_i, inputPhi_i, old_axes[old_jet_i].rap(), old_axes[old_jet_i].phi());
         old_dist = std::pow(old_dist, (0.5*beta-1.0));
      }
      
      // rapidity sum
      new_axis_i.set_rap(new_axis_i.rap() + inputJet_i.perp() * inputRap_i * old_dist);
      // phi sum
      distPhi = inputPhi_i - old_axes[old_jet_i].phi();
      if (fabs(distPhi) <= M_PI){
         new_axis_i.set_phi( new_axis_i.phi() + inputJet_i.perp() * inputPhi_i * old_dist );
      } else if (distPhi > M_PI) {
         new_axis_i.set_phi( new_axis_i.phi() + inputJet_i.perp() * (-2*M_PI + inputPhi_i) * old_dist );
      } else if (distPhi < -M_PI) {
         new_axis_i.set_phi( new_axis_i.phi() + inputJet_i.perp() * (+2*M_PI + inputPhi_i) * old_dist );
      }
      // weights sum
      new_axis_i.set_weight( new_axis_i.weight() + inputJet_i.perp() * old_dist );
      // momentum magnitude sum
      new_jets[old_jet_i] += inputJet_i;
   }
   // normalize sums
   for (int k = 0; k < N; k++) {
      if (new_axes[k].weight() == 0) {
         // no particles were closest to this axis!  Return to old axis instead of (0,0,0,0)
         new_axes[k] = old_axes[k];
      } else {
         new_axes[k].set_rap( new_axes[k].rap() / new_axes[k].weight() );
         new_axes[k].set_phi( new_axes[k].phi() / new_axes[k].weight() );
         new_axes[k].set_phi( std::fmod(new_axes[k].phi() + 2*M_PI, 2*M_PI) );
         new_axes[k].set_mom( std::sqrt(new_jets[k].modp2()) );
      }
   }
   std::vector<LightLikeAxis> new_axes_vec(N);
   for (unsigned k = 0; k < N; ++k) new_axes_vec[k] = new_axes[k];
   return new_axes_vec;
}


// Given starting axes, update to find better axes
// (This is just a wrapper for the templated version above.)
std::vector<LightLikeAxis> UpdateAxes(const std::vector <LightLikeAxis> & old_axes, 
                                      const std::vector <fastjet::PseudoJet> & inputJets, NsubParameters paraNsub, double precision) {
   int N = old_axes.size();
   switch (N) {
      case 1: return UpdateAxesFast<1>(old_axes, inputJets, paraNsub, precision);
      case 2: return UpdateAxesFast<2>(old_axes, inputJets, paraNsub, precision);
      case 3: return UpdateAxesFast<3>(old_axes, inputJets, paraNsub, precision);
      case 4: return UpdateAxesFast<4>(old_axes, inputJets, paraNsub, precision);
      case 5: return UpdateAxesFast<5>(old_axes, inputJets, paraNsub, precision);
      case 6: return UpdateAxesFast<6>(old_axes, inputJets, paraNsub, precision);
      case 7: return UpdateAxesFast<7>(old_axes, inputJets, paraNsub, precision);
      case 8: return UpdateAxesFast<8>(old_axes, inputJets, paraNsub, precision);
      case 9: return UpdateAxesFast<9>(old_axes, inputJets, paraNsub, precision);
      case 10: return UpdateAxesFast<10>(old_axes, inputJets, paraNsub, precision);
      case 11: return UpdateAxesFast<11>(old_axes, inputJets, paraNsub, precision);
      case 12: return UpdateAxesFast<12>(old_axes, inputJets, paraNsub, precision);
      case 13: return UpdateAxesFast<13>(old_axes, inputJets, paraNsub, precision);
      case 14: return UpdateAxesFast<14>(old_axes, inputJets, paraNsub, precision);
      case 15: return UpdateAxesFast<15>(old_axes, inputJets, paraNsub, precision);
      case 16: return UpdateAxesFast<16>(old_axes, inputJets, paraNsub, precision);
      case 17: return UpdateAxesFast<17>(old_axes, inputJets, paraNsub, precision);
      case 18: return UpdateAxesFast<18>(old_axes, inputJets, paraNsub, precision);
      case 19: return UpdateAxesFast<19>(old_axes, inputJets, paraNsub, precision);
      case 20: return UpdateAxesFast<20>(old_axes, inputJets, paraNsub, precision);
      default: std::cout << "N-jettiness is hard-coded to only allow up to 20 jets!" << std::endl;
         return std::vector<LightLikeAxis>();
   }

}


// Go from internal LightLikeAxis to PseudoJet
// TODO:  Make part of LightLikeAxis class.
std::vector<fastjet::PseudoJet> ConvertToPseudoJet(const std::vector <LightLikeAxis>& axes) {
   
   int n_jets = axes.size();
   
   double px, py, pz, E;
   std::vector<fastjet::PseudoJet> FourVecJets;
   for (int k = 0; k < n_jets; k++) {
      E = axes[k].mom();
      pz = (std::exp(2.0*axes[k].rap()) - 1.0) / (std::exp(2.0*axes[k].rap()) + 1.0) * E;
      px = std::cos(axes[k].phi()) * std::sqrt( std::pow(E,2) - std::pow(pz,2) );
      py = std::sin(axes[k].phi()) * std::sqrt( std::pow(E,2) - std::pow(pz,2) );
      fastjet::PseudoJet temp = fastjet::PseudoJet(px,py,pz,E);
      FourVecJets.push_back(temp);
   }
   return FourVecJets;
}

// N-subjettiness pieces
std::vector<double> ConstituentTauValue(const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& axes, const NsubParameters& paraNsub) {// Returns the sub-tau values, i.e. a std::vector of the contributions to tau_N of each Voronoi region (or region within R_0)
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
} 

// N-subjettiness values
double TauValue(const std::vector <fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes,const NsubParameters& paraNsub) {// Calculates tau_N
   std::vector<double> tau_vec = ConstituentTauValue(particles, axes, paraNsub);
   double tau = 0.0;
   for (unsigned j = 0; j < tau_vec.size(); j++) {tau += tau_vec[j];}
   return tau;
}

// Get exclusive kT subjets
std::vector<fastjet::PseudoJet> GetKTAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputJets) {
   fastjet::JetDefinition jet_def = fastjet::JetDefinition(fastjet::kt_algorithm,M_PI/2.0,fastjet::E_scheme,fastjet::Best);
   fastjet::ClusterSequence jet_clust_seq(inputJets, jet_def);
   return jet_clust_seq.exclusive_jets(n_jets);
}

// Get exclusive CA subjets
std::vector<fastjet::PseudoJet> GetCAAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputJets) {
   fastjet::JetDefinition jet_def = fastjet::JetDefinition(fastjet::cambridge_algorithm,M_PI/2.0,fastjet::E_scheme,fastjet::Best);
   fastjet::ClusterSequence jet_clust_seq(inputJets, jet_def);
   return jet_clust_seq.exclusive_jets(n_jets);
}

// Get inclusive anti kT hardest subjets
std::vector<fastjet::PseudoJet> GetAntiKTAxes(int n_jets, double R0, const std::vector <fastjet::PseudoJet> & inputJets) {
   fastjet::JetDefinition jet_def = fastjet::JetDefinition(fastjet::antikt_algorithm,R0,fastjet::E_scheme,fastjet::Best);
   fastjet::ClusterSequence jet_clust_seq(inputJets, jet_def);
   std::vector<fastjet::PseudoJet> myJets = sorted_by_pt(jet_clust_seq.inclusive_jets());
   myJets.resize(n_jets);  // only keep n hardest
   return myJets;
}


// Get minimization axes
std::vector<fastjet::PseudoJet> GetMinimumAxes(const std::vector <fastjet::PseudoJet> & seedAxes, const std::vector <fastjet::PseudoJet> & inputJets, KmeansParameters para, 
                                          NsubParameters paraNsub) {
   int n_jets = seedAxes.size();
   double noise = 0, tau = 10000.0, tau_tmp, cmp;
   std::vector< LightLikeAxis > new_axes(n_jets, LightLikeAxis(0,0,0,0)), old_axes(n_jets, LightLikeAxis(0,0,0,0));
   std::vector<fastjet::PseudoJet> tmp_min_axes, min_axes;
   for (int l = 0; l < para.n_iterations(); l++) { // Do minimization procedure multiple times
      // Add noise to guess for the axes
      for (int k = 0; k < n_jets; k++) {
         if ( 0 == l ) { // Don't add noise on first pass
            old_axes[k].set_rap( seedAxes[k].rap() + noise );
            old_axes[k].set_phi( seedAxes[k].phi() + noise );
         } else {
            noise = ((double)rand()/(double)RAND_MAX) * para.noise_range() * 2 - para.noise_range();
            old_axes[k].set_rap( seedAxes[k].rap() + noise );
            noise = ((double)rand()/(double)RAND_MAX) * para.noise_range() * 2 - para.noise_range();
            old_axes[k].set_phi( seedAxes[k].phi() + noise );
         }
      }
      cmp = 100.0; int h = 0;
      while (cmp > para.precision() && h < para.halt()) { // Keep updating axes until near-convergence or too many update steps
         cmp = 0.0; h++;
         new_axes = UpdateAxes(old_axes, inputJets, paraNsub, para.precision()); // Update axes
         for (int k = 0; k < n_jets; k++) {
            cmp += Distance(new_axes[k].rap(),new_axes[k].phi(),old_axes[k].rap(),old_axes[k].phi());
         }
         cmp = cmp / ((double) n_jets);
         old_axes = new_axes;
      }
      tmp_min_axes = ConvertToPseudoJet(old_axes); // Convert axes directions into four-std::vectors
      tau_tmp = TauValue(inputJets, tmp_min_axes,paraNsub); 
      if (tau_tmp < tau) {tau = tau_tmp; min_axes = tmp_min_axes;} // Keep axes and tau only if they are best so far
   }	
   return min_axes;
}

//Use NsubAxesMode to pick which type of axes to use
void Njettiness::establishAxes(unsigned int n_jets, const std::vector <fastjet::PseudoJet> & inputs) {
   switch (_axes) {
      case kt_axes:
         _currentAxes = GetKTAxes(n_jets,inputs);
         break;
      case ca_axes:
         _currentAxes = GetCAAxes(n_jets,inputs);
         break;
      case antikt_0p2_axes:
         _currentAxes = GetAntiKTAxes(n_jets,0.2,inputs);
         break;
      case onepass_kt_axes:
      case min_axes:
         _currentAxes = GetKTAxes(n_jets,inputs);
         _currentAxes = GetMinimumAxes(_currentAxes, inputs, _paraKmeans, _paraNsub);
         break;
      case onepass_ca_axes:
         _currentAxes = GetCAAxes(n_jets,inputs);
         _currentAxes = GetMinimumAxes(_currentAxes, inputs, _paraKmeans, _paraNsub);
         break;
      case onepass_antikt_0p2_axes:
         _currentAxes = GetAntiKTAxes(n_jets,0.2,inputs);
         _currentAxes = GetMinimumAxes(_currentAxes, inputs, _paraKmeans, _paraNsub);
         break;
      case onepass_manual_axes:
         assert(_currentAxes.size() == n_jets);
         _currentAxes = GetMinimumAxes(_currentAxes, inputs, _paraKmeans, _paraNsub);
         break;
      case manual_axes:
         assert(_currentAxes.size() == n_jets);
         break;
      default:
         assert(false);
         break;
   }      
}

//Constructor sets KmeansParameters from NsubAxesMode input
Njettiness::Njettiness(AxesMode axes, NsubParameters paraNsub) : _axes(axes), _paraNsub(paraNsub), _paraKmeans() {
   switch (_axes) {
      case kt_axes:
      case ca_axes:
      case antikt_0p2_axes:
         _paraKmeans = KmeansParameters(0,0.0,0,0.0);
         break;
      case onepass_kt_axes:
      case onepass_ca_axes:
      case onepass_antikt_0p2_axes:
      case onepass_manual_axes:
         _paraKmeans = KmeansParameters(1,0.0001,1000,0.8);
         break;
      case min_axes:
         _paraKmeans = KmeansParameters(100,0.0001,1000,0.8);
         break;
      case manual_axes:
         _paraKmeans = KmeansParameters(0,0.0,0,0.0);
         break;
      default:
         assert(false);
         break;
   }
}


// Partition a list of particles according to which N-jettiness axis they are closest to.
// Return a vector of length _currentAxes.size() (which should be N).
// Each vector element is a list of ints corresponding to the indices in
// particles of the particles belonging to that jet.
std::vector<std::list<int> > Njettiness::getPartition(const std::vector<fastjet::PseudoJet> & particles) {
   double Rcutoff = _paraNsub.Rcutoff();
   
   std::vector<std::list<int> > partitions(_currentAxes.size());

   int j_min = -1;
   for (unsigned i = 0; i < particles.size(); i++) {
      // find minimum distance
      double minR = 10000.0; // large number
      for (unsigned j = 0; j < _currentAxes.size(); j++) {
         double tempR = std::sqrt(particles[i].squared_distance(_currentAxes[j])); // delta R distance
         if (tempR < minR) {
            minR = tempR;
            j_min = j;
         }
      }
      if (minR > Rcutoff) {
         // do nothing
      } else {
         partitions[j_min].push_back(i);
      }
   }
   return partitions;
}


// Having found axes, assign each particle in particles to an axis, and return a set of jets.
// Each jet is the sum of particles closest to an axis (Njet = Naxes).
std::vector<fastjet::PseudoJet> Njettiness::getJets(const std::vector<fastjet::PseudoJet> & particles) {
   
   std::vector<fastjet::PseudoJet> jets(_currentAxes.size());

   std::vector<std::list<int> > partition = getPartition(particles);
   for (unsigned j = 0; j < partition.size(); ++j) {
      std::list<int>::const_iterator it, itE;
      for (it = partition[j].begin(), itE = partition[j].end(); it != itE; ++it) {
         jets[j] += particles[*it];
      }
   }
   return jets;
}

