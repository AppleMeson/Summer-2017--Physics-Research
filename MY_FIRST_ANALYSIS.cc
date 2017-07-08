// -*- C++ -*-
/*
Aditya Verma
Rutgers University
June 9th, 2017

To analyse (simple analysis) jets in pp files that we just created
*/

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MY_FIRST_ANALYSIS : public Analysis {
  public:

    /// Constructor
    MY_FIRST_ANALYSIS()
      : Analysis("MY_FIRST_ANALYSIS")
    {   }


    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs(-5.0, 5.0, 0.150*GeV);
      addProjection(fs, "FS");
      FastJets fj(fs, FastJets::ANTIKT, 0.4);
      addProjection(fj, "Jets");

      //intializing variables
      _ptCut = 30.0;

      //initializing histograms
      _h_njets = bookHisto1D("njets",10,0,10);
      _h_jetpT = bookHisto1D("jetpT", 100, 0, 300);
      _h_jeteta = bookHisto1D("jeteta", 100, -5, 5);
      _h_jetphi = bookHisto1D("jetphi", 100, -M_PI, M_PI);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // event cross section weight
      const double weight = event.weight();
      // Jet collection of importance
      const FastJets& Ajets = applyProjection<FastJets>(event, "Jets");
      Cut cuts = Cuts::etaIn(-2,2) & (Cuts::pT > _ptCut*GeV);
      const Jets ajets = Ajets.jetsByPt(cuts);
      // Fill histograms

      _h_njets->fill(ajets.size(), weight);
      //Loop pver each jet
      foreach (PseudoJet jet, ajets){
	_h_jetpT->fill(jet.pt(), weight);
	_h_jeteta->fill(jet.eta(), weight);
	_h_jetphi->fill(jet.phi_std(), weight);
      }
      

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_njets, crossSection()/picobarn/sumOfWeights());
      scale(_h_jetpT, crossSection()/picobarn/sumOfWeights());
      scale(_h_jeteta, crossSection()/picobarn/sumOfWeights());
      scale(_h_jetphi, crossSection()/picobarn/sumOfWeights());

    }

    //Declare variables

  private:
    Histo1DPtr _h_njets;
    Histo1DPtr _h_jetpT;
    Histo1DPtr _h_jetphi;
    Histo1DPtr _h_jeteta;
  protected:
    double _ptCut;


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MY_FIRST_ANALYSIS);


}
