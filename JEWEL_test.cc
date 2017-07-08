// -*- C++ -*-

/*
  
  Raghav Kunnawalkam Elayavalli
  Rutgers University
  June 9th, 2017 
  
  Example RIVET Analysis to learn about rivet and analyze hepMC files 

  The way RIVET handles an input hepMC file is the following: 
  0. First you create the analysis as follows: 
     $ rivet-mkanalysis Analysis_Name 
     That creates Analysis_Name.cc, Analysis_Name.plot and Analysis_Name.info 
  1. You define a specific analysis class to do what you want 
  2. Global variables and functions/methods are declared as either private/protected 
  3. The init function is run in the beginning. so initialize your global variables there 
     3a. finalstate particles and jet collections, histograms etc. are defined here
  4. the analyze function is run per event 
     4a. Get the weight and event specific variables 
     4b. Loop over jets and fill histograms 
  5. Finalyze function runs at the end of the events 
     5a. Use that to normalize histograms and other things that you want to do at the end 

 */



#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  //! Analysis class : inherits from public Analysis 
  class JEWEL_test : public Analysis {

  public:
    //! Constructor 
    JEWEL_test()
      : Analysis("JEWEL_test")
    {    }

    //! Initialize global variables 
    void init() {

      //! finalstate particles and the jets with the algorithm and radius definition
      FinalState fs(-5.0, 5.0, 0.150*GeV);
      addProjection(fs, "FS");
      FastJets fj(fs, FastJets::ANTIKT, 0.4);
      addProjection(fj, "Jets");

      //! initialize variables
      _ptCut = 30.0;

      //! initialize histograms       
      _h_njets = bookHisto1D("njets", 10, 0, 10);
      _h_jetpT = bookHisto1D("jetpT", 100, 0, 300);
      _h_jeteta = bookHisto1D("jeteta", 100, -5, 5);
      _h_jetphi = bookHisto1D("jetphi", 100, -M_PI, M_PI);

    }

    //! Analyze function that runs per event 
    void analyze(const Event& event) {

      //! event cross section weight 
      const double weight = event.weight();
      //! Jet collection of importance 
      const FastJets& Ajets = applyProjection<FastJets>(event, "Jets");
      Cut cuts = Cuts::etaIn(-2, 2) & (Cuts::pT > _ptCut*GeV);
      const Jets ajets = Ajets.jetsByPt(cuts);

      //! Fill histograms 
      _h_njets->fill(ajets.size(), weight);
      //! Loop over jets in the event
      foreach ( PseudoJet jet, ajets ){
	_h_jetpT->fill(jet.pt(), weight);
	_h_jeteta->fill(jet.eta(), weight);
	_h_jetphi->fill(jet.phi_std(), weight);	
      }
      
    }

    //! Runs after all the events. normalize histograms etc.
    void finalize() {
      scale(_h_njets, crossSection()/picobarn/sumOfWeights());
      scale(_h_jetpT, crossSection()/picobarn/sumOfWeights());
      scale(_h_jeteta, crossSection()/picobarn/sumOfWeights());
      scale(_h_jetphi, crossSection()/picobarn/sumOfWeights());
    }

    //! Declare variables 
  private:
    Histo1DPtr _h_njets;
    Histo1DPtr _h_jetpT;
    Histo1DPtr _h_jeteta;
    Histo1DPtr _h_jetphi;

  protected:
    double _ptCut;
    
  };

  DECLARE_RIVET_PLUGIN(JEWEL_test);
}
