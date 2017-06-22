//
//Aditya Verma
//June 15, 2017
//This code creates a histigram of Z=pttrack/ptjet which is the fragmentation analysis
//
//
//
//
//
//

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {

  /// @brief Add a short analysis description here
  class delphitest_analysis : public Analysis {

  public:

    /// Constructor
    delphitest_analysis()
      : Analysis("delphitest_analysis")
    {    }

    void init() {

      //! finalstate particles and the jets with the algorithm and radius definition
      FinalState fs(-5.0, 5.0, 0.150*GeV);
      addProjection(fs, "FS");
      FastJets fj(fs, FastJets::ANTIKT, 0.4);
      addProjection(fj, "Jets");

   

      ChargedFinalState cfs(-5.0, 5.0, 0.100*GeV);
      addProjection(cfs, "CFS");

      //! initialize variables
      _ptCut = 30.0;
      _nJets = 0.0;
      
      //! initialize histograms
      _h_jetpT = bookHisto1D("jetpT", 100, 0, 300);
      _h_Z = bookHisto1D("Z", 10, 0, 1);
      _h_Aj = bookHisto1D("Aj", 100, 0, 1);
      _h_delphi = bookHisto1D("Delta_Phi",100, 0, 5);
      _h_Iaa = bookHisto1D("Iaa",100, 0, 150);

    

      
    }

    //intialize global variables
    void analyze(const Event& event) {

      //! event cross section weight 
      const double weight = event.weight();
      //! Jet collection of importance 
      const FastJets& Ajets = applyProjection<FastJets>(event, "Jets");
      Cut cuts = Cuts::etaIn(-2, 2) & (Cuts::pT > _ptCut*GeV);
      const Jets ajets = Ajets.jetsByPt(cuts);
      if(ajets.size()<2)
	vetoEvent;

      
      const ChargedFinalState tracks= applyProjection<ChargedFinalState>(event, "CFS");

      
      //! Fill histograms 
      //! Loop over jets in the event
      foreach ( PseudoJet jet, ajets ){
	_h_jetpT->fill(jet.pt(), weight);

	if(ajets[0].pt()>30.0 && _nJets!=0){
	  _h_Iaa->fill(jet.pt(), weight);
	}
	  
	
	int nTracks =0;
	_nJets+=weight;//Jet counter

	
	foreach ( Particle track, tracks.particles() ){
	  
	  double delR = deltaR(track.eta(),track.phi(MINUSPI_PLUSPI),jet.eta(),jet.phi_std());
	  if(delR<4){
	    nTracks++;
	    _h_Z->fill(track.pt()*cos(delR)/jet.pt(), weight);
	    
	  }
	}
	
      }

      //! aj 
      _h_Aj->fill((ajets[0].pt()-ajets[1].pt())/(ajets[0].pt()+ajets[1].pt()), weight);
      _h_delphi->fill(ajets[0].phi()-ajets[1].phi(),weight);
    
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_jetpT, crossSection()/picobarn/sumOfWeights());
      scale(_h_Z, _nJets);
      scale(_h_Aj, 1/sumOfWeights());
      scale(_h_delphi, 1/sumOfWeights());
      scale(_h_Iaa,crossSection()/picobarn/sumOfWeights());
    }

    //! Declare variables 
  private:
    Histo1DPtr _h_jetpT;
    Histo1DPtr _h_Z;
    Histo1DPtr _h_Aj;
    Histo1DPtr _h_delphi;
    Histo1DPtr _h_Iaa;

  protected:
    double _ptCut;
    double _nJets;//Jet Counter, Track Counter
    
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(delphitest_analysis);


}
