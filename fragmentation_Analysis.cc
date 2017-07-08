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
  class fragmentation_Analysis : public Analysis {
  public:

    /// Constructor
    fragmentation_Analysis()
      : Analysis("fragmentation_Analysis")
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
      _h_njets = bookHisto1D("njets", 10, 0, 10);
      _h_jetpT = bookHisto1D("jetpT", 100, 0, 300);
      _h_jeteta = bookHisto1D("jeteta", 100, -2.5, 2.5);
      _h_jetphi = bookHisto1D("jetphi", 100, -M_PI, M_PI);

      _h_ntracksPJ = bookHisto1D("nTracksPerJet", 10, 0, 50);//
      _h_trackpT = bookHisto1D("trackpT", 100, 0, 300);//
      _h_tracketa = bookHisto1D("tracketa", 100, -2.5, 2.5);//
      _h_trackphi = bookHisto1D("trackphi", 100, -M_PI, M_PI);//
      _h_Z = bookHisto1D("Z", 10, 0, 1);

    }

    //intialize global variables
    void analyze(const Event& event) {

      //! event cross section weight 
      const double weight = event.weight();
      //! Jet collection of importance 
      const FastJets& Ajets = applyProjection<FastJets>(event, "Jets");
      Cut cuts = Cuts::etaIn(-2, 2) & (Cuts::pT > _ptCut*GeV);
      const Jets ajets = Ajets.jetsByPt(cuts);

      const ChargedFinalState tracks= applyProjection<ChargedFinalState>(event, "CFS");
      
      //! Fill histograms 
      _h_njets->fill(ajets.size(), weight);
      //! Loop over jets in the event
      foreach ( PseudoJet jet, ajets ){
	_h_jetpT->fill(jet.pt(), weight);
	_h_jeteta->fill(jet.eta(), weight);
	_h_jetphi->fill(jet.phi_std(), weight);
	
	int nTracks =0;
	_nJets+=weight;//Jet counter

	foreach ( Particle track, tracks.particles() ){

	  
	  double delR = deltaR(track.eta(),track.phi(MINUSPI_PLUSPI),jet.eta(),jet.phi_std());
	  if(delR<4){
	    nTracks++;
	    _h_trackpT->fill(track.pt(), weight);
	    _h_tracketa->fill(track.eta(), weight);
	    _h_trackphi->fill(track.phi(MINUSPI_PLUSPI), weight);
	    _h_Z->fill(track.pt()*cos(delR)/jet.pt(), weight);
	    
	  }
	}
		_h_ntracksPJ->fill(nTracks, weight);
	
      }
      	
      
    }

    /// Normalise histograms etc., after the run
     void finalize() {
      scale(_h_njets, crossSection()/picobarn/sumOfWeights());
      scale(_h_jetpT, crossSection()/picobarn/sumOfWeights());
      scale(_h_jeteta, crossSection()/picobarn/sumOfWeights());
      scale(_h_jetphi, crossSection()/picobarn/sumOfWeights());

      scale(_h_ntracksPJ, _nJets);//Scaling it to number of jets
      scale(_h_trackpT, _nJets);
      scale(_h_tracketa, _nJets);
      scale(_h_trackphi, _nJets);
      
    }

    //! Declare variables 
  private:
    Histo1DPtr _h_njets;
    Histo1DPtr _h_jetpT;
    Histo1DPtr _h_jeteta;
    Histo1DPtr _h_jetphi;

    Histo1DPtr _h_ntracksPJ;
    Histo1DPtr _h_trackpT;
    Histo1DPtr _h_tracketa;
    Histo1DPtr _h_trackphi;

    Histo1DPtr _h_Z;

  protected:
    double _ptCut;
    double _nJets;//Jet Counter, Track Counter
    
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(fragmentation_Analysis);


}
