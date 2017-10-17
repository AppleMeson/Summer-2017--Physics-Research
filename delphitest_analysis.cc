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
      FinalState fs(-5.0, 5.0, 0.0*GeV);
      addProjection(fs, "FS");
      FastJets fj(fs, FastJets::ANTIKT, 0.4);
      addProjection(fj, "Jets");

   
      //! charged particles (tracks) with the appropriate cuts
      ChargedFinalState cfs(-5.0, 5.0, 1.0*GeV);
      addProjection(cfs, "CFS");

      //! initialize variables
      _ptCut = 10.0;// change to something else... 2 doesn't work?
      _nJets = 0.0;
      _nevts = 0.0;
      
      //! initialize histograms
      _h_jetpT = bookHisto1D("jetpT", 20, 0, 100);
      _h_Z = bookHisto1D("Z", 10, 0, 1);
      _h_Aj = bookHisto1D("Aj", 20, 0, 1);
      _h_delphi = bookHisto1D("Delta_Phi", 60, -5, 5);
      _h_Iaa = bookHisto1D("Iaa",20, 0, 100);

    

      
    }

    //intialize global variables
    void analyze(const Event& event) {

      int nJetpE=0;
      //! event cross section weight 
      const double weight = event.weight(); 
      const ChargedFinalState tracks= applyProjection<ChargedFinalState>(event, "CFS");

      //! Jet collection of importance 
      const FastJets& Ajets = applyProjection<FastJets>(event, "Jets"); //
      Cut cuts = Cuts::etaIn(-1, 1) & (Cuts::pT > _ptCut*GeV);
      const Jets ajets = Ajets.jetsByPt(cuts);
      if(ajets.size()==0)
	vetoEvent;

      //! cut for leading jet, used in plotting the Iaa spectrum
      bool leadjetcut = false;
      if(ajets.at(0).pt() > 30)
	leadjetcut = true;
      //! dijet cut used for delphi and Aj
      bool dijetcut = false;
      if(ajets.size()>=2 && ajets.at(0).pt() > 20 && ajets.at(1).pt() > 10)//20 n 10
	{
	dijetcut = true;
	_nevts+=weight;
      }
      
      
      //! Fill histograms 
      //! Loop over jets in the event
      foreach ( PseudoJet jet, ajets ){ 

	_h_jetpT->fill(jet.pt(), weight);

	if(leadjetcut && nJetpE>0)
	  _h_Iaa->fill(jet.pt(), weight);
	
	int nTracks =0;
	_nJets+=weight;//Jet counter

	
	foreach ( Particle track, tracks.particles() ){

	  //! association cut
	  double delR = deltaR(track.eta(),track.phi(MINUSPI_PLUSPI),jet.eta(),jet.phi_std());

	  if(delR<=0.4){
	    nTracks++;
	    //! Z
	    _h_Z->fill(track.pt()*cos(delR)/jet.pt(), weight);
	        
	  }
	  else
	    //cout<<"No"
	       ;
	  
	}//! particle loop
	
	nJetpE++;
      }//! jet loop

      //! aj
      if(dijetcut){

	if(abs(ajets[0].phi()-ajets[1].phi()-PI)<0.4){
	_h_Aj->fill(((ajets[0].pt()-ajets[1].pt())/(ajets[0].pt()+ajets[1].pt())), weight);// put in a phi cut for phi1-phi2-pi<0.4
	}
	
	_h_delphi->fill(abs(ajets[0].phi()-ajets[1].phi()), weight);
      }

      
    }

    

    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_jetpT, crossSection()/picobarn/sumOfWeights());
      scale(_h_Z, 1./_nJets);
      scale(_h_Aj, 1./_nevts);
      scale(_h_delphi, 1./_nevts);
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
    double _nevts;
    
  };


  // The hook for the plugin syst
  DECLARE_RIVET_PLUGIN(delphitest_analysis);


}
