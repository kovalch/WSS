#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include <TrackingTools/TransientTrack/interface/TransientTrackBuilder.h>
int main(int argc, char* argv[]) 
{
  // define what muon you are using; this is necessary as FWLite is not 
  // capable of reading edm::Views
//  using reco::Muon;
  
  // ----------------------------------------------------------------------
  // First Part: 
  //
  //  * enable the AutoLibraryLoader 
  //  * book the histograms of interest 
  //  * open the input file
  // ----------------------------------------------------------------------

  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  // initialize command line parser
  optutl::CommandLineParser parser ("Analyze FWLite Histograms");

  // set defaults
  parser.integerValue ("maxEvents"  ) = 1000;
  parser.integerValue ("outputEvery") =   10;
  parser.stringValue  ("outputFile" ) = "analyzeFWLiteHistograms.root";

  // parse arguments
  parser.parseArguments (argc, argv);
  int maxEvents_ = parser.integerValue("maxEvents");
//  unsigned int outputEvery_ = parser.integerValue("outputEvery");
  std::string outputFile_ = parser.stringValue("outputFile");
  std::vector<std::string> inputFiles_ = parser.stringVector("inputFiles");

  // book a set of histograms
  fwlite::TFileService fs = fwlite::TFileService(outputFile_.c_str());
  TFileDirectory dir = fs.mkdir("");

  TH1F* muonPx_  = dir.make<TH1F>("muonPx"  , "px"  ,   100,   0., 300.);
  TH1F* muonPy_  = dir.make<TH1F>("muonPy"  , "py"  ,   100,   0., 300.);
  TH1F* muonEt_  = dir.make<TH1F>("muonEt"  , "et"  ,   100,   0., 300.);
  TH1F* muonPt_  = dir.make<TH1F>("muonPt"  , "pt"  ,   100,   0., 300.);
  TH1F* muonEta_ = dir.make<TH1F>("muonEta" , "eta" ,   100,  -3.,   3.);
  TH1F* muonPhi_ = dir.make<TH1F>("muonPhi" , "phi" ,   100,  -5.,   5.);  
  TH1F* mumuMass_= dir.make<TH1F>("mumuMass", "mass",    90,  30.,  120.);
  TH1F* electronPt_  = dir.make<TH1F>("electronPt"  , "pt"  ,   100,   0., 300.);
  TH1F* electronEta_ = dir.make<TH1F>("electronEta" , "eta" ,   100,  -3.,   3.);
  TH1F* electronPhi_ = dir.make<TH1F>("electronPhi" , "phi" ,   100,  -5.,   5.);
  TH1F* jetPt_  = dir.make<TH1F>("jetPt"  , "pt"  ,   100,   0., 300.);
  TH1F* jetEta_ = dir.make<TH1F>("jetEta" , "eta" ,   100,  -3.,   3.);
  TH1F* jetPhi_ = dir.make<TH1F>("jetPhi" , "phi" ,   100,  -5.,   5.);
  TH1F* metPt_  = dir.make<TH1F>("metPt"  , "pt"  ,   100,   0., 300.);
  TH1F* metEta_ = dir.make<TH1F>("metEta" , "eta" ,   100,  -3.,   3.);
  TH1F* metPhi_ = dir.make<TH1F>("metPhi" , "phi" ,   100,  -5.,   5.);

  TH1F* photonPt_  = dir.make<TH1F>("photonPt"  , "pt"  ,   100,   0., 300.);
  TH1F* photonEta_ = dir.make<TH1F>("photonEta" , "eta" ,   100,  -3.,   3.);
  TH1F* photonPhi_ = dir.make<TH1F>("photonPhi" , "phi" ,   100,  -5.,   5.);

  //create the output file 
  std:: ofstream file("file.txt");
 
  // loop the events
  int ievt=0;  
  file  <<"{ " <<std::endl;

  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
    // open input file (can be located on castor)
    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
    if( inFile ){
      // ----------------------------------------------------------------------
      // Second Part: 
      //
      //  * loop the events in the input file 
      //  * receive the collections of interest via fwlite::Handle
      //  * fill the histograms
      //  * after the loop close the input file
      // ----------------------------------------------------------------------      
      fwlite::Event ev(inFile);

      for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
	edm::EventBase const & event = ev;
	

	// break loop if maximal number of events is reached 
//	if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
	// simple event counter
//	if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false) 
	std::cout << "  processing event: " << ievt << std::endl;

	file  <<"<| " <<std::endl;

	// Handle to the object collection
	edm::Handle<std::vector<reco::Muon> > muons;
	event.getByLabel(std::string("muons"), muons);

        edm::Handle<std::vector<reco::GsfElectron> > electrons;
        event.getByLabel(std::string("gsfElectrons"), electrons);
    
	edm::Handle<std::vector<reco::PFJet> > jets;
        event.getByLabel(std::string("kt4PFJets"), jets);

        edm::Handle<std::vector<reco::PFMET> > mets;
        event.getByLabel(std::string("pfMet"), mets);

        edm::Handle<std::vector<reco::Photon> > photons;
        event.getByLabel(std::string("photons"), photons);

        edm::Handle<std::vector<reco::Vertex> > vertices;
        event.getByLabel(std::string("offlinePrimaryVertices"), vertices);

	edm::Handle<std::vector<reco::PFTau> > taus;
        event.getByLabel(std::string("hpsPFTauProducer"), taus);


	file  <<" \"Particles\" -> " <<std::endl;  
	file  <<"  <|" <<std::endl;
	file  <<"  \"Muons\"  -> {" <<std::endl;

	// loop object collection and fill histograms
	for(std::vector<reco::Muon>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1){
    	  muonEt_ ->Fill( mu1->et () );
          muonPy_ ->Fill( mu1->py () );
	  muonPx_ ->Fill( mu1->px () );
	  muonPt_ ->Fill( mu1->pt () );
	  muonEta_->Fill( mu1->eta() );
	  muonPhi_->Fill( mu1->phi() );	  

	  //write information about Muons
	  file  <<" <| \"TransverseMomentum\" -> " << mu1->pt() << ", \"Energy\"  -> " << mu1->energy() << ", \"Pseudorapidity\" -> "<< mu1->eta()<<", \"AzimuthalAngle\" -> " << mu1->phi()  << ", \"Charge\" -> " << mu1->charge()  << ", \"Mass\" -> "<< mu1->mass() << " |>, "<< std::endl;
	   
	 
	  if( mu1->pt()>20 && fabs(mu1->eta())<2.1 ){
	    for(std::vector<reco::Muon>::const_iterator mu2=muons->begin(); mu2!=muons->end(); ++mu2){

	      if(mu2>mu1){ // prevent double conting
		if( mu1->charge()*mu2->charge()<0 ){ // check only muon pairs of unequal charge 
		  if( mu2->pt()>20 && fabs(mu2->eta())<2.1 ){
		    mumuMass_->Fill( (mu1->p4()+mu2->p4()).mass() );
		  }
		}
	      }
	    }
	  }
	}

	file  <<" }, " <<std::endl;
        file  <<" \"Electrons\"  -> {" <<std::endl;

	for(std::vector<reco::GsfElectron>::const_iterator el1=electrons->begin(); el1!=electrons->end(); ++el1){
	 
	  //write information about Electrons
          file  <<" <| \"TransverseMomentum\" -> " << el1->pt() << ", \"Energy\"  -> " << el1->energy() << ", \"Pseudorapidity\" -> "<< el1->eta()<<", \"AzimuthalAngle\" -> " << el1->phi()  << ", \"Charge\" -> " << el1->charge()  << ", \"Mass\" -> "<< el1->mass() << " |>, "<< std::endl;
          electronPt_ ->Fill( el1->pt () );
          electronEta_->Fill( el1->eta() );
          electronPhi_->Fill( el1->phi() );
	}
       
	file  <<" }, " <<std::endl;
	file  <<" \"Jets\"  -> {" <<std::endl;

        for(std::vector<reco::PFJet>::const_iterator jet1=jets->begin(); jet1!=jets->end(); ++jet1){
	  //write information about Jets
          file  <<" <| \"TransverseMomentum\" -> " << jet1->pt() << ", \"Energy\"  -> " << jet1->energy() << ", \"Pseudorapidity\" -> "<< jet1->eta()<<",\"AzimuthalAngle\" -> " << jet1->phi()  << ", \"Charge\" -> Missing[\"NotAvailable\"], \"Mass\" -> "<< jet1->mass() << " |>, "<< std::endl;
	  jetPt_ ->Fill( jet1->pt () );
          jetEta_->Fill( jet1->eta() );
          jetPhi_->Fill( jet1->phi() );
        }

	file  <<" }, " <<std::endl;
        file  <<" \"MissingTransverseEnergy\"  -> {" <<std::endl;

	for(std::vector<reco::PFMET>::const_iterator met1=mets->begin(); met1!=mets->end(); ++met1){

	  //write information about Missing Transverse Energy 
          file  <<" <| \"TransverseMomentum\" -> " << met1->pt() << ", \"Energy\"  -> " << met1->energy() << ", \"Pseudorapidity\" ->  Missing[\"NotApplicable\"], \"AzimuthalAngle\" -> " << met1->phi()  << ", \"Charge\" -> " << met1->charge()  << ", \"Mass\" -> "<< met1->mass() << " |>, "<< std::endl;
          metPt_ ->Fill( met1->pt () );
          metEta_->Fill( met1->eta() );
          metPhi_->Fill( met1->phi() );
        }

        file  <<" }, " <<std::endl;
        file  <<" \"Photons\"  -> {" <<std::endl;

        for(std::vector<reco::Photon>::const_iterator ph1=photons->begin(); ph1!=photons->end(); ++ph1){

	  //write information about Photon
          file  <<" <| \"TransverseMomentum\" -> " << ph1->pt() << ", \"Energy\"  -> " << ph1->energy() << ", \"Pseudorapidity\" -> "<< ph1->eta()<<",\
 \"AzimuthalAngle\" -> " << ph1->phi()  << ", \"Charge\" -> " << ph1->charge()  << ", \"Mass\" -> "<< ph1->mass() << " |>, "<< std::endl;
	  photonPt_ ->Fill( ph1->pt () );
          photonEta_->Fill( ph1->eta() );
          photonPhi_->Fill( ph1->phi() );
        }

	file  <<" }, " <<std::endl;
	file  <<"|> " <<std::endl;


        file  <<" , \"ReconstructedObjects\" ->  " <<std::endl;


        file  <<" <|" <<std::endl;

        file  <<"  \"Taus\" -> { " <<std::endl;

        for(std::vector<reco::PFTau>::const_iterator tau1=taus->begin(); tau1!=taus->end(); ++tau1){
      
	  //Write information about Taus
	  file  <<"<| \"TransverseMomentum\" -> " << tau1->pt() << ", \"Energy\"  -> " << tau1->energy() << ", \"Pseudorapidity\" -> "<< tau1->eta()<< " , \"AzimuthalAngle\" -> " << tau1->phi()  << ", \"Charge\" -> " <<tau1->charge()  << ", \"Mass\" -> "<< tau1->mass() << " |>, "<< std::endl;        
	}

        file  <<" }, " <<std::endl;
        file  <<"|> " <<std::endl;

        file  <<" , \"Vertices\" ->  " <<std::endl;


	file  <<" <|" <<std::endl;
        file  <<"  \"PrimaryVertex\" -> { " <<std::endl;

	for(std::vector<reco::Vertex>::const_iterator vtx1=vertices->begin(); vtx1!=vertices->end(); ++vtx1){
	  
	  //write information about primary vertices
	  file  <<"<|  \"Position\"  -> { " << vtx1->x() << " , " << vtx1->y() << " , " << vtx1->z() <<" } |>," <<std::endl;
        }
        file  <<" }, " <<std::endl;

	/*** get beamspot and TransientTrackBuilder from the event/eventSetup ***/
	edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
        event.getByLabel(std::string("offlineBeamSpot"), recoBeamSpotHandle);
	reco::BeamSpot vertexBeamSpot;
        if ( recoBeamSpotHandle.isValid() )
          {
            vertexBeamSpot= *recoBeamSpotHandle;

          } else
          {
	    edm::LogInfo("MyAnalyzer")
              << "No beam spot available from EventSetup \n";
          }
        file  <<" \"VertexBeamSpot\"  -> {" <<std::endl;
	//write information about Vertex Beam Spot
	file  <<"<|  \"Position\" -> {" << vertexBeamSpot.x0() << ", "<< vertexBeamSpot.y0()<< ", "<< vertexBeamSpot.z0()<<"} |>, } "<<std::endl;
                                                                                                                                                  
        file  <<"|> " <<std::endl;


	file  <<" |>, " <<std::endl;
 	
      }  
      file << "}" <<std::endl;
      // close input file
      inFile->Close();
    }
    // break loop if maximal number of events is reached:
    // this has to be done twice to stop the file loop as well
    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
  }
  return 0;
}
