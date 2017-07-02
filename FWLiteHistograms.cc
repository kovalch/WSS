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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"


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

  std:: ofstream file("file.txt");

  // loop the events
  int ievt=0;  
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


file  <<"EventNumber : " << ievt <<std::endl;	
	// loop object collection and fill histograms
	for(std::vector<reco::Muon>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1){
    	  muonEt_ ->Fill( mu1->et () );
          muonPy_ ->Fill( mu1->py () );
	  muonPx_ ->Fill( mu1->px () );
	  muonPt_ ->Fill( mu1->pt () );
	  muonEta_->Fill( mu1->eta() );
	  muonPhi_->Fill( mu1->phi() );	  

	  file  <<"Object : Muon, TransverseMomentum : " << mu1->pt() << ", Momentum  : { " << mu1->px() << " , " << mu1->py() << " , " << mu1->pz() << " } , TransverseEnergy :  " << mu1->et() << ", Pseudorapidity : "<< mu1->eta()<<", AzimuthalAngle : " << mu1->phi()  << ", PolarAngle : "<<mu1->theta() <<std::endl;
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
       
        for(std::vector<reco::GsfElectron>::const_iterator el1=electrons->begin(); el1!=electrons->end(); ++el1){
	  file  <<"Object : Electron, TransverseMomentum : " << el1->pt() << ", Momentum  : { " << el1->px() << " , " << el1->py() << " , " << el1->pz() << " } , TransverseEnergy :  " << el1->et() << ", Pseudorapidity : "<< el1->eta()<<", AzimuthalAngle : " << el1->phi()  << ", PolarAngle : "<<el1->theta() <<std::endl;
          electronPt_ ->Fill( el1->pt () );
          electronEta_->Fill( el1->eta() );
          electronPhi_->Fill( el1->phi() );
	}
       
        for(std::vector<reco::PFJet>::const_iterator jet1=jets->begin(); jet1!=jets->end(); ++jet1){
	  file  <<"Object : Jet, TransverseMomentum : " << jet1->pt() << ", Momentum  : { " << jet1->px() << " , " << jet1->py() << " , " << jet1->pz() << " } , TransverseEnergy :  " << jet1->et() << ", Pseudorapidity : "<< jet1->eta()<<", AzimuthalAngle : " << jet1->phi()  << ", PolarAngle : "<<jet1->theta() <<std::endl;
	  jetPt_ ->Fill( jet1->pt () );
          jetEta_->Fill( jet1->eta() );
          jetPhi_->Fill( jet1->phi() );
        }

	for(std::vector<reco::PFMET>::const_iterator met1=mets->begin(); met1!=mets->end(); ++met1){
	  file  <<"Object : MET, TransverseMomentum : " << met1->pt() << ", Momentum  : { " << met1->px() << " , " << met1->py() << " , " << met1->pz() << " } , TransverseEnergy :  " << met1->et() << ", Pseudorapidity : "<< met1->eta()<<", AzimuthalAngle : " << met1->phi()  << ", PolarAngle : "<<met1->theta() <<std::endl;
          metPt_ ->Fill( met1->pt () );
          metEta_->Fill( met1->eta() );
          metPhi_->Fill( met1->phi() );
        }

        for(std::vector<reco::Photon>::const_iterator ph1=photons->begin(); ph1!=photons->end(); ++ph1){
	  file  <<"Object : Photon, TransverseMomentum : " << ph1->pt() << ", Momentum  : { " << ph1->px() << " , " << ph1->py() << " , " << ph1->pz() << " } , TransverseEnergy :  " << ph1->et() << ", Pseudorapidity : "<< ph1->eta()<<", AzimuthalAngle : " << ph1->phi()  << ", PolarAngle : "<<ph1->theta() <<std::endl;
	  photonPt_ ->Fill( ph1->pt () );
          photonEta_->Fill( ph1->eta() );
          photonPhi_->Fill( ph1->phi() );
        }
	for(std::vector<reco::Vertex>::const_iterator vr1=vertices->begin(); vr1!=vertices->end(); ++vr1){
	  file  <<"Object : Vertex, Position  : { " << vr1->x() << " , " << vr1->y() << " , " << vr1->z() <<" }" <<std::endl;
         
         
         
        }


 	file  <<"---- " <<std::endl;
      }  
      // close input file
      inFile->Close();
    }
    // break loop if maximal number of events is reached:
    // this has to be done twice to stop the file loop as well
    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
  }
  return 0;
}
