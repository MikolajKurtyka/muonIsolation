#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Utilities/interface/InputTag.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/VectorUtil.h"
#include <sstream>
#include <typeinfo>
#include <math.h>
#include <string>
#include <cstring>

//useful for tp
typedef edm::RefVector<std::vector<TrackingParticle>> TrackingParticleContainer;
typedef std::vector<TrackingParticle> TrackingParticleCollection;

//useful fro iterators
typedef TrackingParticleRefVector::iterator tp_iterator;
typedef TrackingVertex::genv_iterator genv_iterator;
typedef TrackingVertex::g4v_iterator g4v_iterator;

/*
bool decaysToMuon(edm::Ptr<TrackingParticle>& tpPtr){
	for(auto& decayVert : tpPtr->decayVertices()){
		for(auto& daughterTrack : decayVert->daughterTracks()){
			if(abs(daughterTrack->pdgId()) == 13) return true;
		}	
	}
}; */


//object definition
class MuonAnalyzer : public edm::EDAnalyzer {
public:

  //constructor, function is called when new object is created
  explicit MuonAnalyzer(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~MuonAnalyzer();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;
  std::string rootFileName;

  TFile *myRootFile;
  TH2D *h_dr_Pt;
  TH2D *h_dr_Pt_Sum;
  TH2D *h_dz_drho;
  TH1D *h_dz;

  //TrackingParticleCollection TrackingParticleCollection_;
  edm::EDGetTokenT< TrackingParticleCollection> vec_TrackingParticle_Token_;


};


MuonAnalyzer::MuonAnalyzer(const edm::ParameterSet& conf) 
  : theConfig(conf), theEventCount(0) 
{
  std::cout <<" CTORXX" << std::endl;
  //TrackingParticleInputTag = theConfig.getParameter<edm::InputTag>("TrackingParticleInputTag");
  //TrackingParticleToken_ = consumes< TrackingParticleCollection>(edm::InputTag("mix", "MergedTrackTruth"));
//edm::InputTag trackingTruth = theConfig.getUntrackedParameter<edm::InputTag>("trackingTruth");
//vec_TrackingParticle_Token_(consumes<TrackingParticleCollection>(trackingTruth));
//vec_TrackingParticle_Token_(consumes<TrackingParticleCollection>( theConfig.getUntrackedParameter<edm::InputTag>("trackingTruth")));
//edm::InputTag trackingTruth = theConfig.getUntrackedParameter<edm::InputTag>("trackingTruth");
vec_TrackingParticle_Token_  = consumes< TrackingParticleCollection>(edm::InputTag("mix", "MergedTrackTruth"));

rootFileName = theConfig.getUntrackedParameter<std::string>("rootFileName");
}


MuonAnalyzer::~MuonAnalyzer() 
{ 
  std::cout <<" DTOR" << std::endl;
}

void MuonAnalyzer::beginJob()
{  
  char crootFileName[rootFileName.length() + 1];
  strcpy(crootFileName, rootFileName.c_str());

  myRootFile = new TFile(crootFileName, "RECREATE");
  h_dz =  new TH1D("h_dz", "h_dz; dz; #n", 1000, -10., 10.);
  h_dr_Pt  = new TH2D("h_dr_Pt", "h_dr_Pt; dr; Pt", 100, 0., 1.,200, 0., 4.);
  h_dr_Pt_Sum  = new TH2D("h_dr_Pt_Sum", "h_dr_Pt_Sum; dr; Pt_Sum", 100, 0., 1., 1000, 0., 5.);
  h_dz_drho = new TH2D("h_dz_drho", "h_dz_drho; dr; drho", 200, -10., 10., 100, 0., 1.);
  std::cout << "HERE MuonAnalyzer::beginJob()" << std::endl;
}

void MuonAnalyzer::endJob()
{ 
  h_dz->Write();
  h_dr_Pt->Write();
  h_dr_Pt_Sum->Write();
  h_dz_drho->Write();
  myRootFile->Write();
  myRootFile->Close();
  
  delete h_dz;
  delete h_dr_Pt;
  delete h_dr_Pt_Sum;
  delete h_dz_drho;
  delete myRootFile;
  std::cout << "HERE MuonAnalyzer::endJob()" << std::endl;
}

void MuonAnalyzer::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
    std::cout << " HERE MuonAnalyzer::analyze "<< std::endl;
    std::cout << "Tracking Particles" << std::endl;

	//Tracking partilce 
	edm::Handle<TrackingParticleCollection>  TruthTrackContainer ;
	ev.getByToken( vec_TrackingParticle_Token_, TruthTrackContainer );
    //TrackingParticleCollection tPC   = *(TruthTrackContainer.product());
	const TrackingParticleCollection *tPC   = TruthTrackContainer.product();
    //const TrackingParticleCollection *tPC   = TruthTrackContainer.product();


	for(std::vector<TrackingParticle>::const_iterator itP = tPC->begin(); itP != tPC->end(); ++itP){
			int bx = itP->eventId().bunchCrossing();

		//if muon with pt bigger than 2
		if(abs(itP->pdgId()) == 13 && itP->pt() >20. && bx == 0){

			//++m;
			

			int nxbins = h_dr_Pt_Sum->GetNbinsX();
			double Pt_Sum[nxbins] = {0.};
			//loop to check relations between our chosen muon and rest of the particles
			for(std::vector<TrackingParticle>::const_iterator jtP = tPC->begin(); jtP != tPC->end(); ++jtP){
				//if not muon with pt bigger that 1 GeV
				if((abs(jtP->pdgId()) != 13) && jtP->pt() > 1){

					
					if((itP->eventId().bunchCrossing() == jtP->eventId().bunchCrossing())){
// && itP->eventId().bunchCrossing()==0){
						double dr = reco::deltaR(itP->eta(), itP->phi(), jtP->eta(), jtP->phi());
						h_dr_Pt->Fill(dr, (jtP->pt())/(itP->pt()));
					
						double drho = sqrt((itP->vx() - jtP->vx())*(itP->vx() - jtP->vx())+(itP->vy() - jtP->vy())*(itP->vy() - jtP->vy()));
						h_dz_drho->Fill(itP->vz() - jtP->vz(), drho);
						h_dz->Fill(itP->vz() - jtP->vz());

						for(int bin = 0; bin < nxbins; ++bin){
							double bin_dr = h_dr_Pt_Sum->GetXaxis()->GetBinLowEdge(bin) + 2* h_dr_Pt_Sum->GetXaxis()->GetBinWidth(bin);
							//std::cout << bin_dr << std::endl;
							if((dr<bin_dr)&&(abs(itP->vz() - jtP->vz()) < 0.2)){
								Pt_Sum[bin] += jtP->pt();
							}
						}
					}
											
				}
			}

			for(int bin = 0; bin < nxbins; ++bin){
				double bin_dr = h_dr_Pt_Sum->GetXaxis()->GetBinLowEdge(bin) + 2* h_dr_Pt_Sum->GetXaxis()->GetBinWidth(bin);
				//std::cout << bin << " " << bin_dr << std::endl;				
				h_dr_Pt_Sum->Fill(bin_dr+0.005, Pt_Sum[bin]/itP->pt());
			}
		}
	}

}




DEFINE_FWK_MODULE(MuonAnalyzer);

