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
  TH2D *h_dr_Pt_Sum_Ch;
  TH2D *h_dr_Pt_Sum_dz05;
  TH2D *h_dr_Pt_Sum_dz05_Ch;
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
//consumes trackingcolections edm::input tag gmtmuons eventualnie L1 + buildfile

//rootFileName = theConfig.getUntrackedParameter<std::string>("rootFileName");
}


MuonAnalyzer::~MuonAnalyzer() 
{ 
  std::cout <<" DTOR" << std::endl;
}

void MuonAnalyzer::beginJob()
{  
  //char crootFileName[rootFileName.length() + 1];
  //strcpy(crootFileName, rootFileName.c_str());

  myRootFile = new TFile(theConfig.getParameter<std::string>("rootFileName").c_str(), "RECREATE");
  h_dz =  new TH1D("h_dz", "h_dz; dz; #n", 1000, -2., 2.);
  h_dr_Pt  = new TH2D("h_dr_Pt", "h_dr_Pt; dr; Pt", 100, 0., 1.,200, 0., 4.);
  h_dr_Pt_Sum  = new TH2D("h_dr_Pt_Sum", "h_dr_Pt_Sum; dr; Pt_Sum", 100, 0., 1., 10000, 0., 50.);
  h_dr_Pt_Sum_Ch  = new TH2D("h_dr_Pt_Sum_Ch", "h_dr_Pt_Sum_Ch; dr; Pt_Sum", 100, 0., 1., 10000, 0., 50.);
  h_dr_Pt_Sum_dz05  = new TH2D("h_dr_Pt_Sum_dz05", "h_dr_Pt_Sum; dr; Pt_Sum", 100, 0., 1., 10000, 0., 50.);
  h_dr_Pt_Sum_dz05_Ch  = new TH2D("h_dr_Pt_Sum_dz05_Ch", "h_dr_Pt_Sum; dr; Pt_Sum", 100, 0., 1., 10000, 0., 50.);
  h_dz_drho = new TH2D("h_dz_drho", "h_dz_drho; dz; drho", 200, -10., 10., 100, 0., 1.);
  std::cout << "HERE MuonAnalyzer::beginJob()" << std::endl;
}

void MuonAnalyzer::endJob()
{ 
  std::cout << "endjob" << std::endl;
  h_dz->Write();
  h_dr_Pt->Write();
  h_dr_Pt_Sum->Write();
  h_dr_Pt_Sum_Ch->Write();
  h_dr_Pt_Sum_dz05->Write();
  h_dr_Pt_Sum_dz05_Ch->Write();
  h_dz_drho->Write();
  myRootFile->Write();
  myRootFile->Close();
  
  delete h_dz;
  delete h_dr_Pt;
  delete h_dr_Pt_Sum;
  delete h_dr_Pt_Sum_Ch;
  delete h_dr_Pt_Sum_dz05;
  delete h_dr_Pt_Sum_dz05_Ch;
  delete h_dz_drho;
  delete myRootFile;
  std::cout << "HERE MuonAnalyzer::endJob()" << std::endl;
}

void MuonAnalyzer::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
    std::cout << " HERE MuonAnalyzer::analyze "<< std::endl;
	std::cout << ev.id() << std::endl;
	//Tracking partilce 
	edm::Handle<TrackingParticleCollection>  TruthTrackContainer ;
	ev.getByToken( vec_TrackingParticle_Token_, TruthTrackContainer );
    //TrackingParticleCollection tPC   = *(TruthTrackContainer.product());
	const TrackingParticleCollection *tPC   = TruthTrackContainer.product();
    //const TrackingParticleCollection *tPC   = TruthTrackContainer.product();
	
	bool debug = 1;
	int nev = 0;	

	for(std::vector<TrackingParticle>::const_iterator itP = tPC->begin(); itP != tPC->end(); ++itP){
			int bx = itP->eventId().bunchCrossing();

		//if muon with pt bigger than 2
		if(abs(itP->pdgId()) == 13 && itP->pt() >3. && bx == 0){

			int nxbins = h_dr_Pt_Sum->GetNbinsX();
			double Pt_Sum[nxbins] = {0.};
            double Pt_Sum_Ch[nxbins] = {0.};
            double Pt_Sum_dz05[nxbins] = {0.};
            double Pt_Sum_dz05_Ch[nxbins] = {0.};
			//loop to check relations between our chosen muon and rest of the particles
			for(std::vector<TrackingParticle>::const_iterator jtP = tPC->begin(); jtP != tPC->end(); ++jtP){
				//if not muon with pt bigger that 1 GeV
				if((abs(jtP->pdgId()) != 13) && jtP->pt() > 3.){
					
					if((itP->eventId().bunchCrossing() == jtP->eventId().bunchCrossing())){
// && itP->eventId().bunchCrossing()==0){
						double dr = reco::deltaR(itP->eta(), itP->phi(), jtP->eta(), jtP->phi());
						h_dr_Pt->Fill(dr, (jtP->pt())/(itP->pt()));
					
						double drho = sqrt((itP->vx() - jtP->vx())*(itP->vx() - jtP->vx())+(itP->vy() - jtP->vy())*(itP->vy() - jtP->vy()));
						h_dz_drho->Fill(itP->vz() - jtP->vz(), drho);
						h_dz->Fill(itP->vz() - jtP->vz());

						double t_ip = sqrt(pow((jtP->vx()+jtP->py()*(jtP->charge()/0.003/3.8)), 2) +  pow((jtP->vy()-jtP->px()*(jtP->charge()/0.003/3.8)), 2) ) - jtP->pt()/(0.003*3.8);

						if(debug){
							std::cout << "------------TrackingParticle------------" << std::endl;
							std::cout << "Nr, eventid: " << nev << " " << itP->eventId().event() << std::endl;
							std::cout << "Muon (pt,eta, phi, z0, x0, y0) " << std::endl;
							std::cout << itP->pt()/0.025 << " " << itP->eta() << " " << itP->phi() << " " << itP->vz() << " " << itP->vx() << " " << itP->vy() << std::endl;							
							std::cout << "Particle (pt, eta, phi, z0, x0, y0) " << std::endl;
							std::cout << jtP->pt()/0.025 << " " << jtP->eta() << " " << jtP->phi() << " " << jtP->vz() << " " << jtP->vx() << " " << jtP->vy() << std::endl;
							std::cout << "t_ip" << std::endl;
							std::cout << t_ip << std::endl;
							std::cout << "dr , dz " << " " << dr << " " <<  itP->vz() - jtP->vz() << std::endl;
							if((dr<0.3)&&(abs(itP->vz() - jtP->vz()) < 0.2)) std::cout << "TRACK PASS" << std::endl;
							else std::cout << "DIDN'T PASS" << std::endl;
						}
					
						for(int bin = 0; bin < nxbins; ++bin){
							double bin_dr = h_dr_Pt_Sum->GetXaxis()->GetBinLowEdge(bin) + 2* h_dr_Pt_Sum->GetXaxis()->GetBinWidth(bin);
							//std::cout << bin_dr << std::endl;
							if((dr<bin_dr)&&(abs(itP->vz() - jtP->vz()) < 0.2)){
								Pt_Sum[bin] += jtP->pt();
							}
							if((dr<bin_dr)&&(abs(itP->vz() - jtP->vz()) < 0.05)){
								Pt_Sum_dz05[bin] += jtP->pt();
							}
							if((dr<bin_dr)&&(abs(itP->vz() - jtP->vz()) < 0.2)&&(jtP->charge() != 0)){

								Pt_Sum_Ch[bin] += jtP->pt();
							}
							if((dr<bin_dr)&&(abs(itP->vz() - jtP->vz()) < 0.05)&&(jtP->charge() != 0)){
								Pt_Sum_dz05_Ch[bin] += jtP->pt();
							}
						}
					}
											
				}
			}

			for(int bin = 0; bin < nxbins; ++bin){
				double bin_dr = h_dr_Pt_Sum->GetXaxis()->GetBinLowEdge(bin) + 2* h_dr_Pt_Sum->GetXaxis()->GetBinWidth(bin);
				//std::cout << bin << " " << bin_dr << std::endl;				
				h_dr_Pt_Sum->Fill(bin_dr+h_dr_Pt_Sum->GetXaxis()->GetBinWidth(bin)/2., Pt_Sum[bin]/itP->pt());
				h_dr_Pt_Sum_dz05->Fill(bin_dr+h_dr_Pt_Sum->GetXaxis()->GetBinWidth(bin)/2., Pt_Sum_dz05[bin]/itP->pt());
				h_dr_Pt_Sum_Ch->Fill(bin_dr+h_dr_Pt_Sum->GetXaxis()->GetBinWidth(bin)/2., Pt_Sum_Ch[bin]/itP->pt());
				h_dr_Pt_Sum_dz05_Ch->Fill(bin_dr+h_dr_Pt_Sum->GetXaxis()->GetBinWidth(bin)/2., Pt_Sum_dz05_Ch[bin]/itP->pt());
			}
		    std::cout << "Isolation (PSUM, RPSUM) dr=0.3 : "  << Pt_Sum_Ch[30] << " : " << Pt_Sum_Ch[30]/itP->pt() << std::endl;
		}
		nev++;
	}

}




DEFINE_FWK_MODULE(MuonAnalyzer);

