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

#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/VectorUtil.h"
#include <sstream>
#include <typeinfo>
#include <math.h>

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

  TFile *myRootFile;
  TH1D *h_Pt_Sum;
  TH1D *h_Pt_Sum_Bx;
  TH1D *h_Pt_Sum_Bx_Vrtx;

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
}


MuonAnalyzer::~MuonAnalyzer() 
{ 
  std::cout <<" DTOR" << std::endl;
}

void MuonAnalyzer::beginJob()
{  
  myRootFile = new TFile("data.root", "RECREATE");
  h_Pt_Sum = new TH1D("h_Pt_Sum", "h_Pt_Sum; Pt_Sum; #muons", 50, 0. , 50.);
  h_Pt_Sum_Bx = new TH1D("h_Pt_Sum_Bx", "h_Pt_Sum_Bx; Pt_Sum; #muons", 30, 0., 30.); 
  h_Pt_Sum_Bx_Vrtx = new TH1D("h_Pt_Sum_Bx_Vrtx", "h_Pt_Sum_Bx_Vrtx; Pt_Sum; #muons", 30, 0., 30.); 
  std::cout << "HERE MuonAnalyzer::beginJob()" << std::endl;
}

void MuonAnalyzer::endJob()
{ 
  h_Pt_Sum->Write();
  h_Pt_Sum_Bx->Write();
  h_Pt_Sum_Bx_Vrtx->Write();
  myRootFile->Write();
  myRootFile->Close();
  
  delete h_Pt_Sum;
  delete h_Pt_Sum_Bx;
  delete h_Pt_Sum_Bx_Vrtx;
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


	int m = 0;
	int m5 = 0;
	int mI = 0;


	for(std::vector<TrackingParticle>::const_iterator itP = tPC->begin(); itP != tPC->end(); ++itP){
			bx = itP->eventId().bunchCrossing();
			event = itP->eventId().event();

		//if muon with pt bigger than 2
		if(abs(itP->pdgId()) == 13 && itP->pt() >2){

			//++m;
			//std::cout << "TP momentum, q, ID, & Event #: " << itP->p4() << " " << itP->charge() << " " << itP->pdgId() << " bx: "     << itP->eventId().bunchCrossing() << " " << itP->eventId().event() << std::endl;

			//std::cout << " TP Vertex " << itP->vertex() << std::endl;
			//std::cout << " Source vertex " << itP->parentVertex()->position() << std::endl ;

			
			//use for sum of pt of particles around muon, with some conditions
			double PtSum = 0, PtSumBx=0, PtSumBxVrtx=0;

			//loop to check relations between our chosen muon and rest of the particles
			for(std::vector<TrackingParticle>::const_iterator jtP = tPC->begin(); jtP != tPC->end(); ++jtP){
				//if not muon with pt bigger that 1 GeV
				if((abs(jtP->pdgId()) != 13) && jtP->pt() > 1){

					//dr calculation between muon and other particle
					dr = ROOT::Math::VectorUtil::DeltaR(itP->p4() , jtP->p4());

					if(dr < 0.2){
						//Ptsum within dr = 0.2
						PtSum += jtP->pt();
	
						if(itP->eventId().bunchCrossing() == jtP->eventId().bunchCrossing()){
							//same bunchcorssing of muon and j particle
							PtSumBx += jtP->pt();
							if(abs(itP->vz() - jtP->vz()) < 0.2){
								//close vertexes, only in z axis
								PtSumBxVrtx += jtP->pt();
							}
						}
					}						
				}
			}
			std::cout << "Pt sum " << PtSum << std::endl;
			h_Pt_Sum->Fill(PtSum/(itP->pt()));
			h_Pt_Sum_Bx->Fill(PtSumBx/(itP->pt()));
			if(PtSumBxVrtx/(itP->pt()) > 1) h_Pt_Sum_Bx_Vrtx->Fill(PtSumBxVrtx/(itP->pt()));
		}
	}
	/*for(std::vector<TrackingParticle>::const_iterator itP=tPC.begin(); itP != tPC.end(); itP++){
		//std::cout << n++ << std::endl;
		//std::cout << itP->pdgId() << " " << itP->mass() << std::endl;
		++m;
		if(abs(itP->pdgId()) == 13 && itP->pt() > 2){
			
			const std::vector<SimTrack> &simTracks = itP->g4Tracks();
			int size = simTracks.size();
			for (std::vector<SimTrack>::const_iterator it=simTracks.begin(); it<simTracks.end(); it++) {
    			const SimTrack & track= *it;
				if(track.type() == itP->pdgId()){
					++m5;
					std::cout << "Correctly clasified muon" << std::endl;
					std::cout << "Energy sum in cone **** calculatations begin" << std::endl;

					int nC = 0;
					int P = 0;
					double phi = itP->phi();
					double eta = itP->eta();
					double sumPt = 0;

					for(std::vector<TrackingParticle>::const_iterator jtP=tPC.begin(); jtP != tPC.end(); jtP++){
						if(abs(jtP->pdgId()) != 13){
					
							double deltaPhi = (jtP->phi()- phi);
							double deltaEta = (jtP->eta()- eta);
							
							double DeltaR = sqrt(pow(deltaPhi,2)+pow(deltaEta, 2));	
							++P;
							if(DeltaR < 0.1){
								sumPt += jtP->pt();
								++nC;
							}
									
						}

					}
					std::cout << "Energy sum in cone *** calculations end" << std::endl;
					std::cout << "Pt sum in cone " <<sumPt << " of " << nC << " particles from" << P <<  std::endl;
					if(sumPt < 5) ++mI;

				}
				else{
					std::cout<< "Mismatch b/t TrackingParticle and Geant types" << std::endl;
				}
					
			}
		}
		


	}
	std::cout << m << " " << m5 << " " << mI << std::endl; */
/*
  std::cout <<" SIMULATED MUONS: "<<std::endl;
  edm::Handle<edm::SimTrackContainer> simTk;
  ev.getByToken(inputSim, simTk);
  std::vector<SimTrack> mySimTracks = *(simTk.product());
*/  

/*
  for (std::vector<SimTrack>::const_iterator it=mySimTracks.begin(); it<mySimTracks.end(); it++) {
    const SimTrack & track = *it;
    if ( track.type() == -99) continue;
    if ( track.vertIndex() != 0) continue;

    //sucessful muon, add to count
    simMuonCount++;

    phi_sim = track.momentum().phi(); //momentum azimutal angle
    pt_sim = track.momentum().pt(); //transverse momentum
    eta_sim = track.momentum().eta(); //pseudorapidity
    std::cout <<" trackId: " <<track.trackId() 
          << " pt_sim: " << pt_sim <<" eta_sim: "<<eta_sim<<" phi_sim: "<<phi_sim
          <<" vtx: "<<track.vertIndex()<<" type: "<<track.type() // 13 or -13 is a muon
          << std::endl; 
  }
  if(simMuonCount!=1) {
     cout<<"    Simulated muon count != 1"<<endl;
     return;
  }
  */

  // tracking particles




/*
  cout << "Tracking particle" << endl;
  edm::Handle< std::vector< TrackingParticle > > trackingParticleHandle;
  ev.getByToken(TrackingParticleToken_, trackingParticleHandle);
  
  cout << "Number of tracked particles " << trackingParticleHandle->size() << endl;
  int nM = 0;  
  int nM2 = 0;
  int nM5 = 0;
  int nM10 = 0;
  for(unsigned int iTP = 0; iTP < trackingParticleHandle->size(); ++iTP){
	edm::Ptr<TrackingParticle> tpPtr(trackingParticleHandle, iTP);
 	
	if(abs(tpPtr->pdgId()) == 13 ){
		histo->Fill(tpPtr->pt());
		++nM; 
		++tnM;
		if(tpPtr->pt() > 2){
			 ++nM2;
			 ++tnM2;
			}
		if(tpPtr->pt() > 5) ++nM5;
		if(tpPtr->pt() > 10){
			 ++nM10;
		   	 ++tnM10;
			 cout << printTrackigParticleShort(tpPtr) << endl;
			}
	} 		
		
  }*/

/*
  edm::Handle<TrackingParticleCollection> mergedPH;
  edm::Handle<TrackingVertexCollection> mergedVH;

  edm::InputTag trackingTruth = theConfig.getUntrackedParameter<edm::InputTag>("trackingTruth");


  ev.getByLabel(trackingTruth, mergedPH);
  ev.getByLabel(trackingTruth, mergedVH);

  std::cout << "fuck" << std::endl;

  if(theConfig.getUntrackedParameter<bool>("dumpVertexes")) {

 
    std::cout << std::endl << "Dumping merged vertices: " << std::endl;
    for (TrackingVertexCollection::const_iterator iVertex = mergedVH->begin(); iVertex != mergedVH->end(); ++iVertex){
	std::cout << std::endl << &iVertex;
	std::cout << "Daughters of this vertex:" << std::endl;
	//std::cout << iVertex<< std::endl;
    		//for (TrackingParticleRefVector::iterator iTrack = iVertex->daughterTracks_begin(); iTrack != iVertex->daughterTracks_end(); ++iTrack) std::cout << "yes " << std::endl;//std::cout << **iTrack;
    }

    std::cout << "in" <<std::endl;
  }
 
  if (theConfig.getUntrackedParameter<bool>("dumpVertexes")) {
    std::cout << std::endl << "Dumping merged vertices: " << std::endl;
    for (TrackingVertexCollection::const_iterator iVertex = mergedVH->begin(); iVertex != mergedVH->end(); ++iVertex) {
      std::cout << std::endl << *iVertex;
      std::cout << "Daughters of this vertex:" << std::endl;
      for (tp_iterator iTrack = iVertex->daughterTracks_begin(); iTrack != iVertex->daughterTracks_end(); ++iTrack)
        std::cout << **iTrack;
    }
    std::cout << std::endl;
  }

  if (theConfig.getUntrackedParameter<bool>("dumpOnlyBremsstrahlung")) {
    std::cout << std::endl << "Dumping only merged tracks: " << std::endl;
    for (TrackingParticleCollection::const_iterator iTrack = mergedPH->begin(); iTrack != mergedPH->end(); ++iTrack)
      if (iTrack->g4Tracks().size() > 1)
        std::cout << &iTrack << std::endl;
  } else {
    std::cout << std::endl << "Dump of merged tracks: " << std::endl;
    for (TrackingParticleCollection::const_iterator iTrack = mergedPH->begin(); iTrack != mergedPH->end(); ++iTrack)
      std::cout << &iTrack << std::endl;
  }*/

}




DEFINE_FWK_MODULE(MuonAnalyzer);

