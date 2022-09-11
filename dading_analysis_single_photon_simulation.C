#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TH2.h>
#include <TMatrixD.h>
#include <TNtuple.h>
#include <TF1.h>

#include "/phenix/plhf/henpdading/HELIOS_simulation/offline/AnalysisTrain/DileptonAnalysis/MyEvent.h"
#include "/phenix/plhf/henpdading/HELIOS_simulation/offline/AnalysisTrain/DileptonAnalysis/Reconstruction.h"
#include "/phenix/plhf/henpdading/scratch/WriteEvent.h"

using namespace std;
using namespace DileptonAnalysis;

static const double pi = TMath::Pi();
static const double Me = 0.000510998918;  // Particle Data Group 7-25-2004
static const double Me2 = Me*Me;

const double DCENTERCUT = 10; // RICH ghost cut
const double PC1_DPHI_CUT = 0.02, PC1_DZ_CUT = 0.5; // PC1 ghost cut


double getDcenter(double phi1, double z1, double phi2, double z2)
{
	double dcenter_phi_sigma = 0.01, dcenter_phi_offset = 0.;
	double dcenter_z_sigma = 3.6, dcenter_z_offset = 0.;

	double dcenter_z = (z1-z2-dcenter_z_offset)/dcenter_z_sigma;
	double dcenter_phi = (phi1-phi2-dcenter_phi_offset)/dcenter_phi_sigma;

	return sqrt(dcenter_phi*dcenter_phi+dcenter_z*dcenter_z);
}


bool identify_electron_from_DC(MyTrack* mytrk)
{
	double _p = sqrt(pow(mytrk->GetPx(),2) + pow(mytrk->GetPy(),2) + pow(mytrk->GetPz(),2));
	double _theta = mytrk->GetThe0();

	double _pt = _p * sin(_theta);
	if (_pt<0.20) return false;

	double _e = mytrk->GetEcore();
	if(_e/_p < 0.5) return false;

	double _n0 = mytrk->GetN0();
	if(_n0 < 2) return false;

	double _zed = mytrk->GetZDC();
	if(fabs(_zed) > 75) return false;

	int _charge = mytrk->GetCharge();
	if(fabs(_charge) != 1) return false;

	double _emcdz = mytrk->GetEmcdz();
	if(fabs(_emcdz) > 20) return false;

	double _emcdphi = mytrk->GetEmcdphi();
	if(fabs(_emcdphi) > 0.05) return false;

	return true;
}


//reconstruction algortihm
bool is_conv_photon(MyTrack* mytrk1, MyTrack* mytrk2, MyPair* mypair, Reconstruction* reco, float zVtx)
{
	float DPHI=0.005, R_LO=1, R_HI=29, DZED_CUT = 4;//dphi,radius,dzed cuts

	double crkz_A1    = mytrk1->GetCrkz();
	double crkphi_A1  = mytrk1->GetCrkphi();
	double crkz_A2    = mytrk2->GetCrkz();
	double crkphi_A2  = mytrk2->GetCrkphi();

	double dcenter = getDcenter(crkphi_A1, crkz_A1, crkphi_A2, crkz_A2);
	if ( dcenter<DCENTERCUT ) return false;//dcenter cut

	float phi_A1 = mytrk1->GetPhi0();
	float zed_A1 = mytrk1->GetZDC();
	float phi_A2 = mytrk2->GetPhi0();
	float zed_A2 = mytrk2->GetZDC();

	if (fabs(phi_A1-phi_A2)<PC1_DPHI_CUT && fabs(zed_A1-zed_A2)<PC1_DZ_CUT) return false;

        /*
	float dzed = -9999;
        
	if ( mytrk1->GetCharge()==1 )//e-e+
	{//1 is positron, 2 is electron
		dzed = zed_A1-zed_A2;
	}

	if ( mytrk1->GetCharge()==-1 )//e+e-
	{//1 is electron, 2 is positron
		dzed = zed_A2-zed_A1;
	}

	if(fabs(dzed) > DZED_CUT) return false;
	*/
	reco->findIntersection(mytrk1, mytrk2, mypair, zVtx);//find the intersection of both tracks
	float radius_r = mypair->GetRPair();
	if ( (radius_r<=R_LO) || (radius_r>=R_HI) ) return false;//radius cut

	float dphi_r = mypair->GetPhiElectron()-mypair->GetPhiPositron();
	if ( dphi_r>TMath::Pi() )  dphi_r = 2*TMath::Pi()-dphi_r; // 1.3pi->0.7pi
	if ( dphi_r<-TMath::Pi() ) dphi_r = -2*TMath::Pi()-dphi_r; // -1.3pi->-0.7pi
	if ( TMath::Abs(dphi_r)>=DPHI ) return false;//dphi cut

	return true;
}

void dading_analysis_single_photon_simulation(const char* inFile = "rdDST_00000.root", const char* outFile = "out_00000.root", int num = 0)
{
	gSystem->Load("/phenix/plhf/henpdading/install/lib/libDileptonAnalysisEvent");
	gSystem->Load("/phenix/plhf/henpdading/install/lib/libDileptonAnalysisReco");
	gSystem->Load("/phenix/plhf/henpdading/HELIOS_simulation/WriteEvent_C.so");
	Reconstruction reco("/phenix/plhf/roli/install/share/DileptonAnalysis/lookup_3D_one_phi.root");

	float npt[21]={0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0};//flat pT ranges
	TH1F* reconstructed=new TH1F("reconstructed","reconstructed pairs with different pt",20,npt);//histogram for reconstructed pairs
	TH1F* total=new TH1F("total","total pairs with different pt",20,npt);//histogram for pairs that are generating without using the reco algorithm
	TH1F* photon_distribution=new TH1F("photon_distribution","number of photons with different pt",20,npt);//histogram for photon distribution along the pt


	float flat_pt[31]={0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0};
	TH1F* pion_distribution=new TH1F("pion_distribution","number of pions with different pt",30,flat_pt);//histogram for pion distribution along the pt





	TLorentzVector pion(0.,0.,0.,0.);//lorentz vecotor for pi0

	TFile* f1 = new TFile(inFile,"READ");//read input file
	if(!(f1))
	{
		cout<<"can't find input file..."<<endl;
		exit(1);
	}

	TTree* T = (TTree*)f1->Get("T");
	TBranch* br = T->GetBranch("MyEvent");
	MyEvent* event=0;
	br ->SetAddress(&event);


	int nevt = T->GetEntries();

	int nfolder = num; //root file number

	int nstart = nfolder * 100000;

	TFile* helios = new TFile("/direct/phenix+u/roli/scratch/HELIOS/simulation/singlePizero_equalmix_flateta_1B.root","READ");//read the file that contains the pi0 decay information

	TTree* helios_T = (TTree*)helios->Get("T");
	TBranch* helios_br = helios_T->GetBranch("MyEvent");

	WriteEvent* helios_event = 0;
	helios_br->SetAddress(&helios_event);

	for (int ievt = 0; ievt < nevt; ++ievt){


		if(ievt%10000 == 0) cout << "Event:  " << ievt << "/" << nevt << endl;//counting the events

		event->ClearEvent();
		br->GetEntry(ievt);

		int ntrk = event->GetNtrack();  
		int nclust = event->GetNcluster();  
		float zVtx = event->GetVtxZ();

		if(ntrk != 2) {//number of tracks need to be 2
			continue;
		}


		bool dalitz = false;//initialize dalitz

		int id = event->GetEvtNo() - 1;
		helios_br->GetEntry(abs(nstart+id));

		// Get all info from HELIOS output

		for(int i = 0; i < helios_event->GetNEntries(); i++){

			WriteTrack helios_track = helios_event->GetWriteTrack(i);//create helios_track as an instance that has the access to the track info in helio events

			if(helios_track.GetID() == 111){ //eta is 221; pi0 is 111
				pion.SetPxPyPzE(helios_track.GetPx(), helios_track.GetPy(), helios_track.GetPz(), helios_track.GetEnergy());//set px,py,pz,e for pi0
				pion_distribution->Fill(pion.Pt());//fill the number of pions in different pt

				
				if(helios_track.GetBranch() == 1) dalitz = true;//if dalitz decay found
			}
		}
		cout<<"check dalitz: "<<dalitz<<endl;		

		if(dalitz) continue;//if it is dalitz decay, ignore them


		//remove ghost tracks

		int ghost_A[40], real_A[40];
		int all_A[40];
		int nghost_A = 0, nreal_A = 0, nall_A = 0;


		if(event->GetNMcTrack() < 2) continue;

		for (int k1 = 0; k1 < ntrk; ++k1)
		{
			if ( !identify_electron_from_DC( &( event->GetEntry(k1)))) continue; // eid cuts
			all_A[nall_A] = k1;
			++nall_A;        
			int ghost_flag = 0;

			for (int k2 = 0; k2 < ntrk; ++k2)
			{ // make sure k2 start from 0
				if(k1 == k2) continue;

				if ( !identify_electron_from_DC( &( event->GetEntry(k2)))) continue;

				// found a RICH ghost pair
				if ( getDcenter((event->GetEntry(k1)).GetCrkphi(), (event->GetEntry(k1)).GetCrkz(), (event->GetEntry(k2)).GetCrkphi(), (event->GetEntry(k2)).GetCrkz()) < DCENTERCUT ) { ghost_flag = 1; break; }
			}
			if ( ghost_flag ) { ghost_A[nghost_A] = k1; ++nghost_A; }
			else { real_A[nreal_A] = k1; ++nreal_A; }
		}

		if(nreal_A < 2){
			continue;
		}


		TLorentzVector p1, p2;
		TLorentzVector convPhoton;

		bool event_moreThanOne=false;//will use this variable here to filter out the event that has more than one convPhoton as the parents for both tracks
		bool event_Conv=false;//will use this varibale here to filter out the event that has convPhoton as the parent for both track.It is a criteria used for dalitz decay

                   //The necessary criteria for convPhoton in pi0 decays is set up here 
		for(int i=0; i<event->GetNMcTrack(); i++){
			McTrack mctrk_1 = event->GetMcEntry(i);//track 1

			for(int j=i+1; j<event->GetNMcTrack(); j++){
				McTrack mctrk_2 = event->GetMcEntry(j);//track 2

				//if track 1 is e- track and track 2 is e+ track   OR     if track 1 is e+ track and track 2 is e- track
				if((mctrk_1.GetParticleID()==3&&mctrk_2.GetParticleID()==2) || (mctrk_1.GetParticleID()==2 && mctrk_2.GetParticleID()==3)){

			//following criteria are set up for pi0 decays to gamma gamma,making sure two tracks has the same parent is very important for later analysis 
						
				float Parent1_pT=sqrt(mctrk_1.GetParentPx()*mctrk_1.GetParentPx()+mctrk_1.GetParentPy()*mctrk_1.GetParentPy());//pT for track1 parent
				float Parent2_pT=sqrt(mctrk_2.GetParentPx()*mctrk_2.GetParentPx()+mctrk_2.GetParentPy()*mctrk_2.GetParentPy());//pT for track2 parent
					cout<<"Parent1_pT: "<<Parent1_pT<<endl;
					cout<<"Parent2_pT: "<<Parent2_pT<<endl;

					if(Parent1_pT!=Parent2_pT){//if the pt between parent1 and parent2 is not the same

					event_moreThanOne=true;//such event has more than one convPhoton as the parents for both tracks
					cout<<"The event that has more than one conversion: "<<ievt<<endl;	

				        }else{


					//if parent has the same pT, fill this same parent(photon) 
					photon_distribution->Fill(Parent2_pT);//you can also fill with parent2 pt


					}				        	
											
					//for dalitz decay, uncomment the following condition
					
					/*
					if(mctrk_1.GetParentID()==1||mctrk_2.GetParentID()==1){//if one of tracks has convPhoton as the parent 
						event_Conv=true;//such event has convPhoton as the parent for one of the tracks
						cout<<"parent1: "<<mctrk_1.GetParentID()<<endl;
						cout<<"parent2: "<<mctrk_2.GetParentID()<<endl;
						cout<<"The event that has conversion photon as a parent for one of the tracks: "<<ievt<<endl;	
							
					}

			                */
							
							   
							   




			}

		}
	}




                      			//if(event_Conv) continue;//filter out the events
  					  if(event_moreThanOne) continue;//fliter out the events



								//looping with PHcentral track							
	

							for (int ireal_A1 = 0; ireal_A1 < nreal_A; ireal_A1++){ 

								MyTrack mytrk_A1 = event->GetEntry(real_A[ireal_A1]);//track 1

								for (int ireal_A2 = ireal_A1+1; ireal_A2 < nreal_A;ireal_A2++){ 

									MyTrack mytrk_A2 = event->GetEntry(real_A[ireal_A2]);//track 2
									if (mytrk_A1.GetCharge() == mytrk_A2.GetCharge()) continue;//should have different charge
									if (mytrk_A1.GetArm() != mytrk_A2.GetArm()) continue;//should be at the same arm
									
									//setting px,py,pz,e for particle at track 1
									p1.SetX(mytrk_A1.GetPx());
									p1.SetY(mytrk_A1.GetPy());
									p1.SetZ(mytrk_A1.GetPz());
									p1.SetE(sqrt(pow(p1.P(),2) + Me2));


									//setting px,py,pz,e for particle at track 2
									p2.SetX(mytrk_A2.GetPx());
									p2.SetY(mytrk_A2.GetPy());
									p2.SetZ(mytrk_A2.GetPz());
									p2.SetE(sqrt(pow(p2.P(),2) + Me2));
									convPhoton = p1 + p2;//combine both information from particles at track1 and track2, should give us the information of convPhoton

									total->Fill(convPhoton.Pt());//now we can get the number of pairs production(no reco algorithm applied so far) 



									//histogram for all ee pairs reconstructed with the Reco algorithm

									MyPair mypair_AA;
									if (is_conv_photon(&mytrk_A1, &mytrk_A2, &mypair_AA, &reco, zVtx) ){//if reco alrogitm identify that it is a convPhoton


										if(convPhoton.M() < 0.04 || convPhoton.M() > 0.12) continue;


										reconstructed->Fill(convPhoton.Pt());//now we can get the number of pairs production with reco algorithm
									}

								

							
						
					


				
			}
		}



	} //event loop

	TFile* out = new TFile(outFile,"RECREATE");
	out->cd();

	total->Write();
	reconstructed->Write();
	photon_distribution->Write();
	pion_distribution->Write();
	out->Close();

	cout << "All done !" << endl;

}
