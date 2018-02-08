#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include "Pythia8/Basics.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace Pythia8;
using namespace std;

Double_t E= 2.71828182845904523536;

TLorentzVector* pToTLV(Vec4 in){
	Double_t px = (double)in.px();
	Double_t py = (double)in.py();
	Double_t pz = (double)in.pz();
	Double_t e =  (double)in.e();
	TLorentzVector *out = new TLorentzVector(px,py,pz,e);
	return out;
}

TVector3 v3fromv4(Vec4 in){
	TVector3 r = TVector3(in.px(),in.py(),in.pz());
}

void check(std::vector<Particle> v, std::vector<float> etas2){
	cout<<"back: "<<v.back().y()<<"\n";
	cout<<"eta diff"<<v.back().y()-v[0].y()<<" : "<<etas2.back()-etas2[0]<<"\n"; // see if the delta eta is boost invarient
}

void makedata(string filename){
	Pythia pythia;
	pythia.readString("Beams:eCM = 5020.");
  	pythia.readString("SoftQCD:nonDiffractive = on");
  	pythia.readString("Random::setSeed = on");
  	pythia.readString("Random::seed =0");
  	pythia.init();

  	TFile* f = new TFile(filename.c_str(),"RECREATE");
  	TTree* t=new TTree("tree100","events");
  	float eta, eta2;
  	t->Branch("eta",&eta);
  	t->Branch("etaboost",&eta2);
  	const float boostB = .434;
  	TVector3 myt3;
  	std::vector<Particle> myparticles(0);
  	Particle ptemp;
  	std::vector<float> eta2s(0);
  	for(int iEvent=0; iEvent<1000; iEvent++){
  		if(iEvent%30==0)  
  			cout<<"Event N: "<<iEvent<<'\n';
    	if (!pythia.next()){
      		cout<<"pythia.next() failed"<<"\n";
      		continue;
      	}
      	for(unsigned i=0; i<pythia.event.size();i++){
      		if(pythia.event[i].isCharged()&&pythia.event[i].status()>80){
      			TLorentzVector tlv;
      			ptemp= pythia.event.at(i);
      			//tlv= pToTLV(ptemp.p());
      			eta = pythia.event.at(i).eta();
      			tlv.Boost(0,0,boostB);
      			eta2 = (float) tlv->Eta();
      			t->Fill();
      			if(iEvent%30==0){
      				cout<<"preboost: "<<eta<<" post: "<<eta2<<"\n";
      			}
      			myparticles.push_back(ptemp);
      			eta2s.push_back(eta2);
      		}
      	}
  	}
  	//check(myparticles,eta2s);
  	t->Write();
  	f->Write();
  	f->Close();
  	delete f;
  	f=NULL;
}


int main(int argc, char *argv[]){
	std::string filename;
	if(argc!=2){
		std::cout<<"accepts 1 arguments: 1. outfile "<<'\n';
		return 1;
	}
	else{
		filename=argv[1];
		makedata(filename);
	}
	return 0;
}