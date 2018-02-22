#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include "Pythia8/Basics.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1F.h"

using namespace Pythia8;
using namespace std;

Double_t E= 2.71828182845904523536;

struct StarParticle
{
  float eta;
  TVector3 p;
};

namespace global{
  const int nBins=53;
  double bins[54]= {-2.65,-2.55,-2.45,-2.35,-2.25,-2.15,-2.05,-1.95,-1.85,-1.75,-1.65,-1.55,-1.45,-1.35,-1.25,-1.15,-1.05,-0.95,-0.85,-0.75,-0.65,-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95,2.05,2.15,2.25,2.35,2.45,2.55,2.65};
  TH1F *binfinder = new TH1F("binner","",nBins,bins);
}

TLorentzVector* pToTLV(Vec4 in){
	Double_t px = (double)in.px();
	Double_t py = (double)in.py();
	Double_t pz = (double)in.pz();
	Double_t e =  (double)in.e();
	TLorentzVector *out = new TLorentzVector(px,py,pz,e);
	return out;
}
inline int calcBinNumber(float in){
  int r = global::binfinder->FindBin((double)in,0,0);
  if(!(r>=0&&r<54)){
    r=-1;
  }
  return r;
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
  int binNumber;
  float beta[global::nBins+1];
  string temp;
  string etastr="eta";
  for(int i=0; i<=global::nBins+1;i++){
    global::bins[i] = global::bins[i]-.05; //accidentally put the bin middles in instead of the bin lows
    temp = etastr+ std::to_string(i);
    t->Branch(temp.c_str(),&beta[i]);
  }
  TLorentzVector *tlv=NULL; // add deletes
  TLorentzVector *tl2 = NULL;
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
      if(tlv!=NULL){
      	delete tlv;
      	tlv=NULL;
      }
      if(tl2!=NULL){
        delete tl2;
        tl2=NULL;
      }
      for(int i=0; i<pythia.event.size();i++){
      	if(pythia.event[i].isCharged()&&pythia.event[i].status()>80){
          for (int j = 0; j <= global::nBins; ++j)
          {
            beta[j] = -10;
          }
      		ptemp= pythia.event.at(i);
      		tlv= pToTLV(ptemp.p());
      		eta = pythia.event.at(i).eta();
          cout<<"pre: "<<eta;
      		tlv->Boost(0,0,boostB);
      		eta2 = (float) tlv->Eta();
          cout<<" post: "<<eta2;
          tl2 = new TLorentzVector(tlv->Px(),tlv->Py(),tlv->Pz(),tlv->E());
          binNumber = calcBinNumber(eta2);
          //binNumber--;
          if (binNumber>=0)
          {
            tl2->Boost(0,0,-1*boostB);
            beta[binNumber]= (float) tl2->Eta();
            cout<<" starred: "<<beta[binNumber]<<'\n';
            t->Fill();
          }
        }
      }
  	}
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