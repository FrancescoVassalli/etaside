#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"
#include "spline.h"
#include <fstream>
#include <array>
#include <stdexcept>

using namespace Pythia8;
using namespace std;

Double_t E= 2.71828182845904523536;


void makedata(string filename){
	Pythia pythia;
	pythia.readString("Beams:eCM = 5020.");
  	pythia.readString("SoftQCD:nonDiffractive = on");
  	pythia.readString("Random::setSeed = on");
  	pythia.readString("Random::seed =0");

  	TFile* f = new TFile(filename.c_str(),"RECREATE");
  	TTree* t=new TTree("tree100","events");
  	float eta;
  	t->Branch("eta",&eta);

  	std::vector<Particles> myparticles(0);
  	Particle ptemp;
  	for(int iEvent=0; iEvent<1000; iEvent++){
  		if(iEvent%30==0)  
  			cout<<"Event N: "<<iEvent<<'\n';
    	if (!pythia.next()){
      		cout<<"pythia.next() failed"<<"\n";
      		continue;
      	for(unsigned i=0; i<pythia.event.size();i++){
      		if(pythia.event[i].isCharged()&&pythia.event[i].isFinal())
      			ptemp= pythia.event.at(i);
      			eta = pythia.event.at(i).eta();
      			t->Fill();
      			myparticles.push_back(ptemp);
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