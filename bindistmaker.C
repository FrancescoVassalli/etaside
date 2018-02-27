#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include "TChain.h"

using namespace std;

double addError(double e1, double e2){
	return TMath::Power((e1*e1+e2*e2),.05);
}

void bindistmaker(){
	TCanvas *tc = new TCanvas();
	TChain *dijet_tree = new TChain("tree100");
	TFile *output = new TFile("bindists.root","recreate");
	string fname = "etadist.root";
	dijet_tree->Add(fname.c_str());
	gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    const int nBins = 53;
    double bins[nBins+1]= {-2.65,-2.55,-2.45,-2.35,-2.25,-2.15,-2.05,-1.95,-1.85,-1.75,-1.65,-1.55,-1.45,-1.35,-1.25,-1.15,-1.05,-0.95,-0.85,-0.75,-0.65,-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95,2.05,2.15,2.25,2.35,2.45,2.55,2.65};
    for(int i=0; i<=nBins;i++){
    	bins[i] = bins[i]-.05; //accidentally put the bin middles in instead of the bin lows
    }
    const int myBinN = 80;
    double mybins[myBinN+1];
    const float binstart = -4;
    for (int i = 0; i <= myBinN; ++i)
    {
    	mybins[i] = binstart+.1*i;
    }

    std::vector<TH1F*> splithists(myBinN+1);
    string en = "eta";
    string star = "#eta*";
    string temp;
    string temp2;
    for (unsigned i = 0; i < splithists.size(); ++i)
    {
    	temp = en+std::to_string(i);
    	temp2  = star+"="+std::to_string(mybins[i])+" #eta distribution";
    	splithists[i] = new TH1F(temp.c_str(),temp2.c_str(),myBinN,mybins);
    	splithists[i]->SetXTitle("#eta");
    	splithists[i]->SetYTitle("Normalized count");
    	splithists[i]->SetMarkerStyle(7);
    	temp = temp+">>"+temp;
    	dijet_tree->Draw(temp.c_str(),"","goff");
    	splithists[i]->Sumw2();
    	if(splithists[i]->Integral()!=0){
    		splithists[i]->Scale(1/splithists[i]->Integral());
            cout<<"Created number: "<<i<<'\n';
    	}
    	else{
    		cout<<"Failed number: "<<i<<'\n';
    	}
    	splithists[i]->SetAxisRange(0,.85,"Y");
    	splithists[i]->Draw();
    	tc->Write();
    }
    output->Close();
    delete output;
}