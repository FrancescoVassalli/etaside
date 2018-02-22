#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include "TChain.h"

using namespace std;

void splithist(){
	TCanvas *tc = new TCanvas();
	TChain *dijet_tree = new TChain("tree100");
	string fname = "etadist.root";
	dijet_tree->Add(fname.c_str());
	gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    const int nBins = 53;
    double bins[nBins+1]= {-2.65,-2.55,-2.45,-2.35,-2.25,-2.15,-2.05,-1.95,-1.85,-1.75,-1.65,-1.55,-1.45,-1.35,-1.25,-1.15,-1.05,-0.95,-0.85,-0.75,-0.65,-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95,2.05,2.15,2.25,2.35,2.45,2.55,2.65};
    for(int i=0; i<=nBins;i++){
    	bins[i] = bins[i]-.05; //accidentally put the bin middles in instead of the bin lows
    }

    std::vector<TH1F*> splithists(nBins+1);
    string en = "eta";
    string star = "#eta*";
    string temp;
    string temp2;
    for (unsigned i = 0; i < splithists.size(); ++i)
    {
    	temp = en+std::to_string(i);
    	temp2  = star+"="+std::to_string(bins[i])+" #eta distribution";
    	splithists[i] = new TH1F(temp.c_str(),temp2.c_str(),nBins,bins);
    	splithists[i]->SetXTitle("#eta");
    	splithists[i]->SetYTitle("Normalized count");
    	splithists[i]->SetMarkerStyle(7);
    	temp = temp+">>"+temp;
    	dijet_tree->Draw(temp.c_str(),"","goff");
    	splithists[i]->Sumw2();
    	if(splithists[i]->Integral()!=0){
    		splithists[i]->Scale(1/splithists[i]->Integral());
    	}
    	splithists[i]->SetAxisRange(0,.85,"Y");
    	splithists[i]->Draw();
    	temp = en+std::to_string(i)+".pdf";
    	tc->SaveAs(temp.c_str());
    }
    TH1F *etastar = new TH1F("eta*","#eta* distribution",nBins,bins);
    etastar->SetXTitle("#eta*");
    etastar->SetYTitle("Normalized count");
	etastar->SetMarkerStyle(7);
	dijet_tree->Draw("etaboost>>eta*","","goff");
	etastar->Sumw2();

	TH1F *etaproper = new TH1F("etaproper","#eta distribution",nBins,bins);
	etaproper->SetXTitle("#eta");
    etaproper->SetYTitle("Normalized count");
	etaproper->SetMarkerStyle(7);
	dijet_tree->Draw("eta>>etaproper","","goff");
	etaproper->Sumw2();

	TH1F *etasolve = new TH1F("solve","calculated #eta distribution",nBins,bins);
	etasolve->SetXTitle("#eta");
    etasolve->SetYTitle("Normalized count");
	etasolve->SetMarkerStyle(7);
	etasolve->Sumw2();
	float caletavals[nBins];
	for (int i = 0; i <= nBins; ++i)
	{
		etasolve->SetBinContent(i,0);
	}
	for (int i = 0; i <= nBins; ++i)
	{
		splithists[i]->Scale(etastar->GetBinContent(i));
		cout<<"bin1: "<<splithists[i]->GetBinContent(1)<<'\n';
		for (int j = 0; j <= nBins; ++j)
		{
			etasolve->SetBinContent(j, etasolve->GetBinContent(j)+splithists[i]->GetBinContent(j));

			etasolve->SetBinError(j,TMath::Power(etasolve->GetBinError(j)*
				etasolve->GetBinError(j)+splithists[i]->GetBinError(j)*splithists[i]->GetBinError(j),0.5));
		}
		cout<<"error: "<<splithists[i]->GetBinError(nBins)<<'\n';
	}
	etasolve->Scale(1/etasolve->Integral());
	//etasolve->SetAxisRange(0.01,.026,"Y");
	etasolve->Draw();
	tc->SaveAs("calced.pdf");
	etaproper->Scale(1/etaproper->Integral());
	etaproper->SetAxisRange(0.01,.026,"Y");
	etaproper->Draw();
	tc->Print("proper.pdf");

	TH1F *residuals = new TH1F("res","calculation residual",nBins,bins);
	residuals->Add(etasolve,etaproper,1,-1);
	residuals->Draw();
	tc->SaveAs("residuals.pdf");
}