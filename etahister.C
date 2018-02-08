#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include "TChain.h"

using namespace std;

void etahister(){
	TCanvas *tc = new TCanvas();
	TChain *dijet_tree = new TChain("tree100");
	string fname = "etadist.root";
	dijet_tree->Add(fname.c_str());
	gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    TLegend *tl=new TLegend(.4,.2,.6,.4);

    TH1F *h_etanorm = new TH1F("enorm","Eta distribution 5.06 TeV charged final state",100,-5,5);
    h_etanorm->SetXTitle("Eta");
    h_etanorm->SetYTitle("norm count");
    dijet_tree->Draw("eta>>enorm","","goff");
    h_etanorm->Sumw2();
    h_etanorm->Scale(1/h_etanorm->Integral(),"width");
    h_etanorm->SetMarkerStyle(20);
    h_etanorm->SetLineWidth(3);
    tl->AddEntry(h_etanorm, "unboosted","p");
    h_etanorm->Draw("Pe");
    tc->Print("etadistnorm.pdf");

    TH1F *h_etaboost = new TH1F("eboost","Eta distribution 5.06",100,-5,5);
    h_etaboost->SetXTitle("Eta of boosted charged final state");
    h_etaboost->SetYTitle("boosted count");
    dijet_tree->Draw("etaboost>>eboost","","goff");
    h_etaboost->Sumw2();
    h_etaboost->Scale(1/h_etaboost->Integral(),"width");
    h_etaboost->SetMarkerStyle(20);
    tl->AddEntry(h_etaboost,"B= (0,0,.434)","p");
    h_etaboost->SetMarkerColor(kCyan+2);
    h_etaboost->SetLineWidth(3);
    h_etaboost->SetLineColor(kRed);
    h_etaboost->Draw("P same");
    tl->Draw("same");
    tc->Print("etadistboost.pdf");

}