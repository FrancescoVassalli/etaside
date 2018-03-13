#include "TFile.h"
#include "TCanvas.h"
#include "AtlasStyle.h"
#include "AtlasUtils.h"
#include "TLegend.h"

using namespace std;

namespace nicehists{
	short colors[7]={kRed,kBlue,kGreen+2,kMagenta+3,kOrange+4,kCyan+1,kGray+1};
	short styles[7]={kFullCircle,kOpenSquare,kFullTriangleUp,kFullDiamond,kFullCross,kFullStar,kOpenFourTrianglesX};
	void makeBins(float* bins, int min, int nBins, float width){
		for(int i=0; i<=nBins;++i){
			bins[i] = min + width*i;
		}
	}

	void makeMarkerNice(TH1F** h, int n){
		for (int i = 1; i < n; ++i)
		{
			(*h)->SetMarkerStyle(styles[i-1]);
			(*h)->SetMarkerColor(colors[i-1]);
			h++;
		}
	}
	void makeLineColors(TH1F** h, int n){
		for (int i = 1; i < n; ++i)
		{
			(*h)->SetLineColor(colors[i-1]);
			h++;
		}
	}
	void makeLegend(TLegend* tl, TH1F** h, int n, std::string *titles){
		for (int i = 0; i < n; ++i)
		{
			tl->AddEntry((*h++),titles++->c_str(),"p");
		}
	}
	void makeNiceHist(TH1* h){
		h->SetMarkerStyle(kFullTriangleUp);
	}
}

void Bzmath(){
	/* set up bins */
	const int nBins = 41;
	const float minvalue = -2.5;
	const float binwidth = .1;
	float bins[nBins+1];
	nicehists::makeBins(bins, minvalue,nBins,binwidth);
	//SetAtlasStyle();

	TCanvas *tc = new TCanvas();
	gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

	TFile *input = new TFile("calout.root");
	const int classes=8;
	TH1F *htemp;
	TH1F *plus[classes];
	TH1F *minus[classes];
	std::string strtemp;
	string front = "cal-eta;";
	/* get data from file*/
	int startbin;
	int mybincounter;
	for (int i = 0; i < classes; ++i)
	{
		strtemp = front+ to_string(i+1);
		htemp = (TH1F*) input->Get(strtemp.c_str());
		startbin = htemp->FindBin(bins[0],0,0);
		plus[i] = new TH1F(strtemp.c_str(),"",nBins,bins);
		plus[i]->Sumw2();
		strtemp = "-"+strtemp;
		minus[i] = new TH1F(strtemp.c_str(),"",nBins,bins);
		minus[i]->Sumw2();
		mybincounter=0;
		for (int j = startbin; j <= startbin+nBins; ++j)
		{
			plus[i]->SetBinContent(mybincounter,htemp->GetBinContent(j));
			plus[i]->SetBinError(mybincounter,htemp->GetBinError(j));
			minus[i]->SetBinError(nBins-mybincounter,htemp->GetBinError(j));
			minus[i]->SetBinContent(nBins-mybincounter++,htemp->GetBinContent(j));
		}
		
		//to check 
		plus[i]->Draw();
		plus[i]->Print("All");
		strtemp = front+ to_string(i+1)+".pdf";
		tc->SaveAs(strtemp.c_str());
		strtemp = "-"+strtemp;
		minus[i]->Draw();
		tc->SaveAs(strtemp.c_str());
		
	}
	/* make the emmision function */
	const float wL = 1;
	/* could make these some kind of histogram/ bin to account uncertainty or do it by hand *
	 might need min bias 0-90% eventually */
	float wR[classes] = {3,6.4,8.8,10.4,12,13.6,15.1,17.2}; // wL is subtracted from Npart
	TH1F *function[classes];
	TH1F *htemp2;
	string centrality[classes] = {"60-90\%","40-60\%","30-40\%","20-30\%","10-20\%","5-10\%","1-5\%","0-1\%"};
	for (int i = 0; i < classes; ++i)
	{
		strtemp = "Emmision Function from p-Pb #eta [-2.1,2.1] " + centrality[i];
		function[i] = new TH1F (centrality[i].c_str(),strtemp.c_str(),nBins,bins);
		function[i] ->SetTitle(strtemp.c_str());
		htemp = (TH1F*)plus[i] ->Clone("temp");
		htemp2 =(TH1F*) plus[i] ->Clone("temp2");
		htemp->Add(minus[i],1);
		htemp2->Add(minus[i],-1);
		htemp->Scale(1/(wR[i]+wL));
		htemp2->Scale(1/(wR[i]-wL));
		htemp->Add(htemp2,1);
		htemp->Scale(.5);
		function[i] = htemp;
		strtemp = "function"+to_string(i)+".pdf";
		function[i]->SetYTitle("F(#eta)");
		function[i]->SetXTitle("#eta");
		nicehists::makeNiceHist(function[i]);
		//SetAtlasStyle(function[i]);
		function[i]->Draw("P hist");
		tc->SaveAs(strtemp.c_str());
	}
	delete tc;
	tc = new TCanvas("tc","Centrality",1280,720);
	tc->Divide(4,2,.001,.001);
	for (int i = 0; i < classes; ++i)
	{
		tc->cd(i+1);
		function[i]->SetTitle(centrality[i].c_str());
    	function[i]->SetTitleOffset(1);
		function[i]->SetTitleSize(.05);
		function[i]->SetLabelSize(.05);
    	function[i]->GetYaxis()->SetTitleSize(.05);
    	function[i]->GetYaxis()->SetLabelSize(.05);
    	function[i]->GetYaxis()->SetTitleOffset(1);
		function[i]->SetLabelSize(.07);
		//function[i]->SetRange()
		function[i]->Draw("P hist");
	}
	tc->SaveAs("Allfunctions.pdf");
	delete tc;
	/* start of the function color plot*/
	tc = new TCanvas("tc","Centrality",720,680);
	//tc->Divide(1,2,.001,.001);
	//tc->cd(1);
	nicehists::makeHistColors(function,classes);
	TLegend *tl=new TLegend(.7,.6,.8,.9,"wounded nucleon model p-Pb 5.02 TeV boost corrected");
	nicehists::makeLegend(tl,function,classes,centrality);
	function[0]->SetAxisRange(0,4.09,"Y");
	function[0]->GetYaxis()->SetTitleOffset(0.6);
	//F(#eta) wounded nucleon model p-Pb 5.02 TeV boost corrected
	function[0]->SetTitle("");
	function[0]->SetTitleSize(.1);
	function[0]->SetTitleOffset(.7);
	function[0]->Draw("l hist");
	for (int i = 1; i < classes; ++i)
	{
		function[i]->Draw("l hist same");
	}
	tl->Draw();
	gPad->SetTicky();
	gPad->SetTickx();
	
	//gPad->SetBottomMargin(.2);
	tc->SaveAs("functioncolor.pdf");
}