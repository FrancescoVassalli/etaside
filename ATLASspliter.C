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
/*std::vector<TH1F*> clonelist(std::vector<TH1F*> v,std::string name){
	std::vector<TH1F*> r(0);
	std::string temp;
	for(unsigned i=0; i<v.size();i++){
		temp = name+std::to_string(i);
		r.push_back((TH1F*)v[i]->Clone(temp.c_str()));
	}
	return r;
}*/

void ATLASspliter(){
	TCanvas *tc = new TCanvas();
	TFile *output = new TFile("bindists.root");
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

    std::queue<TH1F*> splithists;
    /* get the hists from the file*/
    string temp;
    string temp2 ="eta";
    int loadcounter=14;
    TH1F* htemp=NULL;
    while (loadcounter<=67){
    	temp = temp2+std::to_string(loadcounter)+";1";
    	splithists.push((TH1F*) output->Get(temp.c_str()));
    	splithists.back()->Draw();
    	loadcounter++;
    }
   
	float bintemp;
	double errortemp;
	string en = "eta";
	string star = "eta*";
	std::queue<TH1F*> smallsplit;
	en = en+"2-";
	int tempint;
	loadcounter=14;
	/* make the smaller hists*/
	while(loadcounter<=67){
		temp  = en +std::to_string(loadcounter);
		temp2 = star+"="+std::to_string(-2.7+.1*(loadcounter-14));
		htemp = new TH1F(temp.c_str(),temp2.c_str(),nBins,bins);
		htemp->Sumw2();
		for(int bincounter=14;bincounter<=67;bincounter++){
			if (splithists.front()==NULL)
			{
				cout<<"Front is NULL"<<std::endl;
			}
			htemp->SetBinContent(bincounter-13,splithists.front()->GetBinContent(bincounter));
			htemp->SetBinError(bincounter-13,splithists.front()->GetBinError(bincounter));
		}
		splithists.pop();

		htemp->Draw("e1");
		temp = temp+".pdf";
		//tc->SaveAs(temp.c_str());
		smallsplit.push(htemp);
		loadcounter++;
	}
	/*set the ATLAS data*/
	std::vector<TH1F*> atlasstar(8);
	float d90[nBins+1] = {8.681,8.556,8.791,8.889,8.937,8.969,8.993,9.109,9.059,9.093,9.093,9.03,9.016,8.975,8.898,8.868,8.825,8.698,8.672,8.585,8.524,8.456,8.375,8.308,8.325,8.285,8.253,8.214,8.284,8.358,8.415,8.473,8.609,8.686,8.723,8.839,8.955,8.993,9.089,9.212,9.221,9.265,9.34,9.41,9.405,9.405,9.433,9.444,9.372,9.489,9.474,9.34,9.364,9.23};
	float d60[nBins+1] = {17.875,17.76,18.213,18.377,18.472,18.499,18.628,18.596,18.494,18.582,18.48,18.353,18.292,18.125,18.013,17.827,17.725,17.564,17.317,17.067,16.925,16.781,16.599,16.404,16.29,16.246,16.127,16.125,16.062,16.159,16.266,16.337,16.445,16.639,16.697,16.842,16.913,17.014,17.103,17.156,17.251,17.275,17.306,17.331,17.293,17.192,17.143,17.104,17.049,17.036,16.869,16.8,16.576,16.135};
	float d40[nBins+1] = {24.485,24.189,24.759,24.921,25.078,25.159,25.101,25.24,25.184,25.211,24.957,24.86,24.674,24.552,24.339,24.068,23.834,23.536,23.29,23.084,22.681,22.498,22.199,22.062,21.765,21.624,21.528,21.304,21.475,21.436,21.508,21.642,21.785,21.733,21.977,21.947,22.155,22.258,22.296,22.294,22.343,22.393,22.369,22.368,22.077,22.062,21.88,21.743,21.565,21.402,21.196,21.07,20.773,20.64};
	float d30[nBins+1] = {29.805,29.569,30.183,30.321,30.372,30.566,30.435,30.471,30.323,30.511,30.27,30.141,29.773,29.457,29.231,29.046,28.584,28.325,27.967,27.605,27.229,26.995,26.64,26.236,26.064,25.878,25.661,25.43,25.486,25.336,25.48,25.575,25.747,25.762,25.93,26.057,26.028,26.156,26.274,26.053,26.161,26.094,25.825,25.814,25.476,25.436,25.044,24.878,24.585,24.273,24.073,23.859,23.472,23.418};
	float d20[nBins+1] = {36.899,36.447,37.182,37.246,37.477,37.535,37.64,37.514,37.402,37.346,37.036,37.079,36.613,36.275,35.884,35.486,35.233,34.59,34.174,33.723,33.193,32.891,32.394,32.008,31.564,31.321,30.884,30.986,30.761,30.725,30.727,30.705,30.858,30.981,30.961,31.002,31.11,31.168,31.136,31.13,30.994,30.783,30.688,30.425,30.155,29.798,29.52,29.112,28.749,28.395,28.243,27.542,27.18,26.831};
	float d10[nBins+1] ={45.076,44.758,45.486,45.405,45.476,45.737,45.67,45.736,45.468,45.241,45.229,44.734,44.218,43.759,43.496,42.806,42.27,41.797,41.235,40.657,39.928,39.271,39.001,38.332,37.806,37.261,36.838,36.932,36.576,36.64,36.587,36.375,36.542,36.596,36.401,36.606,36.469,36.612,36.459,36.384,36.161,35.915,35.594,35.476,34.928,34.54,34.21,33.818,33.3,32.685,31.913,31.815,31.143,30.302};
	float d5[nBins+1]={54.735,54.505,55.18,55.167,55.363,55.175,55.188,55.37,55.091,54.507,54.545,54.103,53.203,53.042,52.116,51.792,50.974,50.199,49.519,48.815,48.099,46.757,46.5,45.807,45.057,44.566,44.065,43.732,43.334,43.224,43.258,43.009,42.975,42.884,42.871,42.775,42.799,42.744,42.713,42.374,42.285,41.803,41.337,41.002,40.531,39.894,39.241,38.439,37.792,37.272,36.504,35.629,34.976,33.831};
	float d1[nBins+1] = {72.722,72.305,72.463,72.211,72.94,72.589,72.3,71.849,72.646,71.751,70.867,70.164,69.717,68.356,67.299,66.848,64.957,64.473,63.135,62.982,61.85,60.648,59.145,58.71,57.403,56.947,55.909,55.101,54.628,54.621,54.44,54.457,54.194,53.919,53.581,53.583,53.539,53.086,52.688,52.625,51.795,51.467,50.868,50.302,49.617,48.78,47.307,46.313,45.777,44.7,43.795,43.003,41.486,40.29};
	string atlasstr = "ATLAS";
	float e90[nBins+1]= {0.098,0.087,0.088,0.093,0.093,0.096,0.096,0.097,0.097,0.097,0.097,0.097,0.096,0.096,0.095,0.095,0.095,0.094,0.094,0.093,0.093,0.093,0.089,0.089,0.09,0.09,0.09,0.089,0.09,0.091,0.091,0.092,0.092,0.092,0.092,0.093,0.094,0.093,0.094,0.094,0.094,0.094,0.096,0.097,0.097,0.097,0.097,0.098,0.098,0.095,0.094,0.098,0.099,0.172};
	float e60[nBins+1]= {0.136,0.121,0.123,0.129,0.13,0.134,0.135,0.134,0.134,0.135,0.135,0.134,0.133,0.133,0.132,0.131,0.131,0.13,0.129,0.127,0.127,0.126,0.126,0.124,0.124,0.124,0.123,0.124,0.123,0.124,0.124,0.125,0.126,0.127,0.127,0.128,0.128,0.128,0.128,0.128,0.128,0.128,0.131,0.131,0.131,0.13,0.13,0.13,0.131,0.126,0.123,0.13,0.13,0.223};
	float e40[nBins+1]= {0.203,0.18,0.182,0.19,0.192,0.199,0.198,0.199,0.2,0.201,0.2,0.2,0.198,0.198,0.197,0.196,0.195,0.193,0.192,0.191,0.188,0.188,0.186,0.185,0.184,0.183,0.182,0.181,0.183,0.182,0.183,0.184,0.187,0.185,0.187,0.187,0.188,0.188,0.188,0.187,0.188,0.188,0.192,0.192,0.189,0.19,0.189,0.189,0.189,0.181,0.176,0.185,0.185,0.32};
	float e30[nBins+1]= {0.238,0.212,0.214,0.223,0.224,0.233,0.233,0.233,0.232,0.235,0.235,0.235,0.231,0.23,0.229,0.229,0.227,0.226,0.224,0.223,0.22,0.218,0.22,0.216,0.216,0.215,0.214,0.213,0.213,0.211,0.213,0.213,0.218,0.218,0.219,0.219,0.219,0.219,0.219,0.217,0.218,0.217,0.221,0.22,0.218,0.218,0.216,0.216,0.215,0.204,0.199,0.21,0.209,0.366};
	float e20[nBins+1]= {0.286,0.254,0.256,0.267,0.269,0.279,0.281,0.28,0.281,0.282,0.282,0.283,0.279,0.278,0.278,0.275,0.276,0.272,0.27,0.268,0.265,0.263,0.266,0.264,0.26,0.26,0.256,0.257,0.256,0.256,0.255,0.255,0.262,0.262,0.262,0.262,0.263,0.262,0.262,0.262,0.259,0.258,0.266,0.263,0.261,0.259,0.258,0.256,0.255,0.242,0.236,0.245,0.247,0.426};
	float e10[nBins+1]= {0.371,0.334,0.336,0.349,0.35,0.365,0.366,0.366,0.367,0.367,0.37,0.369,0.364,0.363,0.362,0.359,0.358,0.356,0.353,0.35,0.345,0.34,0.347,0.343,0.339,0.334,0.332,0.333,0.33,0.331,0.33,0.328,0.339,0.339,0.337,0.339,0.338,0.34,0.337,0.334,0.333,0.332,0.341,0.34,0.336,0.331,0.331,0.33,0.329,0.31,0.298,0.317,0.315,0.544};
	float e5[nBins+1]= {0.509,0.465,0.467,0.49,0.49,0.507,0.511,0.511,0.513,0.512,0.516,0.515,0.506,0.51,0.503,0.504,0.5,0.497,0.491,0.486,0.483,0.472,0.487,0.48,0.476,0.472,0.467,0.465,0.461,0.46,0.459,0.457,0.472,0.472,0.47,0.471,0.471,0.47,0.47,0.464,0.465,0.46,0.47,0.467,0.464,0.458,0.454,0.445,0.443,0.422,0.407,0.429,0.427,0.733};
	float e1[nBins+1]= {0.529,0.498,0.497,0.519,0.526,0.541,0.543,0.54,0.551,0.549,0.547,0.545,0.544,0.54,0.532,0.535,0.522,0.521,0.513,0.517,0.51,0.502,0.511,0.51,0.501,0.499,0.494,0.49,0.485,0.481,0.483,0.48,0.498,0.495,0.493,0.492,0.491,0.488,0.483,0.483,0.474,0.472,0.489,0.485,0.482,0.474,0.464,0.455,0.456,0.432,0.424,0.449,0.443,0.825};
	/* take the small histograms and make transformations for each centrality class*/
	std::vector<TH1F*> transformations[atlasstar.size()];
	for(unsigned i=0; i<atlasstar.size();i++){
		temp = atlasstr+std::to_string(i);
		temp2 = star+" "+atlasstr;
		htemp = new TH1F(temp.c_str(),temp2.c_str(),nBins,bins);
		htemp->Sumw2();
		atlasstar[i] = htemp;
	}
	int clonecounter =0;
	while(!smallsplit.empty()){ // load the small splits into transformations
		// need to make clones
		temp = "90-"+std::to_string(clonecounter);
		transformations[0].push_back((TH1F*) smallsplit.front()->Clone(temp.c_str()));
		temp = "60-"+std::to_string(clonecounter);
		transformations[1].push_back((TH1F*) smallsplit.front()->Clone(temp.c_str()));
		temp = "40-"+std::to_string(clonecounter);
		transformations[2].push_back((TH1F*) smallsplit.front()->Clone(temp.c_str()));
		temp = "30-"+std::to_string(clonecounter);
		transformations[3].push_back((TH1F*) smallsplit.front()->Clone(temp.c_str()));
		temp = "20-"+std::to_string(clonecounter);
		transformations[4].push_back((TH1F*) smallsplit.front()->Clone(temp.c_str()));
		temp = "10-"+std::to_string(clonecounter);
		transformations[5].push_back((TH1F*) smallsplit.front()->Clone(temp.c_str()));
		temp = "5-"+std::to_string(clonecounter);
		transformations[6].push_back((TH1F*) smallsplit.front()->Clone(temp.c_str()));
		temp = "1-"+std::to_string(clonecounter);
		transformations[7].push_back((TH1F*) smallsplit.front()->Clone(temp.c_str()));
		smallsplit.pop();
		clonecounter++;
	}
	/* load the ATALS data into the hists*/
	for (int i = 0; i <= nBins; ++i) 
	{
		atlasstar[0]->SetBinContent(i,d90[i]);
		atlasstar[1]->SetBinContent(i,d60[i]);
		atlasstar[2]->SetBinContent(i,d40[i]);
		atlasstar[3]->SetBinContent(i,d30[i]);
		atlasstar[4]->SetBinContent(i,d20[i]);
		atlasstar[5]->SetBinContent(i,d10[i]);
		atlasstar[6]->SetBinContent(i,d5[i]);
		atlasstar[7]->SetBinContent(i,d1[i]);
		atlasstar[0]->SetBinError(i,e90[i]);
		atlasstar[1]->SetBinError(i,e60[i]);
		atlasstar[2]->SetBinError(i,e40[i]);
		atlasstar[3]->SetBinError(i,e30[i]);
		atlasstar[4]->SetBinError(i,e20[i]);
		atlasstar[5]->SetBinError(i,e10[i]);
		atlasstar[6]->SetBinError(i,e5[i]);
		atlasstar[7]->SetBinError(i,e1[i]);
	}
	/* clone the unmodified data*/
	std::vector<TH1F*> aClone(8);
	for(unsigned i=0; i<aClone.size();++i){
		aClone[i] = (TH1F*) atlasstar[i]->Clone("unmodified");
	}
	/* weight the transformations */
	for(int i=0; i<=nBins;i++)	
	{
		transformations[0][i]->Scale(atlasstar[0]->GetBinContent(i));
		transformations[1][i]->Scale(atlasstar[1]->GetBinContent(i));
		transformations[2][i]->Scale(atlasstar[2]->GetBinContent(i));
		transformations[3][i]->Scale(atlasstar[3]->GetBinContent(i));
		transformations[4][i]->Scale(atlasstar[4]->GetBinContent(i));
		transformations[5][i]->Scale(atlasstar[5]->GetBinContent(i));
		transformations[6][i]->Scale(atlasstar[6]->GetBinContent(i));
		transformations[7][i]->Scale(atlasstar[7]->GetBinContent(i));
	}

	/*use the transformations to make the untransformed hists for each centrality*/
	std::vector<TH1F*> calAT(8);
	for (int i = 0; i < 8; ++i)
	{
		temp = "solve" + std::to_string(i);
		calAT[i] = new TH1F(temp.c_str(),"Calculated #eta Distribution",nBins,bins);
		calAT[i]->Sumw2();
		calAT[i]->SetXTitle("#eta");
		calAT[i]->SetYTitle("Count");
		for (int j = 0; j <= nBins; ++j)
		{
			bintemp=0;
			errortemp=0;
			/* for each bin in the hist sum the total of that bin from all the transformations*/
			for(int k=0; k<=nBins;k++){
				bintemp+= transformations[i][k]->GetBinContent(j);
				errortemp = addError(errortemp,transformations[i][k]->GetBinError(j));
			}
			calAT[i]->SetBinContent(j,bintemp);
			calAT[i]->SetBinError(j, errortemp);
		}
		calAT[i]->SetMarkerStyle(7);
		calAT[i]->SetAxisRange(-2.7,2.6,"X");
		calAT[i]->Draw();
		temp = temp+".pdf";
		tc->SaveAs(temp.c_str());
	}
	TH1F* hcompare= (TH1F*) calAT[0]->Clone("compare");
	/* draw the untransformed hists at each centrality*/
	delete tc;
	tc = new TCanvas("tc","Centrality",1280,720);
	tc->Divide(4,2,.001,.001);
	calAT[0]->SetTitle("60-90\%");
	calAT[1]->SetTitle("40-60\%");
	calAT[2]->SetTitle("30-40\%");
	calAT[3]->SetTitle("20-30\%");
	calAT[4]->SetTitle("10-20\%");
	calAT[5]->SetTitle("5-10\%");
	calAT[6]->SetTitle("1-5\%");
	calAT[7]->SetTitle("0-1\%");
	TLegend *tl=new TLegend(.35,.2,.66,.4);
	for (int i = 0; i < 8; ++i)
	{
		tc->cd(i+1);
    	calAT[i]->SetTitleOffset(.7);
		calAT[i]->SetTitleSize(.05);
		calAT[i]->SetLabelSize(.05);
    	calAT[i]->GetYaxis()->SetTitleSize(.05);
    	calAT[i]->GetYaxis()->SetLabelSize(.05);
    	calAT[i]->GetYaxis()->SetTitleOffset(1);
		calAT[i]->SetLabelSize(.07);
		calAT[i]->SetMarkerStyle(kFullSquare);
		calAT[i]->Draw("P hist");
		aClone[i]->SetMarkerStyle(kOpenCircle);
		aClone[i]->SetMarkerColor(kRed);
		aClone[i]->Draw("P hist same");
	}
	tc->cd(1);
	tl->AddEntry(calAT[0],"#eta","p");
	tl->AddEntry(aClone[0],"#eta*","p");
	tl->Draw();
	tc->SaveAs("Atlas p-Pb eta.pdf");
	tc->SaveAs("Atlaseta.root");
	for (int i = 0; i < 8; ++i)
	{
    	calAT[i]->GetYaxis()->SetTitle("Ratio");
		tc->cd(i+1);
		calAT[i]->Divide(aClone[i]);
		calAT[i]->Draw("P hist");
	}
}