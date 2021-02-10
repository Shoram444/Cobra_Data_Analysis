#include <iostream>
#include <dirent.h>
#include <string>
#include <math.h> 
#include <vector>
#include <stdio.h>
#include "TROOT.h"  
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include <limits.h>
#include "TSystem.h"

using namespace std;

const int n_of_hist = 7;

void Plot_Spectra() 
{

	TH1F* 		h[n_of_hist];
	TH1F*				 tot;

	vector<const char*>	h_names;

	TCanvas* 	c 	= new TCanvas("c", "Energy", 1000, 600);
	THStack* 	hs  = new THStack("Histograms","Deposited Energy; Energy[keV]; Counts[#]");
	TFile* 		tf  = new TFile("Processed_data.root");
	TCanvas* 	ct 	= new TCanvas("ct", "Tot_Energy", 1000, 600);


	int year_for_title = 2013;
	for(int y = 0; y < n_of_hist; y++)
	{
		stringstream 			 	   h_title;			//names for histograms		
		h_title  <<  "h_" <<   	year_for_title;
		string h_tit_str  = 	 h_title.str();

		h[y] = (TH1F*) tf->Get(h_tit_str.c_str());
		h[y]->SetLineColor(1+y);

		if(y ==0)
		{
			tot = (TH1F*) h[0]->Clone();
			tot->SetTitle("Total");
			tot->SetLineColor(8);
		}
		else
		{
			tot->Add(h[y]);
		}

		hs->Add(h[y]);

		year_for_title++;
	}	

	hs->Add(tot);
	gPad->SetGrid();
	gPad->SetLogy();

	hs->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	ct->cd();
	tot->Draw();

	ct->SaveAs("Total_spectrum_without_flush.png");
  gROOT->ProcessLine(".q");
	
}