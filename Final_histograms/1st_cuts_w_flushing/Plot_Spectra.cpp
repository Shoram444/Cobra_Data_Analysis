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

	TH1F* 				h[n_of_hist];
	TH1F* tot = new TH1F();

	vector<const char*>	h_names;

	TCanvas* 	c 	= new TCanvas("c", "Energy", 1000, 600);
	THStack* 	hs  = new THStack("Histograms","Deposited Energy; Energy[keV]; Counts[#]");
	TFile* 		tf  = new TFile("Final_histograms/1st_cuts_w_flushing/CO_event-Tchain_years.root");
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
			tot->GetYaxis()->SetTitle("Counts [# /2keV]");

		}
		else
		{
			tot->Add(h[y]);
		}

		hs->Add(h[y]);

		year_for_title++;
	}	

	c->cd();
	hs->Add(tot);
	hs->SetNameTitle("Deposited Energy","Deposited Energy");
	gPad->SetGrid();
	gPad->SetLogy();

	hs->Draw("nostack");
	hs->GetXaxis()->SetTitle("Energy [ keV ]"); 
	hs->GetYaxis()->SetTitle("Counts [ # ]"); 
	c->Modified();
	c->Update();

	gPad->BuildLegend(0.85,0.55,0.98,0.75,"");

	ct->cd();
	gPad->SetGrid();
	gPad->SetLogy();
	ct->Update();

	tot->Draw();
	tot->GetXaxis()->SetRangeUser(0.01, 3000);


	ct->SaveAs("Total.png");

	TFile* tff = new TFile("Plots_from_CO_Event-1.root", "NEW");
	c->Write();
	hs->Write();
	c->SaveAs("Hstack.png");

	TCanvas* cy[n_of_hist];
	int j = 0;

	for(int i = 2013; i < 2020; i++ )
	{

		stringstream ss_str_name;
		ss_str_name << "_" << i;
		string strName = ss_str_name.str();

		TString png_name = TString::Format("%s.png", strName.c_str() );

		cy[j]  = new TCanvas(strName.c_str(), strName.c_str(), 800, 600);

		cy[j]->cd();
		h[j]->Draw();

		cy[j]->SaveAs(png_name);

		j++;
	}




	// delete tf;


  // gROOT->ProcessLine(".q");

}