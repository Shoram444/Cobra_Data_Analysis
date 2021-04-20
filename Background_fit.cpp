#include <iostream>
#include <string>
#include <math.h> 
#include <vector>
#include "TROOT.h"  
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"

using namespace std;

Double_t F_Exponential(Double_t* x, Double_t* par)
{
	// par[0] - N0
	// par[1] - Lambda
	// par[2] - G1N0
	// par[3] - G1mu
	// par[4] - G1sigma
	// par[5] - G2N0
	// par[6] - G2mu
	// par[7] - G2sigma


	Double_t expo    = 0.0;
	Double_t gauss1  = 0.0;
	Double_t gauss2  = 0.0;

	Double_t arg1 = (par[4] != 0.0) ? (x[0] - par[3]) / (par[4]) : 0.0;
	Double_t arg2 = (par[7] != 0.0) ? (x[0] - par[6]) / (par[7]) : 0.0;
	
	gauss1 = exp(-0.5 * arg1 * arg1) / (par[4] * sqrt(2.0 * TMath::Pi()));
	gauss2 = exp(-0.5 * arg2 * arg2) / (par[7] * sqrt(2.0 * TMath::Pi()));

	expo = par[0]*exp(-par[1]*x[0]) + par[2]*gauss1 + par[5]*gauss2; 
	return expo;
}	

void Background_fit()
{

	double Q_116Cd 			= 2813.4;
	double Sigma_at_116Cd 	=  48.40;  //average resolution at Q
	double k 				=   50.0;  //average resolution at Q
	int    merged_bins 		=     50; 
	double E_Fit_min 		= 1800.0; 
	double E_Fit_max 		= 3200.0; 



	const char* fname 	= "/home/shoram/Work/Diploma_Thesis/Cobra_Data_Analysis/Final_histograms/1st_cuts_w_flushing/Total_years_root.root";
	TFile* tf 			= new TFile(fname);
	TH1F*   h 			= (TH1F*) tf->Get("h_2013");

	h->Rebin(merged_bins);
	h->GetXaxis()->SetRangeUser(E_Fit_min, E_Fit_max);

	double bin_width = h->GetBinWidth(1);
	cout << bin_width << endl;


	TF1* EFit = new TF1("ExpoFit", F_Exponential, E_Fit_min, E_Fit_max, 2);
	EFit->SetParNames("N0", "Lambda", "G1N0", "G1mu", "G1sigma", "G2N0", "G2mu", "G2sigma");
	EFit->SetParameters( 23500 , 2.62e-3, 17000, 511, 10, 41000, 610, 10);

	// EFit->SetParLimits(2, 15000, 19000);
	// EFit->SetParLimits(3, 510, 520);
	// EFit->SetParLimits(4, 8, 12);

	// EFit->SetParLimits(5, 38000, 42000);
	// EFit->SetParLimits(6, 590, 620);
	// EFit->SetParLimits(7, 8, 12);

	EFit->SetNpx(2000);
	h->Fit("ExpoFit", "E");
	cout <<" Xi2/ndf: " << EFit->GetChisquare() / EFit->GetNDF() << endl;

	h->Draw();
	EFit->Draw("SAME");

	// f->SetNpx(2000);


}



////////////FIT FOR TOTAL SPECTRUM BETWEEN 400 keV to 1100keV

	// double Q_116Cd 			= 2813.4;
	// double Sigma_at_116Cd 	=  48.40;  //average resolution at Q
	// double k 				=   50.0;  //average resolution at Q
	// int    merged_bins 		=     10; 
	// double E_Fit_min 		=  400.0; 
	// double E_Fit_max 		= 1100.0; 



	// const char* fname 	= "/home/shoram/Work/Diploma_Thesis/Cobra_Data_Analysis/Final_histograms/1st_cuts_w_flushing/Total_years_root.root";
	// TFile* tf 			= new TFile(fname);
	// TH1F*   h 			= (TH1F*) tf->Get("h_2013");

	// h->Rebin(merged_bins);
	// h->GetXaxis()->SetRangeUser(E_Fit_min, E_Fit_max);

	// double bin_width = h->GetBinWidth(1);


	// TF1* EFit = new TF1("ExpoFit", F_Exponential, E_Fit_min, E_Fit_max, 8);
	// EFit->SetParNames("N0", "Lambda", "G1N0", "G1mu", "G1sigma", "G2N0", "G2mu", "G2sigma");
	// EFit->SetParameters( 6.3 , 1.13e-3, 17000, 511, 10, 41000, 610, 10);

	// EFit->SetParLimits(2, 15000, 19000);
	// EFit->SetParLimits(3, 510, 520);
	// EFit->SetParLimits(4, 8, 12);

	// EFit->SetParLimits(5, 38000, 42000);
	// EFit->SetParLimits(6, 590, 620);
	// EFit->SetParLimits(7, 8, 12);

	// EFit->SetNpx(2000);
	// h->Fit("ExpoFit", "E");
	// cout <<" Xi2/ndf: " << EFit->GetChisquare() / EFit->GetNDF() << endl;

	// h->Draw();
	// EFit->Draw("SAME");

	// // f->SetNpx(2000);

////////////FIT FOR TOTAL SPECTRUM BETWEEN 400 keV to 1100keV


////////////FIT FOR TOTAL SPECTRUM BETWEEN 1800 keV to 3200keV


	// double Q_116Cd 			= 2813.4;
	// double Sigma_at_116Cd 	=  48.40;  //average resolution at Q
	// double k 				=   50.0;  //average resolution at Q
	// int    merged_bins 		=     50; 
	// double E_Fit_min 		= 1800.0; 
	// double E_Fit_max 		= 3200.0; 



	// const char* fname 	= "/home/shoram/Work/Diploma_Thesis/Cobra_Data_Analysis/Final_histograms/1st_cuts_w_flushing/Total_years_root.root";
	// TFile* tf 			= new TFile(fname);
	// TH1F*   h 			= (TH1F*) tf->Get("h_2013");

	// h->Rebin(merged_bins);
	// h->GetXaxis()->SetRangeUser(E_Fit_min, E_Fit_max);

	// double bin_width = h->GetBinWidth(1);


	// TF1* EFit = new TF1("ExpoFit", F_Exponential, E_Fit_min, E_Fit_max, 2);
	// EFit->SetParNames("N0", "Lambda", "G1N0", "G1mu", "G1sigma", "G2N0", "G2mu", "G2sigma");
	// EFit->SetParameters( 23500 , 2.62e-3, 17000, 511, 10, 41000, 610, 10);

	// // EFit->SetParLimits(2, 15000, 19000);
	// // EFit->SetParLimits(3, 510, 520);
	// // EFit->SetParLimits(4, 8, 12);

	// // EFit->SetParLimits(5, 38000, 42000);
	// // EFit->SetParLimits(6, 590, 620);
	// // EFit->SetParLimits(7, 8, 12);

	// EFit->SetNpx(2000);
	// h->Fit("ExpoFit", "E");
	// cout <<" Xi2/ndf: " << EFit->GetChisquare() / EFit->GetNDF() << endl;

	// h->Draw();
	// EFit->Draw("SAME");

	// // f->SetNpx(2000);

////////////FIT FOR TOTAL SPECTRUM BETWEEN 1800 keV to 3200keV
