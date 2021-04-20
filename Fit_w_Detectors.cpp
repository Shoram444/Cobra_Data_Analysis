#include <iostream>
#include <dirent.h>
#include <fstream>

#include <string>
#include <math.h> 
#include <vector>
#include <stdio.h>
#include "TROOT.h"  
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include <limits.h>
#include </home/shoram/ProgramFiles/boost_1_75_0/boost/algorithm/string/find.hpp>
#include "TSystem.h"
#include "/home/shoram/Work/Diploma_Thesis/CO_event/include/CO_event.hh"
#include "/home/shoram/Work/Diploma_Thesis/CO_detector/include/CO_detector.hh"


using namespace boost;
using namespace std;

R__LOAD_LIBRARY(/home/shoram/Work/Diploma_Thesis/CO_event/lib/libCO_event.so);
R__LOAD_LIBRARY(/home/shoram/Work/Diploma_Thesis/CO_detector/lib/libCO_detector.so);




const char* 		   runs = "runs";

const double     total_expo = 1.25194;


const double 	E_Fit_min_1 =  		400.0;  ///PARAMETERS OF FIT Region 1
const double 	E_Fit_max_1 = 	   1100.0; 
const int 		Bin_Width_1 = 		   20;
const double 		   N0_1 = 2.35078e+04;   
const double 	   lambda_1 = 2.61635e-03;


const double 	E_Fit_min_2 =  	   1800.0;  ///PARAMETERS OF FIT Region 1
const double 	E_Fit_max_2 = 	   3200.0; 
const int 		Bin_Width_2 = 		  100;
const double 		   N0_2 = 1.95084e+04;
const double 	   lambda_2 = 1.03512e-03;


vector<string> ReadFiles(bool _defaultPath = true);
TChain* Make_TChain(vector<string>  	 root_file_path, const char* _ttree);
char* str_to_char(string _str);
vector<string>* ListFiles(string _path, string _key);
vector<string> ReadFiles(bool _defaultPath = true);
vector<CO_detector>  Fill_CO_detector(TChain* _T_Chain);




TChain* Make_TChain(vector<string>  	 root_file_path, const char* _ttree)
{

	TChain* cData = new TChain(_ttree);
  
	for(unsigned int i = 0; i < root_file_path.size(); i++)
		{
			char* root_file;

			root_file = str_to_char(root_file_path.at(i));

			if(i%1000 == 0 )
			{
				cout<< i <<" of " 		<< root_file_path.size() <<  " Files Read" 	<< endl;
			}

			cData->Add(root_file);
		}

  return cData;  
}


char* str_to_char(string _str)
{
	int 	n 	= _str.length();
	char* _char = new char[n+1];

    strcpy(_char, _str.c_str());
	
	return _char;
}


vector<string>* ListFiles(string _path, string _key)
{
	vector<string> *file_list = new vector<string>();

	DIR 		   *di;
	struct dirent *dir;
	char* 		c_path; //conversion of string to char. opendir function takes only char. 

	c_path 	= str_to_char(_path);
	di 		= opendir(c_path);

	if (di)
    {
        while ((dir = readdir(di)) != NULL)
        {
            stringstream 			  dir_name;		//backward conversion of char to string.
    		dir_name  	  		<< dir->d_name;
    		string strname  = 	dir_name.str();

    		if(strname == ".." || strname==".")
    		{
    			continue;
    		}

            size_t found = strname.find(_key);	 //checking if "years_ksey" is present in the string. 
		  	if (found!=std::string::npos){file_list->push_back(dir->d_name);}
        }
        closedir(di);
    }
    else 
    {
        perror ("Failed to read dictionaries! Check whether the folders are correctly named!");
        return file_list;
    }

    return file_list;
} 
vector<string> ReadFiles(bool _defaultPath = true)
{
    string 		years_key;
	string 	  folders_key;
	string 		files_key;
	string 	   parent_dir;			// request for path in the console.
	string     fullpath;

	vector<string>  v_fullpath;

	vector<string>*   years;
	vector<string>* folders;
	vector<string>*   files;

	if(_defaultPath)
	{
		parent_dir = "../COBRA_Data";
		years_key  = 			"20";
		folders_key= 		   "cpg";
		files_key  = 		  "root";
	}
	else
	{
	   	cout<< "Please specify the path to Data parent directory: ";
	    cin >> parent_dir;

	    cout<< "Please specify Year of desired data (Use \"20\" for all years): ";
	    cin >> years_key;

		cout<< "Please specify keyword for which folders to scan (Use \"cpg\" for all years): ";
		cin >> folders_key;
		
		cout<< "Please specify keyword for which files to scan (Use \".root\" for all years): ";
		cin >> files_key;
	}

	years = ListFiles(parent_dir, years_key);

	for (int j = 0; j<years->size(); j++)
	{
		folders = ListFiles(parent_dir + "/" + years->at(j), folders_key);

		for(int i = 0; i < folders->size();i++) 
		{
			files = ListFiles(parent_dir + "/" + years->at(j) + "/" + folders->at(i) , files_key);

			for (int k = 0; k < files->size(); k++)
			{
				fullpath        = parent_dir + "/" + years->at(j) + "/" + folders->at(i) + "/" + files->at(k);

				v_fullpath.push_back(fullpath);

			}
		}
	}

	sort(v_fullpath.begin(),v_fullpath.end());

	return v_fullpath;
}


double Get_Exposure(int _d_id, vector<CO_detector>   _v_det, bool _all = false)
{
	CO_detector 	a_det[64];

	double 	total_exposure = 0.0;
	int    	det_Id		   =  _d_id;

	for( unsigned int i = 0 ; i < _v_det.size() ; i++ )
	{
		if(!_all)
		{
			if(_v_det.at(i).get_c_dno() == det_Id)
			{
				double 		c_mas = _v_det.at(i).get_c_mas();
				Double_t	c_dur = _v_det.at(i).get_c_dur();

				a_det[det_Id-1].Add( det_Id, c_mas , c_dur);
			}
		}
		else
		{
			int 		c_dno = _v_det.at(i).get_c_dno();
			double 		c_mas = _v_det.at(i).get_c_mas();
			Double_t	c_dur = _v_det.at(i).get_c_dur();

			a_det[c_dno-1].Add( c_dno, c_mas, c_dur);
		}
	}

	if(!_all)
	{
  		total_exposure = a_det[det_Id-1].calc_exposure_kgy();
	}
	else
	{
  		for( int d = 0; d < 64; d++)			///WATCH OUT! Detector numbering starts from 1, so there has to be n+1 n_of_det
		{
			total_exposure += a_det[d-1].calc_exposure_kgy();
		}
	}

  	return total_exposure;

}

vector<CO_detector>  Fill_CO_detector(TChain* _T_Chain)
{
	vector<CO_detector> 		v_det;


	vector<int>* 	dno = new vector<int>;
	vector<double>* mas = new vector<double>;
	Double_t		dur;
	Double_t		tim;
	vector<double>* p0  = new vector<double>;
	vector<double>* p1  = new vector<double>;
	vector<double>* p2  = new vector<double>;


	_T_Chain->SetBranchAddress("det_no", &dno);
	_T_Chain->SetBranchAddress("det_mass", &mas);
	_T_Chain->SetBranchAddress("duration", &dur);
	_T_Chain->SetBranchAddress("startTime", &tim);
	_T_Chain->SetBranchAddress("det_res_p0", &p0);
	_T_Chain->SetBranchAddress("det_res_p1", &p1);
	_T_Chain->SetBranchAddress("det_res_p2", &p2);

	for( unsigned int i = 0 ; i < _T_Chain->GetEntries() ; i++ )
	{
		if(i%10000==0) cout << i << " of " << _T_Chain->GetEntries() << " Read!" << endl;
		_T_Chain->GetEntry(i);

		for( unsigned int j = 0 ; j < dno->size() ; j++ )
		{
			v_det.emplace_back( dno->at(j), mas->at(j), dur, tim, p0->at(j), p1->at(j), p2->at(j) );
		}

		dno->clear();
		mas->clear();
		p0->clear();
		p1->clear();
		p2->clear();
	}

	return v_det;

}

double chi2NDF(TH1F* _h, int _E_min, int _E_max, double _C, double _lambda)
{
	/// chi2 = \frac{sum_{n=1}^{n}(h_i - C*exp(-\lambda*E_i))}{\sigma^2}

	double  chi2 		= 0.0;
	double  temp_hi    	= 0.0;
	double  temp_Ei    	= 0.0;
	double  temp_sum    = 0.0;
	int     NDF 		=   0;
	double  chi2NDF     = 0.0;
	double  sig2 		= 0.0;

	int 	min_bin_ID = ceil(_E_min / _h->GetBinWidth(1));
	int 	max_bin_ID = ceil(_E_max / _h->GetBinWidth(1));

	NDF = max_bin_ID - min_bin_ID - 3;

	for( int b = min_bin_ID; b <= max_bin_ID; b++)
	{

		temp_hi =     _h->GetBinContent(b) ;
		temp_Ei = 	  _h->GetBinCenter (b) ;
		sig2 	= 	  temp_hi ;
		if( temp_hi != 0 )
		{
			chi2 += pow(( temp_hi - _C * exp( -1 * _lambda * temp_Ei ) ), 2) / sig2;
		}
		if( temp_hi == 0 )
		{
			NDF -= 1;
		}

		// cout 	<< "_h->GetBinContent(b) =  " 	<< _h->GetBinContent(b) << endl
		// 		<< " temp_hi = "				<< temp_hi 				<< endl
		// 		<< " temp_Ei = " 				<< temp_Ei 				<< endl
		// 		<< " temp_sum = "				<< chi2 				<< endl		
		// 		<< " sig2 " 					<< sig2 				<< endl
		// 		<< " NDF = "					<< NDF 					<< endl
		// 		<< " min_bin_ID = "				<< min_bin_ID 			<< endl	
		// 		<< " max_bin_ID = "				<< max_bin_ID 			<< endl	
		// 		<< " -----------------------------------------------"	<< endl	;

	}

	// chi2 	= temp_sum / sig2 ; 
	chi2NDF	= chi2     / NDF  ;

	return chi2NDF;

}


Double_t F_Exponential(Double_t* x, Double_t* par)
{
	// par[0] - "N0_1"
	// par[1] - "Lambda_1"
	// par[2] - "exp_1"
	// par[3] - "Bin_Width_d"

	Double_t expo    = 0.0;

	expo = (par[0]*par[2]*par[3])/exp(par[1]*x[0]); 
	return expo;
}	

void Fit_w_Detectors()
{
	int d_tot = 64;
	TH1F* 	  h[d_tot];

	TCanvas* 	c 	= new TCanvas("c", "c", 1000, 600);
	TFile* 		tf  = new TFile("Final_histograms/1st_cuts_w_flushing/CO_event-TChain_detectors.root");

	vector<string> 			root_file_path 	= 		ReadFiles();
	TChain* 			 	T_Chain_det		=  		Make_TChain(root_file_path, runs);
	vector<CO_detector>     v_det 			= 		Fill_CO_detector(T_Chain_det);

	double 					Exposures[d_tot];

	double tot = 0;

	for(int d = 0; d < 64; d++)
	{
		Exposures[d] = Get_Exposure(d+1, v_det);
		cout << "Exposure of detector " << d+1 << " = " << Exposures[d] << endl;
		tot += Exposures[d];
	}

	for(int d = 0; d < d_tot; d++)
	{
		stringstream 			 	   h_title;			//names for histograms		
		h_title  <<  "h_" <<   	1+d;
		// h_title  <<  "h_" <<   		detector_n;
		string h_tit_str  = 	 h_title.str();

		h[d] = (TH1F*) tf->Get(h_tit_str.c_str());
	}

	double C_d[d_tot];   //C for each detector
	double l_d[d_tot];	 //lambda for each detector

	double temp_A;  //matrix elements
	double temp_B;
	double temp_D;
	double temp_E;
	double temp_F;
	double temp_G;

	double temp_s2; //sigma^2

	double E_min   = 1000;
	double E_max   = 3200;

	int min_bin_ID; 
	int max_bin_ID; 

	for(int d = 0; d < d_tot; d++)
	{
		temp_A = 0.0;  //matrix elements
		temp_B = 0.0;
		temp_D = 0.0;
		temp_E = 0.0;
		temp_F = 0.0;
		temp_G = 0.0;

		min_bin_ID = 0;
		max_bin_ID = 0;

		h[d]->Rebin(50);

		min_bin_ID = ceil(E_min / h[d]->GetBinWidth(1));
		max_bin_ID = ceil(E_max / h[d]->GetBinWidth(1));

		temp_s2 = Exposures[d] * Exposures[d] ;

		for(int b = min_bin_ID; b < max_bin_ID; b++)
		{
			if( h[d]->GetBinContent(b) != 0 )
			{
				temp_A += h[d]->GetBinCenter(b)  	  / temp_s2;                      	// E_i / sigma^2
				temp_B += 					   		1 / temp_s2;					  	// 1   / sigma^2
				temp_D += log(h[d]->GetBinContent(b)) / temp_s2;						// h_i / sigma^2
				temp_E += h[d]->GetBinCenter(b) * temp_A ;    							// E_i^2 / sigma^2
				temp_F += temp_A; 														// E_i / sigma^2
				temp_G += log(h[d]->GetBinContent(b)) * temp_A;							// h_i*E_i/sigma^2				
			}

		} 

		l_d[d] = ( (temp_B*temp_G) - (temp_D*temp_F) ) / ( (temp_A*temp_F) - (temp_B*temp_E) );

		C_d[d] = ( (temp_G*temp_A) - (temp_D*temp_E) ) / ( (temp_A*temp_F) - (temp_B*temp_E) );
		C_d[d] = exp(C_d[d]) / ( Exposures[d] * h[d]->GetBinWidth(1) ) ;

		


		// if(d == 5)
		// {
		// 	cout << "temp_A = " << temp_A << endl
		// 		 <<	"temp_B = " << temp_B << endl
		// 		 <<	"temp_D = " << temp_D << endl
		// 		 <<	"temp_E = " << temp_E << endl
		// 		 <<	"temp_F = " << temp_F << endl
		// 		 <<	"temp_G = " << temp_G << endl
		// 		 <<	"temp_s2 = "<< temp_s2<< endl
		// 		 <<	"l_d = " 	<< l_d 	  << endl
		// 		 <<	"C_d = " 	<< C_d 	  << endl;

		// }



		// cout << "detector: "<< d  	  << endl	
		// 	 << "l_d[d] = " << l_d[d] << endl	
		// 	 << "C_d[d] = " << C_d[d] << endl
		// 	 << "--------------------"<< endl;
	}

	TF1*     f1[64];
	TCanvas* c1[64];

	gROOT->SetBatch(kTRUE);

	for(int d = 0; d < 64; d++)
	{
		stringstream 	 	   f1_title;				
		f1_title  <<  "f" <<   d+1;
		string f1_tit_str  =   f1_title.str();

		stringstream 	 	   c_title;				
		c_title  <<  "./Figures/Detector_fits/c" <<   d+1 << ".png";
		string c_tit_str  =   c_title.str();

	 	f1[d] = new TF1(f1_tit_str.c_str(), F_Exponential, E_min, E_max, 4  ) ;
	 	c1[d]  = new TCanvas(c_tit_str.c_str(), c_tit_str.c_str(), 1000, 600) ;

	
		f1[d]->SetParNames("N0_1", "Lambda_1","exp_1" ,"Bin_Width_d");
		f1[d]->SetParameters( C_d[d], l_d[d], Exposures[d], h[d]->GetBinWidth(1) );

		h[d]->GetXaxis()->SetRangeUser( E_min, E_max );

		c1[d]->cd();

		h[d]->Draw();
		f1[d]->Draw("SAME");
		c1[d]->SaveAs(c_tit_str.c_str());

	}

	double d_chi2DNF = 0.0;
	double avg_chi   = 0.0;
	double min_chi   = 1e7;
	double max_chi   = 0.0;

	for(int d = 0; d < 64; d++)
	{
		d_chi2DNF = chi2NDF(h[d], E_min, E_max, (C_d[d] *( Exposures[d] * h[d]->GetBinWidth(1) )), l_d[d]);
		avg_chi  += d_chi2DNF / 64;

		cout << "detector:   " << d+1 		<< endl
			 << "d_chi2DNF = " << d_chi2DNF << endl
			 << "-------------------------" << endl;

		if( d_chi2DNF >  max_chi)
		{
			max_chi = d_chi2DNF;
		}
		if( d_chi2DNF < min_chi )
		{
			min_chi = d_chi2DNF;
		}
	}
	cout << "Average chi = " 	<< avg_chi << endl;
	cout << "max chi = " 		<< max_chi << endl;
	cout << "min chi = " 		<< min_chi << endl;


	cout << "---------LAMBDAS----------" << endl;
	for (int d = 0; d < 64; d++)
	{
		cout << l_d[d] << ", ";

	}
	cout << endl;
	cout << "---------C_s ----------" << endl;
	for (int d = 0; d < 64; d++)
	{
		cout << C_d[d] << ", ";

	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// 																//////////////////// 
//////////////////// 																//////////////////// 
//////////////////// 	Fit Each detector by multiplying global fit with exposure 	//////////////////// 
//////////////////// 																//////////////////// 
//////////////////// 																//////////////////// 
//////////////////// 																//////////////////// 

// vector<string> 			root_file_path 	= 		ReadFiles();
// 	TChain* 			 	T_Chain_det		=  		Make_TChain(root_file_path, runs);
// 	vector<CO_detector>     v_det 			= 		Fill_CO_detector(T_Chain_det);

// 	double 					Exposures[64];

// 	double tot = 0;

// 	for(int d = 0; d < 64; d++)
// 	{
// 		Exposures[d] = Get_Exposure(d+1, v_det);
// 		cout << "Exposure of detector " << d+1 << " = " << Exposures[d] << endl;
// 		tot += Exposures[d];
// 	}

// 	TFile* 	tf = new TFile("Final_histograms/1st_cuts_w_flushing/CO_event-TChain_detectors.root");

// 	TH1F* h_d[64];
// 	int detector_n = 1;
// 	for(int d = 0; d < 64; d++)
// 	{
// 		stringstream 			 	   h_title;			//names for histograms		
// 		h_title  <<  "h_" <<   		detector_n;
// 		string h_tit_str  = 	 h_title.str();

// 		h_d[d] = (TH1F*) tf->Get(h_tit_str.c_str());

// 		detector_n++;
// 	}

// 	// TF1* f1[64];
// 	// TCanvas* c1[64];

// 	// for(int d = 0; d < 64; d++)
// 	// {
// 	// 	stringstream 	 	   f1_title;				
// 	// 	f1_title  <<  "f" <<   d+1;
// 	// 	string f1_tit_str  =   f1_title.str();

// 	// 	stringstream 	 	   c_title;				
// 	// 	c_title  <<  "c1" <<   d+1;
// 	// 	string c_tit_str  =   c_title.str();

// 	//  	f1[d] = new TF1(f1_tit_str.c_str(), F_Exponential, E_Fit_min_1, E_Fit_max_1, 5);
// 	//  	c1[d]  = new TCanvas(c_tit_str.c_str(), c_tit_str.c_str(), 1000, 600);

// 	// 	h_d[d]->Rebin(10);

	
// 	// 	f1[d]->SetParNames("N0_1", "Lambda_1", "exp_1", "Bin_Width_1", "Bin_Width_d");
// 	// 	f1[d]->SetParameters( N0_1, lambda_1, Exposures[d], Bin_Width_1, h_d[d]->GetBinWidth(1));

// 	// 	h_d[d]->GetXaxis()->SetRangeUser( E_Fit_min_1, E_Fit_max_1 );

// 	// 	c1[d]->cd();

// 	// 	h_d[d]->Draw();
// 	// 	f1[d]->Draw("SAME");

// 	// }

// 	// TFile* tf2 = new TFile("Final_histograms/Fit_detectors/20210413-1537_Fit_detectors.root","RECREATE");

// 	// TF1* f2[64];
// 	// TCanvas* c2[64];

// 	// for(int d = 0; d < 64; d++)
// 	// {
// 	// 	stringstream 	 	   f1_title;				
// 	// 	f1_title  <<  "f1" <<   d+1;
// 	// 	string f1_tit_str  =   f1_title.str();

// 	// 	stringstream 	 	   c_title;				
// 	// 	c_title  <<  "c2" <<   d+1;
// 	// 	string c_tit_str  =   c_title.str();

// 	//  	f2[d] = new TF1(f1_tit_str.c_str(), F_Exponential, E_Fit_min_2, E_Fit_max_2, 5);
// 	//  	c2[d]  = new TCanvas(c_tit_str.c_str(), c_tit_str.c_str(), 1000, 600);

// 	// 	h_d[d]->Rebin(50);

	
// 	// 	f2[d]->SetParNames("N0_2", "Lambda_2", "exp_2", "Bin_Width_2", "Bin_Width_d");
// 	// 	f2[d]->SetParameters( N0_2, lambda_2, Exposures[d], Bin_Width_2, h_d[d]->GetBinWidth(1));

// 	// 	h_d[d]->GetXaxis()->SetRangeUser( E_Fit_min_2, E_Fit_max_2 );

// 	// 	c2[d]->cd();


// 	// 	h_d[d]->Draw();
// 	// 	f2[d]->Draw("SAME");

// 	// 	f2[d]->Write();
// 	// 	h_d[d]->Write();

// 	cout << "k = " << N0_2 / tot / Bin_Width_2 << endl;
// 	// }
// 	// delete tf2;
// 	TH1F* h_e = new TH1F("h_e", "h_e", 100, 0, 1500);
// 	double events_d[64];

// 	for(int d = 0; d < 64; d++)
// 	{
// 		events_d[d] = 0;
// 		int min_bin_ID = ceil(E_Fit_min_2 / h_d[d]->GetBinWidth(1));
// 		int max_bin_ID = ceil(E_Fit_max_2 / h_d[d]->GetBinWidth(1));

// 		for(int b = min_bin_ID; b <= max_bin_ID; b++)
// 		{
// 			events_d[d] += h_d[d]->GetBinContent(b);
// 		}
// 		// if( events_d[d] / Exposures[d] > )
// 		h_e->Fill(events_d[d] / Exposures[d]);
// 		// cout << "events  " << d +1 << " = " << events_d[d] << endl;
// 		if( (events_d[d] / Exposures[d]) > 50)
// 		{
// 			cout << "events in detector " << d +1 << " = " << events_d[d] / Exposures[d] << endl;

// 		}
// 	}
// 	h_e->Draw();

//////////////////// 																//////////////////// 
//////////////////// 																//////////////////// 
//////////////////// 	Fit Each detector by multiplying global fit with exposure 	//////////////////// 
//////////////////// 																//////////////////// 
//////////////////// 																//////////////////// 
//////////////////// 																//////////////////// 
////////////////////////////////////////////////////////////////////////////////////////////////////////