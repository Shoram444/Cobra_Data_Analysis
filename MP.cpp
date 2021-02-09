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

const int n_of_hist =  7;

char* str_to_char(string _str)
{
	int 	n 	= _str.length();
	char* _char = new char[n+1];

    strcpy(_char, _str.c_str());
	
	return _char;
}

struct paths
{
	string 		directory;
	string 		year;
	string 		folder;
	string 		file;
};

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

vector<paths> ReadFiles()
{
    string 		years_key;
	string 	  folders_key;
	string 		files_key;
	string 	   parent_dir;			// request for path in the console.

	vector<paths> files_read;

	vector<string>*   years;
	vector<string>* folders;
	vector<string>*   files;

	struct paths *r  = new struct paths; 

    cout<< "Please specify the path to Data parent directory: ";
    cin >> parent_dir;

    cout<< "Please specify Year of desired data (Use \"20\" for all years): ";
    cin >> years_key;

	cout<< "Please specify keyword for which folders to scan (Use \"cpg\" for all years): ";
	cin >> folders_key;
	
	cout<< "Please specify keyword for which files to scan (Use \".root\" for all years): ";
	cin >> files_key;

	years = ListFiles(parent_dir, years_key);

	for (int j = 0; j<years->size(); j++)
	{
		folders = ListFiles(parent_dir + "/" + years->at(j), folders_key);

		for(int i = 0; i < folders->size();i++) 
		{
			files = ListFiles(parent_dir + "/" + years->at(j) + "/" + folders->at(i) , files_key);

			for (int k = 0; k < files->size(); k++)
			{
				r->directory	=     parent_dir;
				r->year 		=   years->at(j);
				r->folder 		= folders->at(i);
				r->file 		=   files->at(k);

				files_read.push_back(*r);
			}
		}
	}

	return files_read;
}

TH1F* hist_style(string _title)
{
	TH1F* _h = new TH1F(_title.c_str(), _title.c_str(), 1500, 0.0, 3000.0);
	
	_h->SetLineWidth(2);
    _h->GetXaxis()->SetTitle("Energy [keV]");
    _h->GetXaxis()->SetRange(1,1500);
    _h->GetXaxis()->SetLabelFont(42);
    _h->GetXaxis()->SetLabelSize(0.035);
    _h->GetXaxis()->SetTitleSize(0.035);
    _h->GetXaxis()->SetTitleOffset(1);
    _h->GetXaxis()->SetTitleFont(42);
    _h->GetYaxis()->SetTitle("Counts [#]");
    _h->GetYaxis()->SetLabelFont(42);
    _h->GetYaxis()->SetLabelSize(0.035);
    _h->GetYaxis()->SetTitleSize(0.035);
    _h->GetYaxis()->SetTitleFont(42);
    _h->GetZaxis()->SetLabelFont(42);
    _h->GetZaxis()->SetLabelSize(0.035);
    _h->GetZaxis()->SetTitleSize(0.035);
    _h->GetZaxis()->SetTitleOffset(1);
    _h->GetZaxis()->SetTitleFont(42);

    return _h;
}

void MP() 
{ 
	vector<paths>  	 root_file_path; //initialize vector to hold paths for each root file
	root_file_path = 	ReadFiles(); //Fills the vector with patch to each file. The struct holds: parent directory, years, folders, files
	
	TH1F* 		h[n_of_hist];
	TCanvas* 	c[n_of_hist];

	int year_for_title = 2013;
	for(int y = 0; y < n_of_hist; y++)
	{
		stringstream 			 	   h_title;			//names for histograms		
		h_title  <<  "h_" <<   	year_for_title;
		string h_tit_str  = 	 h_title.str();

		stringstream 				   c_title;			//names for canvases
		c_title  <<  "c_" <<    year_for_title;
		string c_tit_str  = 	 c_title.str();	

		h[y] = hist_style(h_tit_str);
		c[y] = new TCanvas(c_tit_str.c_str(), c_tit_str.c_str());

		year_for_title++;
	}	

	for(int i = 0; i < root_file_path.size(); i++)
	{
		char* root_file;

		root_file = str_to_char(root_file_path.at(i).directory + "/" + root_file_path.at(i).year + "/" + root_file_path.at(i).folder + "/" + root_file_path.at(i).file );

		if(i%1000 == 0 ){cout<< i <<" of " << root_file_path.size() <<  " Read" << endl;}

		TFile* f = new TFile(root_file);
		TTree* t = (TTree*) f->Get("merged_cal");

		vector<double>* ene	= new vector<double>();
		vector<double>* ztc	= new vector<double>();
		vector<double>* aoe	= new vector<double>();
		vector<bool>*	fip	= new vector<bool>();
		vector<bool>*	fbp	= new vector<bool>();
		vector<double>*	tim	= new vector<double>();

		t->SetBranchAddress("cal_edep", &ene);
		t->SetBranchAddress("flag_injected_pulse", &fip);
		t->SetBranchAddress("flag_bad_pulse", &fbp);
		t->SetBranchAddress("cal_ipos_ztc", &ztc);
		t->SetBranchAddress("cal_cpg_diff_AoE", &aoe);
		t->SetBranchAddress("info_systime", &tim);

		int n1 	= t->GetBranch("flag_bad_pulse")->GetEntries();
		int n2 	= t->GetBranch("flag_injected_pulse")->GetEntries();
		int n3 	= t->GetBranch("cal_ipos_ztc")->GetEntries();
		int n4 	= t->GetBranch("cal_cpg_diff_AoE")->GetEntries();
		int n5 	= t->GetBranch("info_systime")->GetEntries();

		if(!(n1 == n2 && n2 == n3 && n3 == n4  && n4 == n5)) // Check if the sizes of leaves are the same. The cut function iterates through size of vector, so they must be same. 
		{
			cout << "ERROR: The leaves in this tree are not of the same size!" << endl;
			cout << "Size of flag_bad_pulse = " 		<< n1				   << endl;
			cout << "Size of flag_injected_pulse = " 	<< n2				   << endl; 
			cout << "Size of cal_ipos_ztc = " 			<< n3				   << endl; 
			cout << "Size of cal_cpg_diff_AoE = " 		<< n4				   << endl; 
			cout << "Size of info_systime = " 			<< n5				   << endl; 

			return;
		}

		for(int j = 0; j < n1 ; j++)
		{
			t->GetEntry(j);

			for(int k = 0; k < ztc->size(); k++)
			{
				TTimeStamp* tts = new TTimeStamp(tim->at(k));

				if(   													//apply cuts
					 !(fip->at(k)) 									&&
					 !(fbp->at(k)) 									&&
					  (ztc->at(k) > 0.2 	&& ztc->at(k) < 0.95) 	&&
					  (aoe->at(k) > 0.872 	&& aoe->at(k) < 1.2 )   
				   )					  
				{
					switch(tts->GetDate()/10000) //tts.GetDate() is in format YYYYMMDD 
					{
						case 2013: h[0]->Fill(ene->at(k));
							       break;
						case 2014: h[1]->Fill(ene->at(k));
								   break;
						case 2015: h[2]->Fill(ene->at(k));
								   break;
						case 2016: h[3]->Fill(ene->at(k));
							       break;
						case 2017: h[4]->Fill(ene->at(k));
								   break;
						case 2018: h[5]->Fill(ene->at(k));
								   break;
						case 2019: h[6]->Fill(ene->at(k));
								   break;
					}
				}
				delete tts;
			}
		}
		delete t;
		delete f;
		delete ene;
		delete ztc;
		delete aoe;
		delete fip;
		delete fbp;
	}

	TFile* tf = new TFile("Processed_data.root", "RECREATE");
	for (int d = 0; d < 7; d++)
	{
		// c[d]->cd();
		// h[d]->Draw();
		h[d]->Write();
	}
	delete tf;

}