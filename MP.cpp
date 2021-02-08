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

// vector<string> files;

char* str_to_char(string _str)
{
	int n = _str.length();
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
	struct 	dirent *dir;

	char* 		 c_path; //conversion of string to char. opendir function takes only char. 
	c_path = str_to_char(_path);

	di = opendir(c_path);

	if (di)
    {
        while ((dir = readdir(di)) != NULL)
        {
            stringstream ss;					//backward conversion of char to string.
    		ss  <<  dir->d_name;
    		string strname  = ss.str();
    		if(strname == ".." || strname==".")
    		{
    			continue;
    		}

            std::size_t found = strname.find(_key);			//checking if "key" is present in the string. For Years, the key si 20, for folders - cpg, for root files it's .root. 
		  	if (found!=std::string::npos)
		    	{
		    		file_list->push_back(dir->d_name);
				}
        }
        closedir(di);
    }
    else {
        perror ("Failed to read dictionaries! Check whether the folders are correctly named!");
        return file_list;
    }
    return file_list;
} 

vector<paths> ReadFiles()
{
	string folds;
	string fls;

	vector<paths> files_read;
	string d;			// request for path in the console. 
    cout<< "Please specify the path to Data directory: ";
    cin >> d;

    string key;
    cout<< "Please specify Year of desired data (Use \"20\" for all years): ";
    cin >> key;

	cout<< "Please specify keyword for which folders to scan (Use \"cpg\" for all years): ";
	cin >> folds;
	
	cout<< "Please specify keyword for which files to scan (Use \".root\" for all years): ";
	cin >> fls;



	vector<string>* years;
	vector<string>* folders;
	vector<string>* files;

	struct paths *r 	= new struct paths; 

	years = ListFiles(d, key);

	for (int j = 0; j<years->size(); j++)
	{
		// if(j == 0 )
		// {
		//     cout<< "Please specify keyword for which folders to scan (Use \"cpg\" for all years): ";
  //   		cin >> folds;
		// }

		folders = ListFiles(d + "/" + years->at(j), folds);

		for(int i = 0; i < folders->size();i++) //files.size()
		{
			// cout << years->at(j) << endl;
			// if(i == 0 )
			// {
			//     cout<< "Please specify keyword for which files to scan (Use \".root\" for all years): ";
   //  			cin >> fls;
			// }

			files = ListFiles(d + "/" + years->at(j) + "/" + folders->at(i) , fls);

			for (int k = 0; k < files->size(); k++)
			{
				r->directory	= d;
				r->year 		= years->at(j);
				r->folder 		= folders->at(i);
				r->file 		= files->at(k);

				// cout<< d 					<< "/" << years->at(j) << "/" ;
				// cout<< folders->at(i)	 	<< "/" << files->at(k) << endl;
				files_read.push_back(*r);
			}
		}
	}
	return files_read;
}


void MP() 
{ 
	vector<paths>  root_file_path;
	root_file_path = ReadFiles();
	
	TH1F* h = new TH1F("h", "All Years", 1500, 0.0, 3000.0);
	TH1F* g[6];
	TCanvas* c[6];

	int i = 2013;

	for(int y = 0; y < 7; y++)
	{
		stringstream s_title;					
		s_title  <<  "E_" << i;
		string strname  = s_title.str();

		g[y] = new TH1F(strname.c_str(), strname.c_str(), 1500, 0.0, 3000.0);

	   g[y]->SetLineWidth(2);
	   g[y]->GetXaxis()->SetTitle("Energy [keV]");
	   g[y]->GetXaxis()->SetRange(1,1500);
	   g[y]->GetXaxis()->SetLabelFont(42);
	   g[y]->GetXaxis()->SetLabelSize(0.035);
	   g[y]->GetXaxis()->SetTitleSize(0.035);
	   g[y]->GetXaxis()->SetTitleOffset(1);
	   g[y]->GetXaxis()->SetTitleFont(42);
	   g[y]->GetYaxis()->SetTitle("Counts [#]");
	   g[y]->GetYaxis()->SetLabelFont(42);
	   g[y]->GetYaxis()->SetLabelSize(0.035);
	   g[y]->GetYaxis()->SetTitleSize(0.035);
	   g[y]->GetYaxis()->SetTitleFont(42);
	   g[y]->GetZaxis()->SetLabelFont(42);
	   g[y]->GetZaxis()->SetLabelSize(0.035);
	   g[y]->GetZaxis()->SetTitleSize(0.035);
	   g[y]->GetZaxis()->SetTitleOffset(1);
	   g[y]->GetZaxis()->SetTitleFont(42);

		i++;
	}	

	
	for(int i = 0; i < root_file_path.size(); i++)
	{
		char* root_file;

		root_file = str_to_char(root_file_path.at(i).directory + "/" + root_file_path.at(i).year + "/" + root_file_path.at(i).folder + "/" + root_file_path.at(i).file );

		if(i%1000 == 0 ){cout<< i <<" of " << root_file_path.size() <<  " Read" << endl;}

		// cout<< root_file << endl;

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
				// cout<< tts->GetDate()/10000 << endl;
				// TTimeStamp* maxTime = new TTimeStamp();
				// TTimeStamp* minTime = new TTimeStamp();
				if(   													//apply cuts
					 !(fip->at(k)) 									&&
					 !(fbp->at(k)) 									&&
					  (ztc->at(k) > 0.2 	&& ztc->at(k) < 0.95) 	&&
					  (aoe->at(k) > 0.872 	&& aoe->at(k) < 1.2 )   
				   )					  // (tts->AsDouble() < maxTime->AsDouble())		&&  // (tts->AsDouble() > minTime->AsDouble())
				{
					
					switch(tts->GetDate()/10000)
					{
						case 2013: g[0]->Fill(ene->at(k));
							       break;
						case 2014: g[1]->Fill(ene->at(k));
								   break;
						case 2015:
						{
							g[2]->Fill(ene->at(k));
							break;
						}		
						case 2016:
						{
							g[3]->Fill(ene->at(k));
							break;
						}		
						case 2017:
						{
							g[4]->Fill(ene->at(k));
							break;
						}		
						case 2018:
						{
							g[5]->Fill(ene->at(k));
							break;
						}	
						case 2019:
						{
							g[6]->Fill(ene->at(k));
							break;
						}			
					}
					// cout<< "Time of event: " << tim->at(k) << " ======== "  ;
					// tts->Print() ;
					h->Fill(ene->at(k));
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

	// h->Draw();

	TFile* year_hist = new TFile("histogram.root", "RECREATE");
	// g[0]->Draw();

	g[0]->Write();
	g[1]->Write();
	g[2]->Write();
	g[3]->Write();
	g[4]->Write();
	g[5]->Write();
	g[6]->Write();

}