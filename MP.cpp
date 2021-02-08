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
	vector<paths> files_read;
	string d;			// request for path in the console. 
    cout<< "Please specify the path to Data directory: ";
    cin >> d;

    string key;
    cout<< "Please specify Year of desired data (Use \"20\" for all years): ";
    cin >> key;

	vector<string>* years;
	vector<string>* folders;
	vector<string>* files;

	struct paths *r 	= new struct paths; 

	years = ListFiles(d, key);

	for (int j = 0; j<years->size(); j++)
	{
		if(j == 0 )
		{
		    cout<< "Please specify keyword for which folders to scan (Use \"cpg\" for all years): ";
    		cin >> key;
		}

		folders = ListFiles(d+"/"+years->at(j), key);

		for(int i = 0; i < folders->size();i++) //files.size()
		{
			if(i == 0 )
			{
			    cout<< "Please specify keyword for which files to scan (Use \".root\" for all years): ";
    			cin >> key;
			}

			files = ListFiles(d + "/" + years->at(j) + "/" + folders->at(i) , key);

			for (int k = 0; k < files->size(); k++)
			{
				r->directory	= d;

				r->year 		= years->at(j);
				r->folder 		= folders->at(i);
				r->file 		= files->at(k);

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
	
	TH1F* h = new TH1F("h", "hh", 1500, 0.0, 3000.0);
	// TH1F* h[6];
	// for(int y = 0; y < 7; y++)
	// {
	// 	char* c_year;
	// 	c_year = str_to_char(root_file_path.at(0).year);

	// 	h[y] = new TH1F("h", c_year, 200, 0.0, 2000.0);
	// }	

	cout << "========== Files Read! ==============" <<endl;
	

	for(int i = 0; i < root_file_path.size(); i++)
	{
		char* root_file;

		root_file = str_to_char(root_file_path.at(i).directory + "/" + root_file_path.at(i).year + "/" + root_file_path.at(i).folder + "/" + root_file_path.at(i).file );

		cout<< root_file << endl;

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

		for(int i = 0; i < n1 ; i++)
		{
			t->GetEntry(i);

			for(int j = 0; j < ztc->size(); j++)
			{
				TTimeStamp* tts = new TTimeStamp(tim->at(j));
				TTimeStamp* maxTime = new TTimeStamp( 2013 , 11 , 25 , 8 , 35 , 0 , 0);
				TTimeStamp* minTime = new TTimeStamp( 2013 , 11 , 25 , 5 , 35 , 0 , 0);

				if(   													//apply cuts
					 !(fip->at(j)) 									&&
					 !(fbp->at(j)) 									&&
					  (ztc->at(j) > 0.2 	&& ztc->at(j) < 0.95) 	&&
					  (aoe->at(j) > 0.872 	&& aoe->at(j) < 1.2 )   &&
					  (tts->AsDouble() < maxTime->AsDouble())		&&
					  (tts->AsDouble() > minTime->AsDouble())
				   )
				{


					cout<< "Time of event: " << tim->at(j) << " ======== "  ;
					tts->Print() ;
					h->Fill(ene->at(j));
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
		
		h->Draw();

	}

<<<<<<< HEAD
=======
	h->Draw();
>>>>>>> Develop
	// TFile* year_hist = new TFile("histogram_years.root", "NEW");
	// h->Write();

}