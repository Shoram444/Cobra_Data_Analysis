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
	string 		year;
	string 		folder;
	string 		file;

};

vector<string>* ReadFiles(string _path, string _key)
{
	vector<string> *f = new vector<string>();
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
		    		f->push_back(dir->d_name);
				}
        }
        closedir(di);
    }
    else {
        perror ("Failed to read dictionaries! Check whether the folders are correctly named!");
        return f;
    }
    return f;
} 


void MP() 
{ 
	vector<string>* years;
	vector<paths>  root_file_path;

	string p;			// request for path in the console. 
    cout<< " p : ";
    cin >> p;

    string key = "20";	// first key, used to save only the folders for years. 20XX. 

	struct paths *r 	= new struct paths; 

	years = ReadFiles(p, key);

	for (int j = 0; j<years->size(); j++)
	{
		vector<string>* folders;
		folders = ReadFiles(p+"/"+years->at(j), "cpg");
		cout << " =================Year: " << years->at(j) <<"=============" << endl;
		cout << endl;

		for(int i = 0; i < folders->size();i++) //files.size()
		{
			cout << " =================folder: " << folders->at(i) << "=============" << endl;
			cout << endl;

			vector<string>* files;
			files = ReadFiles(p + "/" + years->at(j) + "/" + folders->at(i) , ".root");

			for (int k = 0; k < files->size(); k++)
			{
				cout << "File names in folders: " << files->at(k) << endl;

				
				r->year 	= years->at(j);
				r->folder 	= folders->at(i);
				r->file 	= files->at(k);

				cout << "r->year  " << r->year << endl;

				root_file_path.push_back(*r);
				cout<< " root_file_path files" << root_file_path.at(k).year << endl;
			}
		}
	}
	


	// for(int i = 0; i < years->size();i++) //files.size()
	// {
	// 	cout << "File names in files: " << years->at(i) << endl;
	// }

	// char* 	root_file_path;
	// string 	s; 

	// TH1F* h = new TH1F("hh", "haha", 400, 0.0, 4000);

	// for(int i = 0; i < files.size(); i++)
	// {
	// 	s = files.at(i);
	// 	root_file_path = str_to_char(s);

	// 	TFile* f = new TFile(root_file_path);
	// 	TTree* t = (TTree*) f->Get("merged_cal");

	// 	vector<double>* e	= new vector<double>();
	// 	vector<double>* ztc	= new vector<double>();
	// 	vector<double>* aoe	= new vector<double>();
	// 	vector<bool>*	fip	= new vector<bool>();
	// 	vector<bool>*	fbp	= new vector<bool>();

	// 	t->SetBranchAddress("cal_edep", &e);
	// 	t->SetBranchAddress("flag_injected_pulse", &fip);
	// 	t->SetBranchAddress("flag_bad_pulse", &fbp);
	// 	t->SetBranchAddress("cal_ipos_ztc", &ztc);
	// 	t->SetBranchAddress("cal_cpg_diff_AoE", &aoe);

	// 	int n1 	= t->GetBranch("flag_bad_pulse")->GetEntries();
	// 	int n2 	= t->GetBranch("flag_injected_pulse")->GetEntries();
	// 	int n3 	= t->GetBranch("cal_ipos_ztc")->GetEntries();
	// 	int n4 	= t->GetBranch("cal_cpg_diff_AoE")->GetEntries();

	// 	if(!(n1 == n2 && n2 == n3 && n3 == n4)) // Check if the sizes of leaves are the same. The cut function iterates through size of vector, so they must be same. 
	// 	{
	// 		cout << "ERROR: The leaves in this tree are not of the same size!" << endl;
	// 		cout << "Size of flag_bad_pulse = " 		<< n1				   << endl;
	// 		cout << "Size of flag_injected_pulse = " 	<< n2				   << endl; 
	// 		cout << "Size of cal_ipos_ztc = " 			<< n3				   << endl; 
	// 		cout << "Size of cal_cpg_diff_AoE = " 		<< n4				   << endl; 

	// 		return;
	// 	}

	// 	for(int i = 0; i < n1 ; i++)
	// 	{
	// 		t->GetEntry(i);

	// 		for(int j = 0; j < ztc->size(); j++)
	// 		{
	// 			if(   													//apply cuts
	// 				 !(fip->at(j)) 									&&
	// 				 !(fbp->at(j)) 									&&
	// 				  (ztc->at(j) > 0.2 	&& ztc->at(j) < 0.95) 	&&
	// 				  (aoe->at(j) > 0.872 	&& aoe->at(j) < 1.2 )
	// 			   )
	// 			{
	// 				h->Fill(e->at(j));
	// 			}

	// 		}
			
	// 	}
	// }

	// h->Draw();

}