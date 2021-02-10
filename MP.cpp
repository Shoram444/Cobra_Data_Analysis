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

using namespace boost;
using namespace std;

const int n_of_hist 		=    7;
const bool removeBadPeriods = true;
const bool global 			= true;
const bool layer 			= true;
const bool Det 				= true;

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

vector<vector<double>> read_times_from_rootrc(const char* cutFile, const char* cut )
{
  vector<vector<double>> vTime;
  vector<double> time_pair;
  string line;
  
  ifstream inFile;
  inFile.open(cutFile);
  if(inFile.is_open())
  {
  	cout << "rootrc opened successfully!" << endl;
  }
  else
  {
  	cout << "WARNING: Could not open rootrc file!!" << endl;
  }
  
  // read lines of cutFile
  while(getline(inFile, line)) 
  { 
    if (line.find("#") != std::string::npos){continue;}
    
    // find line with cut specifier, extract start/end time from line and store them in vector "time_pair"  
    if (line.find(TString::Format(".%s.",cut)) != std::string::npos)
    {
      if (line.find("low")!= std::string::npos)
      {
        stringstream buffer(line);
        string str;
        double time;
    
        buffer >> str >> time;

        time_pair.push_back(time);
        // cout << "low " << time << "==>  " <<tsl->GetDate() <<endl;

       }
      
      if (line.find("high")!= std::string::npos)
      {
        stringstream buffer(line);
        string str;
        double time;
    
        buffer >> str >> time;

        time_pair.push_back(time);
        // cout << "high " << time << "==>  " <<tsh->GetDate() << endl;
      }
    }
    
    // store all times of specific cut in one vector ("vTime")
    if (time_pair.size() == 2)
    {  
      // expand time_pair vector with detector number, if line contain "det"  
      if (line.find("det")!=std::string::npos) 
      {
        // extract detector number from string   
        iterator_range<string::iterator> second_dot = find_nth(line, ".", 1);
        iterator_range<string::iterator> third_dot  = find_nth(line, ".", 2);
      
        int start = distance(line.begin(), second_dot.begin());
        int stop  = distance(line.begin(), third_dot.begin() );
      
        string det_string =  line.substr(start+1, stop-start-1);
        stringstream det_no(det_string); 
        int det;
        det_no >> det;
 
        time_pair.push_back(det); 
      }
      
      // expand time_pair vector with first and last detector number of the layer, if line contain "layer" 
      if (line.find("layer")!=std::string::npos)
      {
        // extract layer number from string   
        iterator_range<string::iterator> second_dot = find_nth(line, ".", 1);
        iterator_range<string::iterator> third_dot = find_nth(line, ".", 2);
      
        int start = distance(line.begin(), second_dot.begin());
        int stop = distance(line.begin(), third_dot.begin());
      
        string layer_string =  line.substr(start+1, stop-start-1);
        stringstream layer_no(layer_string); 
        int layer;
        layer_no >> layer;

        time_pair.push_back(16*layer-15);
        time_pair.push_back(16*layer);
      }
      
      vTime.push_back(time_pair);
      time_pair = {};
    }
    
  }
 
  inFile.close();
  return vTime;
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
    _h->GetYaxis()->SetTitle("Counts [#/2keV]");
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
	// int  = 10000;
	vector<bool>* 		b_fip = new vector<bool>();			//true if the pulse was not injected. 
	vector<bool>* 		b_fbp = new vector<bool>();			//true if the pulse was not bad.
	vector<bool>* 		b_ztc = new vector<bool>();			//true if the event passed z criterion.
	vector<bool>* 		b_aoe = new vector<bool>();			//true if the event passed aoe crit.
	vector<bool>* 		b_fsh = new vector<bool>();			//true if the event is outside flushing period.
	vector<double>* 	d_ene = new vector<double>();	
	vector<TTimeStamp>* t_tim = new vector<TTimeStamp>();   //Vector that holds TTimestamp type for each hit



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

	// TH1F* g = new TH1F("g","g", 1500, 0, 3000);
	// TCanvas* cc = new TCanvas("cc", "cc");


	vector<vector<double>> 	 vTime_global; 
	vector<vector<double>> 	  vTime_layer; 
	vector<vector<double>> 		vTime_det; 

	if (global) {vTime_global = read_times_from_rootrc("merged_time_cuts.rootrc", "glob"  ); }
	if (layer ) {vTime_layer  = read_times_from_rootrc("merged_time_cuts.rootrc", "layer" ); }
	if (Det   ) {vTime_det    = read_times_from_rootrc("merged_time_cuts.rootrc", "det"   ); }

	for(int i = 0; i < root_file_path.size(); i++)
	{
		char* root_file;

		root_file = str_to_char(root_file_path.at(i).directory + "/" + 
								root_file_path.at(i).year      + "/" + 
								root_file_path.at(i).folder    + "/" + 
								root_file_path.at(i).file );

		// cout<< root_file << endl;

		if(i%1000 == 0 ){cout<< i <<" of " << root_file_path.size() <<  " Files Read" << endl;}

		TFile* f = new TFile(root_file);
		TTree* t = (TTree*) f->Get("merged_cal");

		vector<double>* ene	= new vector<double>();
		vector<double>* ztc	= new vector<double>();
		vector<double>* aoe	= new vector<double>();
		vector<bool>*	fip	= new vector<bool>();
		vector<bool>*	fbp	= new vector<bool>();
		vector<double>*	tim	= new vector<double>();
		vector<int>*	det = new vector<int>();
		vector<int>*	eid = new vector<int>();

		t->SetBranchAddress("cal_edep", &ene);
		t->SetBranchAddress("cal_ipos_ztc", &ztc);
		t->SetBranchAddress("cal_cpg_diff_AoE", &aoe);
		t->SetBranchAddress("flag_injected_pulse", &fip);
		t->SetBranchAddress("flag_bad_pulse", &fbp);
		t->SetBranchAddress("info_systime", &tim);
		t->SetBranchAddress("cal_det", &det);
		t->SetBranchAddress("info_idx", &eid);

		int n1 	= t->GetBranch("flag_bad_pulse")->GetEntries();
		int n2 	= t->GetBranch("flag_injected_pulse")->GetEntries();
		int n3 	= t->GetBranch("cal_ipos_ztc")->GetEntries();
		int n4 	= t->GetBranch("cal_cpg_diff_AoE")->GetEntries();
		int n5 	= t->GetBranch("info_systime")->GetEntries();
		int n6 	= t->GetBranch("cal_det")->GetEntries();
		int n7 	= t->GetBranch("info_idx")->GetEntries();

		if(!(n1 == n2 && n2 == n3 && n3 == n4  && n4 == n5 && n5 == n6 && n6 == n7)) // Check if the sizes of leaves are the same. The cut function iterates through size of vector, so they must be same. 
		{
			cout << "ERROR: The leaves in this tree are not of the same size!" << endl;
			cout << "Size of flag_bad_pulse = " 		<< n1				   << endl;
			cout << "Size of flag_injected_pulse = " 	<< n2				   << endl; 
			cout << "Size of cal_ipos_ztc = " 			<< n3				   << endl; 
			cout << "Size of cal_cpg_diff_AoE = " 		<< n4				   << endl; 
			cout << "Size of info_systime = " 			<< n5				   << endl; 
			cout << "Size of cal_det = " 				<< n5				   << endl; 
			cout << "Size of info_idx = " 				<< n5				   << endl; 

			return;
		}
		
			
		for(int j = 0; j < n1 ; j++)
		{
			t->GetEntry(j);



			for(int k = 0; k < ene->size(); k++)
			{
				TTimeStamp* tts = new TTimeStamp(tim->at(k));
				t_tim->push_back(*tts);
				d_ene->push_back(ene->at(k));
				b_fsh->push_back(true);
				
				if(fip->at(k))								
				{
					b_fip->push_back(false);
				}
				else
				{
					b_fip->push_back(true);
				}
				

				if(fbp->at(k))								
				{
					b_fbp->push_back(false);
				}
				else
				{
					b_fbp->push_back(true);
				}
				

				if( ztc->at(k) <  0.2  	 || 
					ztc->at(k) >  0.95   ||
					isnan(ztc->at(k))   )	
				{
					b_ztc->push_back(false);
				}
				else
				{
					b_ztc->push_back(true);
				}

				if(aoe->at(k) < 0.872 	||
				   aoe->at(k) > 1.3        )
				{
					b_aoe->push_back(false);
				}
				else
				{
					b_aoe->push_back(true);
				}
				
				for(unsigned int l=0; l<vTime_global.size(); l++)
				{ 
					if( tim->at(k) > vTime_global[l][0] && 
						tim->at(k) < vTime_global[l][1] ) 
					{
						b_fsh->back() = false;
					}
				}

				for(unsigned int l=0; l<vTime_layer.size(); l++)
	            { 
					if( ( tim->at(k) >  vTime_layer[l][0] && 
						  tim->at(k) <  vTime_layer[l][1]    )   &&
					    ( det->at(k) >= vTime_layer[l][2] && 
					      det->at(k) <= vTime_layer[l][3]    )   && 
					      b_fsh->back() == true )
					{
						b_fsh->back() = false;
					} 
	            }
				for(unsigned int l=0; l<vTime_det.size(); l++)
				{ 
					if( (tim->at(k) > vTime_det[l][0]  && 
						 tim->at(k) < vTime_det[l][1]     )   && 
						 det->at(k)== vTime_det[l][2]  && 
						 b_fsh->back() == true  ) 
					{
						b_fsh->back() = false;
					}
				}
				// cout<<"event ID :" 	<<   eid->back()				<< endl;
				// cout<<"det ID :" 	<<   det->back()				<< endl;
				// cout<<"fsh " 		<<   b_fsh->back()		<< endl;
				// cout<<"fip " 		<<   b_fip->back()		<< " " <<fip->at(k)<< endl;
				// cout<<"fbp " 		<<   b_fbp->back()		<< " " <<fbp->at(k)<< endl;
				// cout<<"ztc " 		<<   b_ztc->back()		<< " " <<ztc->at(k)<< endl;
				// cout<<"aoe " 		<<   b_aoe->back()		<< " " <<aoe->at(k)<< endl;


				// if(b_fsh->back() &&
				//    b_fip->back() &&
				//    b_fbp->back() &&
				//    b_ztc->back() &&
				//    b_aoe->back()   )
				// {
				// 	switch(tts->GetDate()/10000) //tts.GetDate() is in format YYYYMMDD 
				// 		{
				// 			case 2013: h[0]->Fill(ene->at(k));
				// 				       break;
				// 			case 2014: h[1]->Fill(ene->at(k));
				// 					   break;
				// 			case 2015: h[2]->Fill(ene->at(k));
				// 					   break;
				// 			case 2016: h[3]->Fill(ene->at(k));
				// 				       break;
				// 			case 2017: h[4]->Fill(ene->at(k));
				// 					   break;
				// 			case 2018: h[5]->Fill(ene->at(k));
				// 					   break;
				// 			case 2019: h[6]->Fill(ene->at(k));
				// 					   break;
				// 		}
				// 	// g->Fill(ene->at(k));
				// }
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
	
	// TFile* tf = new TFile("full_spec.root", "RECREATE");

	// g->Draw();
	// cc->Write();

	// delete tf;
	for(int i = 0; i < d_ene->size(); i++)
	{
		if(b_fsh->at(i) &&
		   b_fip->at(i) &&
		   b_fbp->at(i) &&
		   b_ztc->at(i) &&
		   b_aoe->at(i)   ) //all cuts passed 
		{
			switch(t_tim->at(i).GetDate()/10000) //tts.GetDate() is in format YYYYMMDD 
				{
					case 2013: h[0]->Fill(d_ene->at(i));
						       break;
					case 2014: h[1]->Fill(d_ene->at(i));
							   break;
					case 2015: h[2]->Fill(d_ene->at(i));
							   break;
					case 2016: h[3]->Fill(d_ene->at(i));
						       break;
					case 2017: h[4]->Fill(d_ene->at(i));
							   break;
					case 2018: h[5]->Fill(d_ene->at(i));
							   break;
					case 2019: h[6]->Fill(d_ene->at(i));
							   break;
				}
			// g->Fill(ene->at(k));
		}
	}


	TFile* tf = new TFile("Processed_data_w_flush.root", "RECREATE");
	for (int d = 0; d < 7; d++)
	{
		c[d]->cd();
		h[d]->Draw();
		h[d]->Write();
	}
	delete tf;

}