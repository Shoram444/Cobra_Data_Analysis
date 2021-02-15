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


using namespace boost;
using namespace std;

R__LOAD_LIBRARY(/home/shoram/Work/Diploma_Thesis/CO_event/lib/libCO_event.so);

const int n_of_hist 		=    64;
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
	TH1F* _h = new TH1F(_title.c_str(), _title.c_str(), 5000, 0.0, 10000.0);
	
	_h->SetLineWidth(2);
    _h->GetXaxis()->SetTitle("Energy [keV]");
    _h->GetYaxis()->SetTitle("Counts [ # /2keV]");

    return _h;
}

void MP() 
{ 
	vector<paths>  	 root_file_path; //initialize vector to hold paths for each root file
	root_file_path = 	ReadFiles(); //Fills the vector with patch to each file. The struct holds: parent directory, years, folders, files
	
	TH1F* 		h[n_of_hist];
	TCanvas* 	c[n_of_hist];

	// int year_for_title = 2013;
	int detector_number = 1;
	for(unsigned int y = 0; y < n_of_hist; y++)
	{
		stringstream 			 	   h_title;			//names for histograms		
		h_title  <<  "h_" <<   	detector_number;
		string h_tit_str  = 	 h_title.str();

		stringstream 				   c_title;			//names for canvases
		c_title  <<  "c_" <<    detector_number;
		string c_tit_str  = 	 c_title.str();	

		h[y] = hist_style(h_tit_str);
		c[y] = new TCanvas(c_tit_str.c_str(), c_tit_str.c_str());

		detector_number++;
	}	

	vector<CO_event> v_eve;

	for(unsigned int i = 0; i < root_file_path.size(); i++)
	{
		char* root_file;

		root_file = str_to_char(root_file_path.at(i).directory + "/" + 
								root_file_path.at(i).year      + "/" + 
								root_file_path.at(i).folder    + "/" + 
								root_file_path.at(i).file );

		if(i%10 == 0 ){cout<< i <<" of " << root_file_path.size() <<  " Files Read" << endl;}

		// cout<< root_file << endl;

		TFile* f = new TFile(root_file);
		TTree* t = (TTree*) f->Get("merged_cal");

		vector<double>* 	ene = new vector<double>();
		vector<double>* 	ztc = new vector<double>();
		vector<double>* 	aoe = new vector<double>();
		vector<bool>*		fip = new vector<bool>();
		vector<bool>*		fbp = new vector<bool>();
		vector<double>*		tim = new vector<double>();
		vector<int>*		det = new vector<int>();
		vector<int>*		eid = new vector<int>();

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
		

		for(unsigned int j = 0; j < n1 ; j++) //Extracting cuts and energies from TTree
		{
			t->GetEntry(j);

			for(unsigned int k = 0; k < ene->size(); k++)
			{
				CO_event* e = new CO_event( ene->at(k), ztc->at(k), aoe->at(k), tim->at(k)  ,
							                det->at(k), eid->at(k), fip->at(k), fbp->at(k)	);
				e->InitCuts(0.2 , 0.95 , 0.872 , 1.3 , false , false ); // double _ztc_min , double _ztc_max , double _aoe_min , double _aoe_max ,
								                                        // bool   _fip , bool   _fbp ;
				v_eve.push_back(*e);

				delete	e;
			}
		}
		ene->clear();
		ztc->clear();
		aoe->clear();
		fip->clear();
		fbp->clear();
		delete t;
		delete f;
	}

	for(int i = 0; i < v_eve.size(); i++) //Filling Histograms
	{
		// v_eve.at(i).Print();

		// cout<< "cal_det size : " << v_eve.at(i).Get_c_det() << endl ;

		if( v_eve.at(i).Passed() ) //all cuts passed 
		{
			// v_eve.at(i).Print();
			h[v_eve.at(i).Get_c_det() - 1]->Fill(v_eve.at(i).Get_c_ene());
			// cout<< "cal_det size : " << v_eve.at(i).Get_c_det() << endl ;



			// TTimeStamp* tts = new TTimeStamp(v_eve.at(i).Get_c_tim());

			// switch(tts->GetDate()/10000) //tts.GetDate() is in format YYYYMMDD 
			// {
			// 	case 2013: h[0]->Fill(v_eve.at(i).Get_c_ene());
			// 		       break;
			// 	case 2014: h[1]->Fill(v_eve.at(i).Get_c_ene());
			// 			   break;
			// 	case 2015: h[2]->Fill(v_eve.at(i).Get_c_ene());
			// 			   break;
			// 	case 2016: h[3]->Fill(v_eve.at(i).Get_c_ene());
			// 		       break;
			// 	case 2017: h[4]->Fill(v_eve.at(i).Get_c_ene());
			// 			   break;
			// 	case 2018: h[5]->Fill(v_eve.at(i).Get_c_ene());
			// 			   break;
			// 	case 2019: h[6]->Fill(v_eve.at(i).Get_c_ene());
			// 			   break;
			// }

			// delete tts;
		}
	}

	TFile* tf = new TFile("CO_event-detectors.root", "RECREATE");
	for (int d = 0; d < n_of_hist; d++)
	{
		c[d]->cd();
		h[d]->Draw();
		h[d]->Write();
	}
	delete tf;
}