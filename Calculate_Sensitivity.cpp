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

#include "/home/shoram/Work/Diploma_Thesis/MPFC/include/MPFeldman_Cousins.hh"
#include "/home/shoram/Work/Diploma_Thesis/CO_detector/include/CO_detector.hh"
#include "/home/shoram/Work/Diploma_Thesis/CO_event/include/CO_event.hh"

using namespace boost;
using namespace std;


R__LOAD_LIBRARY(/home/shoram/Work/Diploma_Thesis/MPFC/lib/libMPFC.so);
R__LOAD_LIBRARY(/home/shoram/Work/Diploma_Thesis/CO_detector/lib/libCO_detector.so);
R__LOAD_LIBRARY(/home/shoram/Work/Diploma_Thesis/CO_event/lib/libCO_event.so);


const double ln2 	= log(2);
const double N_A 	= 6.02e23;
const double W [9]	= { 107.90418, 113.903365, 127.904461, 69.92532,				//molar masses of isotopes: 108Cd, 114Cd, 128Te, 
						63.929142, 119.90406 , 129.906223, 105.90646, 115.904763} ; // 				70Zn, 64Zn, 120Te, 130Te, 106Cd, 116Cd
																					
const double a [9]  = { 0.89, 28.73, 31.69, 0.61, 49.17, 0.10, 33.8, 1.25, 7.5};	//natural abundances of isotopes: 108Cd, 114Cd, 128Te, 
																					// 				70Zn, 64Zn, 120Te, 130Te, 106Cd, 116Cd
const double Q [9]  = { 271.8, 542.5, 866.5, 997.1, 1094.7, 						//Q of isotopes: 108Cd, 114Cd, 128Te, 
						1730.4, 2527.5, 2775.4, 2813.4};							// 				70Zn, 64Zn, 120Te, 130Te, 106Cd, 116Cd




///////HEADER Functions /////////

struct data_partition
{
	int 			p_dno;
	Double_t 		p_sta_tim;
	Double_t 		p_end_tim;
	vector<double> 	p_ene;
	double 			p_exp;
	double 			p_p0;
	double 			p_p1;
	double 			p_p2;


	double Get_p_res(double _E)
	{
		double 		p_res 			=  0; 
	    int         power           =  2;
	    double      p1_times_sqrt_E = p_p1*sqrt(_E);  
	    double      p2_times_E      = p_p2*_E;  

	    p_res      = sqrt( pow(p_p0,power) + pow(p1_times_sqrt_E,power) + pow(p2_times_E, power) );
	    return p_res;
	}
	
};

double 					Get_Resolution(int _d_id, double _E, vector<CO_detector>  _v_det); 
vector<CO_detector>  	Fill_CO_detector(TChain* _T_Chain);
vector<data_partition>  create_partitions(TChain* _T_Chain_det, TChain* _T_Chain_eve , int _dno);
TH1F* hist_style(string _title);


///////FUNCTIONS - from MP.cpp ////////////
const int n_of_hist 		=   64;
const int n_of_det			= 	65;   ///WATCH OUT! Detector numbering starts from 1, so there has to be n+1 n_of_det
const bool removeBadPeriods = true;
const bool global 			= true;
const bool layer 			= true;
const bool Det 				= true;
const char* 		   runs = "runs";
const char*      merged_cal = "merged_cal";
// TChain* T_Chain;

struct calib_date
{
	int year;
	int month;
	int day;
	int hour = 0;
	int minute = 0;
	int second = 0;
	int nanosecond = 0;
};

TChain* 					Make_TChain(vector<string>  	 root_file_path, const char* _ttree);
char* 						str_to_char(string _str);
vector<string>* 			ListFiles(string _path, string _key);
vector<vector<double>> 		read_times_from_rootrc(const char* cutFile, const char* cut );
vector<string> 				ReadFiles(bool _defaultPath = true);

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


void sanity_check(TChain* tc)
{
	int n1 	= tc->GetBranch("flag_bad_pulse")->GetEntries();
	int n2 	= tc->GetBranch("flag_injected_pulse")->GetEntries();
	int n3 	= tc->GetBranch("cal_ipos_ztc")->GetEntries();
	int n4 	= tc->GetBranch("cal_cpg_diff_AoE")->GetEntries();
	int n5 	= tc->GetBranch("info_systime")->GetEntries();
	int n6 	= tc->GetBranch("cal_det")->GetEntries();
	int n7 	= tc->GetBranch("info_idx")->GetEntries();

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
}


double Get_Total_Exposure(int _d_id, vector<CO_detector>   _v_det, bool _all = false)
{
	CO_detector 	a_det[n_of_det];

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

				a_det[det_Id].Add( det_Id, c_mas , c_dur);
			}
		}
		else
		{
			int 		c_dno = _v_det.at(i).get_c_dno();
			double 		c_mas = _v_det.at(i).get_c_mas();
			Double_t	c_dur = _v_det.at(i).get_c_dur();

			a_det[c_dno].Add( c_dno, c_mas, c_dur);
		}
	}

	if(!_all)
	{
  		total_exposure = a_det[det_Id].calc_exposure();
	}
	else
	{
  		for( int d = 1; d < n_of_det; d++)			///WATCH OUT! Detector numbering starts from 1, so there has to be n+1 n_of_det
		{
			total_exposure += a_det[d].calc_exposure();
		}
	}

  	return total_exposure;

}
///////FUNCTIONS - from MP.cpp ////////////

double Get_Resolution(int _d_id, double _E, vector<CO_detector>  _v_det)
{
	double d_res;

	for(int i = 0; i < _v_det.size(); i++)
	{
		for(int j = 0; j<64 ; j++)
		{
			if(_v_det.at(i).get_c_dno() == _d_id)
			{
				d_res = _v_det.at(i).calc_sensitivity(_E);
			}
		}
	}

  	return d_res;

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
			CO_detector* d = new CO_detector( dno->at(j), mas->at(j), dur, tim, p0->at(j), p1->at(j), p2->at(j) );
			v_det.push_back(*d);
			delete	d;
		}

		dno->clear();
		mas->clear();
		p0->clear();
		p1->clear();
		p2->clear();
	}

	return v_det;

}

//////////////////////////READING THE DATES FOR CALIBRATION PERIODS //////////


vector<string> lines_from_file(string fileName)
{
  string strLine = "";
  vector<string> vLine;

  ifstream infile;
  infile.open(fileName.c_str());
    
  while( infile.good() )
  {  
    if(infile.is_open())
    {
      infile >> strLine;
      if(strLine=="") continue;
      vLine.push_back(strLine);
    }
  }
  
  infile.close();
  
  return vLine;
}

vector<calib_date> line_split(vector<string> _dates)
{
	vector<calib_date> v_cal_date;
	calib_date cd;
	for(int i = 0; i < _dates.size()-1; i++)
	{
		size_t found = _dates.at(i).find("-");
		if(found == string::npos)
			{	
				continue;
			}
		if(found != string::npos)
		{
			string date = _dates.at(i);
	        int _year = stoi(date.substr(0, found));
	        int _month = stoi(date.substr(found+1, 2));
	        int _day = stoi(date.substr(found+4, 2));

	        cd.year = _year;
	        cd.month = _month;
	        cd.day = _day;

	        v_cal_date.push_back(cd);
		}
	}

	return v_cal_date;
}


vector<data_partition> create_partitions(TChain* _T_Chain_det, TChain* _T_Chain_eve , int _dno)
{
	vector<data_partition> v_part;
	int    	  sec_per_day = 86400;	//this is used in calculating exposure in units kgd
	int 	  k 		  = 0;  	//partition indexing



	vector<int>* 	dno = new vector<int>;
	vector<double>* mas = new vector<double>;
	Double_t		dur;
	Double_t		tim;
	vector<double>* p0  = new vector<double>;
	vector<double>* p1  = new vector<double>;
	vector<double>* p2  = new vector<double>;


	_T_Chain_det->SetBranchAddress("det_no", &dno);
	_T_Chain_det->SetBranchAddress("det_mass", &mas);
	_T_Chain_det->SetBranchAddress("duration", &dur);
	_T_Chain_det->SetBranchAddress("startTime", &tim);
	_T_Chain_det->SetBranchAddress("det_res_p0", &p0);
	_T_Chain_det->SetBranchAddress("det_res_p1", &p1);
	_T_Chain_det->SetBranchAddress("det_res_p2", &p2);

	data_partition* 	   p_part = new data_partition();

	p_part->p_dno = _dno;
	p_part->p_p0  =  -1.0;
	p_part->p_p1  =  -1.0;
	p_part->p_p2  =  -1.0;
	p_part->p_exp =   0  ;
	
	v_part.push_back(*p_part);


	for(unsigned int i = 0; i < _T_Chain_det->GetEntries(); i++ ) //
	{
		// if(i%10000==0) cout << i << " of " << _T_Chain_det->GetEntries() << " Read!" << endl;
		_T_Chain_det->GetEntry(i);

		for( unsigned int j = 0 ; j < dno->size() ; j++ )
		{
			data_partition* 	   temp_part = new data_partition();

			if(isnan(p0->at(j))){continue;} //there are entries where p0 is nan. These will be skipped over. 

			if(dno->at(j) == _dno)
			{
				if( k 				   == 0 		&& 			//for the first partition (initiated before) fill it up with the first entry in tchain
				    v_part.at(k).p_p0  == -1.0		&&
				    v_part.at(k).p_p1  == -1.0		&&
				    v_part.at(k).p_p2  == -1.0		
				  )
				{
					v_part.at(k).p_dno 		= dno->at(j);
					v_part.at(k).p_sta_tim 	= tim;
					v_part.at(k).p_exp 		= dur * mas->at(j) / sec_per_day;
					v_part.at(k).p_p0 		= p0->at(j);
					v_part.at(k).p_p1 		= p1->at(j);
					v_part.at(k).p_p2 		= p2->at(j);


					// stringstream 			 	   h_title;			//names for histograms		
					// h_title  <<  "h_" <<   	_dno << "-" << k;
					// string h_tit_str  = 	 h_title.str();


					// v_part.push_back(*temp_part);
				}

				else if( 	v_part.at(k).p_p0 == p0->at(j) 		&& //if all the parameters are the same in a partition, add exposure and end time as start of the partition + its duration
							v_part.at(k).p_p1 == p1->at(j)     	&&
							v_part.at(k).p_p2 == p2->at(j)
					   )
				{
					v_part.at(k).p_exp		+= dur * mas->at(j) / sec_per_day;
					v_part.at(k).p_end_tim 	= tim + dur; 			
				}

				else if(	v_part.at(k).p_p0 != p0->at(j)     || // if either of the parameters has changed, create new partition
							v_part.at(k).p_p1 != p1->at(j)	   ||
							v_part.at(k).p_p2 != p2->at(j)
					   )
				{
					k+=1;

					temp_part->p_dno 		= dno->at(j);
					temp_part->p_sta_tim 	= tim;
					temp_part->p_exp 		= dur * mas->at(j) / sec_per_day;
					temp_part->p_p0 		= p0->at(j);
					temp_part->p_p1 		= p1->at(j);
					temp_part->p_p2 		= p2->at(j);

					// stringstream 			 	   h_title;			//names for histograms		
					// h_title  <<  "h_" <<   	_dno << "-" << k;
					// string h_tit_str  = 	 h_title.str();


					v_part.push_back(*temp_part);
				}
			}
			delete temp_part;
		}
	}

	// vector<float>* 		ene = new vector<float>();
	// vector<float>* 		ztc = new vector<float>();
	// vector<float>* 		aoe = new vector<float>();
	// vector<bool>*		fip = new vector<bool>();
	// vector<bool>*		fbp = new vector<bool>();
	// vector<double>*		eti = new vector<double>(); //time of the event had to be renamed to eti, cause tim is already present in detector tree
	// vector<int>*		det = new vector<int>();
	// vector<int>*		eid = new vector<int>();

	// _T_Chain_eve->SetBranchAddress("cal_edep", &ene);
	// _T_Chain_eve->SetBranchAddress("cal_ipos_ztc", &ztc);
	// _T_Chain_eve->SetBranchAddress("cal_cpg_diff_AoE", &aoe);
	// _T_Chain_eve->SetBranchAddress("flag_injected_pulse", &fip);
	// _T_Chain_eve->SetBranchAddress("flag_bad_pulse", &fbp);
	// _T_Chain_eve->SetBranchAddress("info_systime", &eti);
	// _T_Chain_eve->SetBranchAddress("cal_det", &det);
	// _T_Chain_eve->SetBranchAddress("info_idx", &eid);

	// sanity_check(_T_Chain_eve);

	// int en = 0;

	// for(unsigned int j = 0; j < _T_Chain_eve->GetEntries(); j++)
	// {
	// 	// if(j%10000==0) cout << j << " of " << _T_Chain_eve->GetEntries() << " Read!" << endl;
	// 	_T_Chain_eve->GetEntry(j);

	// 	for(unsigned int k = 0; k < ene->size(); k++) ///MAROS EVENT
	// 	{
	// 		CO_event* e = new CO_event( ene->at(k), ztc->at(k), aoe->at(k), eti->at(k)  ,
	// 					                det->at(k), eid->at(k), fip->at(k), fbp->at(k)	);
	// 		e->InitCuts(0.2 , 0.95 , 0.872 , 1.3 , false , false ); // double _ztc_min , double _ztc_max , double _aoe_min , double _aoe_max ,
	// 						                                        // bool   _fip , bool   _fbp ;

	// 		if(e->Passed())
	// 		{
	// 			for(int l = 0; l < v_part.size() ; l++)
	// 			{
	// 				if(e->Get_c_tim() > v_part.at(l).p_sta_tim &&
	// 				   e->Get_c_tim() < v_part.at(l).p_end_tim )
	// 				{
	// 					v_part.at(l).p_ene.push_back(e->Get_c_ene());
	// 					// cout << " energy in partition: " << v_part.at(l).p_ene.at(en);
	// 					// en++;
	// 				}
	// 			}
	// 		}
	// 		// v_eve.push_back(*e);

	// 		delete	e;
	// 	}

	// 	ene->clear();
	// 	ztc->clear();
	// 	aoe->clear();
	// 	fip->clear();
	// 	fbp->clear();
	// 	eti->clear();
	// 	det->clear();
	// 	eid->clear();
	// }



	return v_part;
}


TH1F* hist_style(string _title)
{
	TH1F* _h = new TH1F(_title.c_str(), _title.c_str(), 5000, 0.0, 10000.0);
	
	_h->SetLineWidth(2);
    _h->GetXaxis()->SetTitle("Energy [keV]");
    _h->GetYaxis()->SetTitle("Counts [ # /2keV]");

    return _h;
}


void Calculate_Sensitivity() 
{
	int 	detector_ID  = 	  1;
	int 	gaus_cutoff  =  320;

	double 	entries 		[9];
	double	sensitivities 	[9];
	double  t_half 			[9];
	int 	min_bin_ID 		[9];
	int 	max_bin_ID 		[9];
	int 	Q_bin_ID		[9];
	double 	resolution 		[9];

	vector<string> 			root_file_path 	= 		ReadFiles();
	TChain* 			 	T_Chain_det		=  		Make_TChain(root_file_path, runs);
	TChain* 			 	T_Chain_eve	    =  		Make_TChain(root_file_path, merged_cal);

	//// CREATION OF PARTITIONS 					 	/////////////////////

	int total_partitions = 0;
	double exposure_partitions = 0;
	vector<data_partition>       v_part[64];

	cout << " =====================================" << endl;
	cout << " dno \t | start \t | end \t | exp \t | p0 \t | p1 \t | p2 \t | " << endl;
	for(int i = 0; i < 64; i++)
	{
		v_part[i]	= 		create_partitions(T_Chain_det, T_Chain_eve, i+1);
		for(int j = 0 ; j < v_part[i].size(); j++)
		{
			cout 	<< v_part[i].at(j).p_dno 		<< " \t|" 
					<< v_part[i].at(j).p_sta_tim 	<< " \t|"
					<< v_part[i].at(j).p_end_tim 	<< " \t|"
					<< v_part[i].at(j).p_exp	 	<< " \t|"
					<< v_part[i].at(j).p_p0	 		<< " \t|"
					<< v_part[i].at(j).p_p1	 		<< " \t|"
					<< v_part[i].at(j).p_p2	 		<< " \t|" << endl;
			exposure_partitions += v_part[i].at(j).p_exp;
		}
		total_partitions += v_part[i].size();
		cout << " =====================================" << endl;

	}
	cout << " Number of partitions created : " 	<< total_partitions << endl;
	cout << " Total Exposure : " 				<< exposure_partitions << endl;

	TH1F* 		h_p[total_partitions];
	TCanvas* 	c_p[total_partitions];
	
	ofstream myfile;
	myfile.open ("Data_Partitions/Detector_Partitions_info.txt");
	myfile << "detector number; partition number; start time; end time; no of events inside; exposure; p0;p1;p2 \n";
	for(unsigned int i = 0; i < 64; i++)
	{

		for(int p = 0; p < v_part[i].size(); p++)
		{
			stringstream 			 	   h_title;			//names for histograms		
			h_title  <<  "h_" <<   	i << "-" << p;
			string h_tit_str  = 	 h_title.str();

			stringstream 				   c_title;			//names for canvases
			c_title  <<  "c_" <<    i << "-" << p;
			string c_tit_str  = 	 c_title.str();	

			h_p[i] = hist_style(h_tit_str);
			c_p[i] = new TCanvas(c_tit_str.c_str(), c_tit_str.c_str());

			myfile  << v_part[i].at(p).p_dno 		<< "; "
					<< p 							<< "; "
					<< v_part[i].at(p).p_sta_tim	<< "; "
					<< v_part[i].at(p).p_end_tim	<< "; "
					<< v_part[i].at(p).p_ene.size()	<< "; "
					<< v_part[i].at(p).p_exp		<< "; "
					<< v_part[i].at(p).p_p0 		<< "; "
					<< v_part[i].at(p).p_p1 		<< "; "
					<< v_part[i].at(p).p_p2 		<< "\n";

			// for(int j = 0; j < v_part[i].at(p).p_ene.size(); j++)
			// {
			// 	h_p[i]->Fill(v_part[i].at(p).p_ene.at(j));

			// 	// stringstream 	partition_info;
				
			// 	// myfile << partition_info ;
			// }
		}
	}	
	myfile.close();



	// TFile* tf = new TFile("CALCULATE_SENSITIVITY_ALL_HISTOGRAMS.root", "RECREATE");
	// for (int i = 0; i < total_partitions; i++)
	// {
	// 	// c_p[i]->cd();
	// 	// h_p[i]->Draw();
	// 	h_p[i]->Write();
	// }
	// delete tf;

	//// CREATION OF PARTITIONS 					 	/////////////////////

	//// ----------------------------------------- 		/////////////////////

	//// CALCULATING T1/2 FROM HISTOGRAMS	 		 	/////////////////////


	// MPFeldman_Cousins* 		obj 			= 		new MPFeldman_Cousins(gaus_cutoff, 0.9);
	// TFile* 					tf 				= 		new TFile("./Final_histograms/1st_cuts_w_flushing/CO_event-TChain_detectors.root");

	// for(int n = 1; n <= 64; n++)
	// {
	// 	stringstream 			 	   h_title;
	// 	h_title  <<  "h_" <<   		detector_ID;
	// 	string h_tit_str  = 	 h_title.str();

	// 	TH1F* 	h = (TH1F*) tf->Get(h_tit_str.c_str()); // getting all detectors

	// 	for(int i =0; i < 9; i++) 		//Determining ROI. Q_bin_ID is the bin number where Q lies. min_bin_ID is the bin number of -3 sigma from Q. 
	// 	{
	// 		Q_bin_ID[i] 	= ceil(Q[i]/2);  //the assumption is that the binning is 2kev/bin (this might change!)
	// 		resolution[i] 	= Get_Resolution(detector_ID, Q[i], v_det);

	// 		min_bin_ID[i]	= ceil ((Q[i] - 3 * resolution[i])/ 2);  	
	// 		max_bin_ID[i]   = ceil ((Q[i] + 3 * resolution[i] )/ 2);	
	// 	}


	// 	for(int j = 0; j <9; j++)
	// 	{
	// 		entries[j] = 0;
			
	// 		for(int i = min_bin_ID[j]; i <= max_bin_ID[j]; i++)
	// 		{
	// 			if(i == min_bin_ID[j]) //proportion of the amount of entries in relation to position of -3sigma within bin
	// 			{
	// 				double ratio_factor = (min_bin_ID[j] - (Q[j] - 3 * resolution[j] ) / 2) / 2;  
	// 				double entries_full = h->GetBinContent(i);

	// 				entries[j] += ratio_factor * entries_full;
	// 			}
	// 			else if(i == max_bin_ID[j]) //proportion of the amount of entries in relation to position of +3sigma within bin
	// 			{
	// 				double ratio_factor = (max_bin_ID[j] - (Q[j] - 3 * resolution[j] ) / 2) / 2; 
	// 				double entries_full = h->GetBinContent(i);

	// 				entries[j] += ratio_factor * entries_full;
	// 			}
	// 			else
	// 			{
	// 				entries[j] += h->GetBinContent(i); 				//sum of all events in ROI
	// 			}
	// 		}
	// 		if(entries[j] > gaus_cutoff)
	// 		{
	// 			sensitivities[j] = 3*sqrt(entries[j]);			//if there are more events in ROI than 320, gaussian distribution is used for sensitivity, else FC
	// 		}
	// 		else
	// 		{
	// 			sensitivities[j] = obj->get_sensitivity(entries[j]);
	// 		}

	// 	}
	// 	cout << " =========== Detector number =  " << detector_ID << " ==============" << endl;

	// 	double exposure 	= 	Get_Total_Exposure(detector_ID, v_det);

	// 	for( int i = 0 ; i < 9 ; i++)
	// 	{
	// 		t_half[i] = ln2 * N_A * a[i] * exposure / ( W[i] * sensitivities[i] );

	// 		cout << "T_half : " << t_half[i] << endl; 

	// 	}

	// 	cout << " ===================================================================" << endl;

	// 	detector_ID += 1;

	// }

	//// CALCULATING T1/2 FROM HISTOGRAMS	 		 	/////////////////////

	//// ----------------------------------------- 		/////////////////////

 	//////////////////////////////////////// Detector Parameters ////////////////////	

	// vector<string> 	    dates 	= lines_from_file("Calib_Dates.txt");
	// vector<calib_date>  v_cd    = line_split(dates);

	// Double_t 			d_det_res_p0;
	// Double_t 			d_det_res_p1;
	// Double_t 			d_det_res_p2;
	// Double_t			d_time;

	// TCanvas* 	c0 = new TCanvas("c0","p0 parameter in time",900,600);
	// TCanvas* 	c1 = new TCanvas("c1","p1 parameter in time",900,600);
	// TCanvas* 	c2 = new TCanvas("c2","p2 parameter in time",900,600);
	// TGraph* 	tg0 = new TGraph();//points, d_time, d_det_res_p0);
	// TGraph* 	tg1 = new TGraph();//points, d_time, d_det_res_p0);
	// TGraph* 	tg2 = new TGraph();//points, d_time, d_det_res_p0);


	// int points = 0;
	// for(int i = 0; i < v_det.size(); i++)
	// {

	// 	if(v_det.at(i).get_c_dno() == 1)
	// 	{
	// 		d_det_res_p0 = (Double_t) v_det.at(i).get_c_p0();
	// 		d_det_res_p1 = (Double_t) 10*v_det.at(i).get_c_p1();
	// 		d_det_res_p2 = (Double_t) 1000*v_det.at(i).get_c_p2();
	// 		d_time = v_det.at(i).get_c_tim();

	// 		tg0->SetPoint(points,d_time, d_det_res_p0);
	// 		tg1->SetPoint(points,d_time, d_det_res_p1);
	// 		tg2->SetPoint(points,d_time, d_det_res_p2);
	// 		points +=1;

	// 	}

	// }

 // 	tg0->SetMarkerStyle(8);
 // 	tg0->SetMarkerColor(3);

	// tg1->SetMarkerStyle(2);
 // 	tg1->SetMarkerColor(4);

	// tg2->SetMarkerStyle(22);
 // 	tg2->SetMarkerColor(5);



	// tg0->SetTitle("Resolution Fit Parameter 0 in Time, Detector ID: 1");
	// tg0->GetXaxis()->SetTitle("Unix Time [ns]");
	// tg0->GetYaxis()->SetRangeUser(0, 55);
	// tg0->GetYaxis()->SetTitle("Parameter 0 [keV]");

	

	// // TLine* tl0[i] = new TLine(	1389271413, 0 ,1389271413 ,55 )	;


	// tg1->SetTitle("Resolution Fit Parameter 1 in Time, Detector ID: 1");
	// tg1->GetXaxis()->SetTitle("Unix Time [ns]");
	// tg1->GetYaxis()->SetRangeUser(0, 1);
	// tg1->GetYaxis()->SetTitle("Parameter 1 [keV]");	

	// tg2->SetTitle("Resolution Fit Parameter 2 in Time, Detector ID: 1");
	// tg2->GetXaxis()->SetTitle("Unix Time [ns]");
	// tg2->GetYaxis()->SetRangeUser(-0.01, 0.015);
	// tg2->GetYaxis()->SetTitle("Parameter 2 [keV]");

	// c0->cd();
	// tg0->Draw("ap");
	// TLine* tl0[v_cd.size()];
	// TLine* tl1[v_cd.size()];
	// TLine* tl2[v_cd.size()];
	// for(int i = 0; i < v_cd.size() - 1 ; i++ )
	// {
	// 	TTimeStamp* tts = new TTimeStamp(v_cd.at(i).year,v_cd.at(i).month, v_cd.at(i).day, v_cd.at(i).hour, v_cd.at(i).minute, v_cd.at(i).second, v_cd.at(i).nanosecond );
	// 	tl0[i] = new TLine(	tts->AsDouble(), 0 ,tts->AsDouble() ,55 )	;
	// 	tl0[i]->Draw("SAME");
	// }
	// // tg0->Draw("al");
	// // tl0->Draw();

	// c1->cd();
	// tg1->Draw("ap");
	// for(int i = 0; i < v_cd.size() - 1 ; i++ )
	// {
	// 	TTimeStamp* tts = new TTimeStamp(v_cd.at(i).year,v_cd.at(i).month, v_cd.at(i).day, v_cd.at(i).hour, v_cd.at(i).minute, v_cd.at(i).second, v_cd.at(i).nanosecond );
	// 	tl1[i] = new TLine(	tts->AsDouble(), 0 ,tts->AsDouble() ,1 )	;
	// 	tl1[i]->Draw("SAME");
	// }

	// c2->cd();
	// tg2->Draw("ap");
	// for(int i = 0; i < v_cd.size() - 1 ; i++ )
	// {
	// 	TTimeStamp* tts = new TTimeStamp(v_cd.at(i).year,v_cd.at(i).month, v_cd.at(i).day, v_cd.at(i).hour, v_cd.at(i).minute, v_cd.at(i).second, v_cd.at(i).nanosecond );
	// 	tl2[i] = new TLine(	tts->AsDouble(), -0.01,tts->AsDouble() , 0.015 )	;
	// 	tl2[i]->Draw("SAME");
	// }



	// TCanvas* 	c3 = new TCanvas("c3","Calibration Parameters in time",900,600);
	// c3->cd();
	// TMultiGraph *mg = new TMultiGraph();
	// mg->Add(tg0);
	// mg->Add(tg1);
	// mg->Add(tg2);
	// // for(int i = 0; i < v_cd.size() - 1 ; i++ )
	// // {
	// // 	TTimeStamp* tts = new TTimeStamp(v_cd.at(i).year,v_cd.at(i).month, v_cd.at(i).day, v_cd.at(i).hour, v_cd.at(i).minute, v_cd.at(i).second, v_cd.at(i).nanosecond );
	// // 	tl2[i] = new TLine(	tts->AsDouble(), -0.01,tts->AsDouble() , 0.015 )	;
	// // 	mg->Add(tl2[i]);
	// // }
	// mg->Draw("ap");

}