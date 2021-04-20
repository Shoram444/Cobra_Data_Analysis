#include <iostream>
#include <dirent.h>
#include <string>
#include <math.h> 
#include <vector>
#include <stdio.h>
// #include "TROOT.h"  
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
const double W [5]	= { 113.903365, 127.904461, 69.92532,				//molar masses of isotopes: 108Cd, 114Cd, 128Te, 
						129.906223, 115.904763} ; // 				70Zn, 64Zn, 120Te, 130Te, 106Cd, 116Cd
																					
const double a [5]  = { 28.73, 31.69, 0.61, 33.8, 7.5};	//natural abundances of isotopes: 108Cd, 114Cd, 128Te, 
																					// 				70Zn, 64Zn, 120Te, 130Te, 106Cd, 116Cd
const double Q [5]  = { 542.5, 866.5, 997.1, 						//Q of isotopes: 108Cd, 114Cd, 128Te, 
						2527.5, 2813.4};							// 				70Zn, 64Zn, 120Te, 130Te, 106Cd, 116Cd

const char *Iso[5]  = { "114Cd", "128Te", "70Zn",							//Names of isotopes: 108Cd, 114Cd, 128Te, 
						"130Te","116Cd"	};								// 				70Zn, 64Zn, 120Te, 130Te, 106Cd, 116Cd

const double l_d[64]  = { 0.00158074, 0.0014663, 0.00152984, 0.00224158, 0.00243734, 0.00214086, 0.00135858,
						  0.00236419, 0.00164415, 0.00134696, 0.0018021, 0.000654638, 0.0018542, 0.00145634,
						  0.00127227, 0.00147545, 0.00257315, 0.0020255, 0.00216535, 0.00154364, 0.00253963, 
						  0.00180284, 0.00111421, 0.00173198, 0.0010829, 0.0022111, 0.00155617, 0.00166385, 
						  0.00251095, 0.00269775, 0.0019247, 0.00159094, 0.00240441, 0.00118512, 0.00125396, 
						  0.00124772, 0.00100205, 0.00125668, 0.000190051, 0.00102673, 0.00123903, 0.00146379, 
						  0.0011961, 0.00102498, 0.000977185, 0.00118687, 0.00120023, 0.000834888, 0.00145179, 
						  0.000992457, 0.000624736, 0.00104622, 0.00233785, 0.00136586, 0.00163178, 0.00197096, 
						  0.00178988, 0.0013888, 0.00143418, 0.00143628, 0.000901274, 0.000624002, 0.000562265, 0.00088992	};	

const double C_d[64] = { 527.325, 197.642, 298.409, 307.837, 321.424, 213.956, 187.341, 282.945, 123.539, 89.9317, 
						 109.042, 82.6086, 127.132, 172.602, 101.422, 160.598, 378.357, 309.254, 733.889, 156.869, 
						 736.482, 133.94, 111.833, 244.343, 159.792, 754.77, 195.365, 82.2015, 484.799, 322.938, 
						 179.4, 178.122, 686.881, 595.642, 230.749, 340.665, 410.682, 677.329, 204.974, 216.652, 
						 189.433, 492.09, 395.591, 279.557, 50.9068, 212.356, 131.343, 101.704, 295.492, 145.803, 61.1561, 
						 121.932, 505.175, 662.733, 343.429, 247.893, 346.578, 958.394, 100.605, 199.643, 148.148, 151.356, 
						 174.474, 215.472 };	



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
	bool 			p_hiz;


	double Get_p_res(double _E)
	{
		if(p_p0 != -1.0 || p_p1 != -1.0 || p_p2 != -1.0 )
		{
			double 		p_res 			=  0; 
		    int         power           =  2;
		    double      p1_times_sqrt_E = p_p1*sqrt(_E);  
		    double      p2_times_E      = p_p2*_E;  

		    p_res      = sqrt( pow(p_p0,power) + pow(p1_times_sqrt_E,power) + pow(p2_times_E, power) );
		    return p_res;			
		}
		else
		{
			double p_res      = 0 ;
		    return p_res;
		}
	}
};

// double 					Get_Resolution(int _d_id , double _E, data_partition _p , vector<CO_detector>  _v_det = vector<CO_detector>()); 
vector<CO_detector>  	Fill_CO_detector(TChain* _T_Chain);
vector<data_partition>  create_partitions(TChain* _T_Chain_det, int _dno);
TH1F* 					hist_style(string _title);


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

				a_det[det_Id -1 ].Add( det_Id, c_mas , c_dur);
			}
		}
		else
		{
			int 		c_dno = _v_det.at(i).get_c_dno();
			double 		c_mas = _v_det.at(i).get_c_mas();
			Double_t	c_dur = _v_det.at(i).get_c_dur();

			a_det[c_dno - 1].Add( c_dno, c_mas, c_dur);
		}
	}

	if(!_all)
	{
  		total_exposure = a_det[det_Id - 1].calc_exposure_kgy();
	}
	else
	{
  		for( int d = 0; d < n_of_det; d++)			///WATCH OUT! Detector numbering starts from 1, so there has to be n+1 n_of_det
		{
			total_exposure += a_det[d - 1].calc_exposure_kgy();
		}
	}

  	return total_exposure;
}
///////FUNCTIONS - from MP.cpp ////////////

double Get_Resolution(int _d_id , double _E, data_partition _p , vector<CO_detector>  _v_det = vector<CO_detector>())
{
	double d_res;
	if(_v_det.size() > 1)
	{
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
	else
	{
		d_res = _p.Get_p_res(_E);
		return d_res;
	}

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


vector<data_partition> create_partitions(TChain* _T_Chain_det, int _dno)
{
	vector<data_partition> v_part;
	int    	  sec_per_year = 86400*365;	//this is used in calculating exposure in units kgy
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

	for(int i = 0; i <2; i++)
	{
		p_part->p_dno = _dno;
		p_part->p_p0  =  -1.0;
		p_part->p_p1  =  -1.0;
		p_part->p_p2  =  -1.0;
		p_part->p_exp =   0  ;

		if(i%2==0)
		{
			p_part->p_hiz =   true  ;
		}
		else
		{
			p_part->p_hiz =   false  ;
		}

	
		v_part.push_back(*p_part);
	}

	


	for(unsigned int i = 0; i < _T_Chain_det->GetEntries(); i++ ) //
	{
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
					v_part.at(k).p_exp 		= (dur * mas->at(j) / sec_per_year) * 0.5 ; //since there are always 2 partitions, one with high z and one with low z. The exposure must be halved.
					v_part.at(k).p_p0 		= p0->at(j);
					v_part.at(k).p_p1 		= p1->at(j);
					v_part.at(k).p_p2 		= p2->at(j);

					v_part.at(k+1).p_dno 		= dno->at(j);
					v_part.at(k+1).p_sta_tim 	= tim;
					v_part.at(k+1).p_exp 		= (dur * mas->at(j) / sec_per_year) * 0.5;
					v_part.at(k+1).p_p0 		= p0->at(j);
					v_part.at(k+1).p_p1 		= p1->at(j);
					v_part.at(k+1).p_p2 		= p2->at(j);

				}

				else if( 	v_part.at(k).p_p0 == p0->at(j) 		&& //if all the parameters are the same in a partition, add exposure and end time as start of the partition + its duration
							v_part.at(k).p_p1 == p1->at(j)     	&&
							v_part.at(k).p_p2 == p2->at(j)
					   )
				{
					v_part.at(k).p_exp		+= (dur * mas->at(j) / sec_per_year) * 0.5;
					v_part.at(k).p_end_tim 	= tim + dur; 		

					v_part.at(k+1).p_exp		+= (dur * mas->at(j) / sec_per_year) * 0.5;
					v_part.at(k+1).p_end_tim 	= tim + dur; 

				}

				else if(	v_part.at(k).p_p0 != p0->at(j)     || // if either of the parameters has changed, create new partition
							v_part.at(k).p_p1 != p1->at(j)	   ||
							v_part.at(k).p_p2 != p2->at(j)
					   )
				{
					k+=2;

					temp_part->p_dno 		= dno->at(j);
					temp_part->p_sta_tim 	= tim;
					temp_part->p_exp 		= (dur * mas->at(j) / sec_per_year) * 0.5;
					temp_part->p_p0 		= p0->at(j);
					temp_part->p_p1 		= p1->at(j);
					temp_part->p_p2 		= p2->at(j);
					temp_part->p_hiz 		= true;

					v_part.push_back(*temp_part);

					temp_part->p_hiz 		= false;

					v_part.push_back(*temp_part);
				}
			}
			delete temp_part;
		}
	}

	return v_part;
}

vector<data_partition> Read_partitions_from_file(string _fileName) //
{
	vector<string> lines =  lines_from_file(_fileName);
	vector<data_partition> v_part;
	data_partition p_part;

	vector<string> temp_data;

	for(int i = 14; i < lines.size()-1; i++)
	{
		size_t found = lines.at(i).find(";");
		if (found != std::string::npos) 
		{
			temp_data.push_back(lines.at(i).substr(0, found));
		}	
		else
		{
			temp_data.push_back(lines.at(i));
		}
	}
	for( int i = 0; i < temp_data.size() ; i++ )
	{
		switch(i%9)
		{
			case 0:
				p_part.p_dno 		= stoi(temp_data.at(i));
				break;
			case 2:
				p_part.p_sta_tim 	= stod(temp_data.at(i));
				break;
			case 3:				
				p_part.p_end_tim 	= stod(temp_data.at(i));
				break;
			case 5:				
				p_part.p_exp 		= stod(temp_data.at(i));
				break;
			case 6:				
				p_part.p_p0 		= stod(temp_data.at(i));
				break;
			case 7:				
				p_part.p_p1 		= stod(temp_data.at(i));
				break;
			case 8:				
				p_part.p_p2 		= stod(temp_data.at(i));
				break;
		}
		if(i%9==0 && i > 0)
		{
			v_part.push_back(p_part);
		}
	}
	return v_part;

}

vector<data_partition> Fill_Partitions(TChain* _T_Chain_eve, vector<data_partition> _v_part)
{

	vector<float>* 		ene = new vector<float>();
	vector<float>* 		ztc = new vector<float>();
	vector<float>* 		aoe = new vector<float>();
	vector<bool>*		fip = new vector<bool>();
	vector<bool>*		fbp = new vector<bool>();
	vector<double>*		eti = new vector<double>(); //time of the event had to be renamed to eti, cause tim is already present in detector tree
	vector<int>*		det = new vector<int>();
	vector<int>*		eid = new vector<int>();

	_T_Chain_eve->SetBranchAddress("cal_edep", &ene);
	_T_Chain_eve->SetBranchAddress("cal_ipos_ztc", &ztc);
	_T_Chain_eve->SetBranchAddress("cal_cpg_diff_AoE", &aoe);
	_T_Chain_eve->SetBranchAddress("flag_injected_pulse", &fip);
	_T_Chain_eve->SetBranchAddress("flag_bad_pulse", &fbp);
	_T_Chain_eve->SetBranchAddress("info_systime", &eti);
	_T_Chain_eve->SetBranchAddress("cal_det", &det);
	_T_Chain_eve->SetBranchAddress("info_idx", &eid);

	sanity_check(_T_Chain_eve);


	for(unsigned int j = 0; j < _T_Chain_eve->GetEntries(); j++)
	{
		if(j%10000==0) cout << j << " of " << _T_Chain_eve->GetEntries() << " Read!" << endl;
		_T_Chain_eve->GetEntry(j);

		for(unsigned int k = 0; k < ene->size(); k++) ///MAROS EVENT
		{
			if( !fip->at(k) && !fbp->at(k) && ztc->at(k) < 1.0 && aoe->at(k) < 1.4)
			{
				CO_event e( ene->at(k), ztc->at(k), aoe->at(k), eti->at(k)  ,
							                det->at(k), eid->at(k), fip->at(k), fbp->at(k)	);
				e.InitCuts(0.2 , 0.95 , 0.872 , 1.3 , false , false ); // double _ztc_min , double _ztc_max , double _aoe_min , double _aoe_max ,
								                                        // bool   _fip , bool   _fbp ;

				if(e.Passed())
				{
					for(int l = 0; l < _v_part.size() ; l++)
					{
						// cout << " z = " << e.Get_c_ztc() << endl;
						if(	e.Get_c_tim() >  _v_part.at(l).p_sta_tim &&
						   	e.Get_c_tim() <  _v_part.at(l).p_end_tim &&
						   	e.Get_c_det() == _v_part.at(l).p_dno 	 &&
						   	e.Get_c_ztc() >= 	0.6					 &&
						   	_v_part.at(l).p_hiz								//z-cut partitioning to high and low partitions. 
	 					   )
						{
							_v_part.at(l).p_ene.push_back(e.Get_c_ene());
							_v_part.at(l).p_hiz = true;
							// cout << e.Get_c_ztc() << " z in partition: " << _v_part.at(l).p_hiz << endl;
						}
						else if( e.Get_c_tim() >  _v_part.at(l).p_sta_tim &&
							   	 e.Get_c_tim() <  _v_part.at(l).p_end_tim &&
							   	 e.Get_c_det() == _v_part.at(l).p_dno 	  &&
							   	 e.Get_c_ztc() < 	0.6					  &&
							   	 !_v_part.at(l).p_hiz
								)
						{
							_v_part.at(l).p_ene.push_back(e.Get_c_ene());
							_v_part.at(l).p_hiz = false;
							// cout << e.Get_c_ztc() << " z in partition: " << _v_part.at(l).p_hiz << endl;
						}
					}
				}
			}
		}

		ene->clear();
		ztc->clear();
		aoe->clear();
		fip->clear();
		fbp->clear();
		eti->clear();
		det->clear();
		eid->clear();
	}

	return _v_part;
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
	int 	detector_ID  = 	   1;
	int 	gaus_cutoff  =   5;
	double 		  k_min  =     0;
	double 	      k_stp  =   0.1;
	double 		  k_max  =    10;
	int 		  k_tot  = int( ( k_max - k_min ) / k_stp ); 
	double 		  k_sig;

	double 		  min_E  = 0;
	double 		  max_E  = 0;

	const char*  txt_outf = "Data_Partitions/Detector_Partitions_info-w_zcut-20210419_1826.txt";
	const char*  h_p_outf = "20210419_1826-Histograms_ksigma.root";
	const char*  M3__outf = "20210419_1826-Method_3.root";
	const char*  M2__outf = "20210419_1826-Graphs_ksigma.root";
	const char*  g_b_outf = "20210419_1826-Bounds-Graphs.root";

	const char*  M3_png   = "20210419_1826-Method_3-dfit.png";

	double 	entries 		[9][k_tot];
	double	sensitivities 	[9][k_tot];
	double  t_half 			[9][k_tot];
	int 	min_bin_ID 		[9];
	int 	max_bin_ID 		[9];
	int 	Q_bin_ID		[9];
	double 	resolution 		[9]; 

	int 	d_num			[64];

	vector<string> 			root_file_path 	= 		ReadFiles(false);
	TChain* 			 	T_Chain_det		=  		Make_TChain(root_file_path, runs 	  );
	TChain* 			 	T_Chain_eve	    =  		Make_TChain(root_file_path, merged_cal);

	// CREATION OF PARTITIONS 					 	/////////////////////

	vector<data_partition> temp_v_part[64];
	vector<data_partition> v_part;

	cout << " CREATING PARTITIONS " << endl;

	for(int i = 0; i < 64; i++)
	{
		temp_v_part[i] = create_partitions(T_Chain_det, i+1);
	}

	for (int i = 0; i < 64; i++)
	{
		for(int j = 0; j < temp_v_part[i].size(); j++)
		{
			v_part.push_back(temp_v_part[i].at(j));
		}
	}

	cout << " PARTITIONS CREATED" << endl;

	cout << " FILLING PARTITIONS " << endl;

	v_part 						  = Fill_Partitions(T_Chain_eve, v_part);


	int total_partitions = 0;
	double total_p_exp = 0;
	// vector<data_partition>       v_part[64];

	cout << " =====================================" << endl;
	cout << " p \t | dno \t | start \t | end \t | exp \t | p0 \t | p1 \t | p2 \t | " << endl;

	for (int j = 0; j < v_part.size(); j++)
	{
		cout 	<< j 						<< " \t|"
				<< v_part.at(j).p_dno 		<< " \t|" 
				<< v_part.at(j).p_sta_tim 	<< " \t|"
				<< v_part.at(j).p_end_tim 	<< " \t|"
				<< v_part.at(j).p_exp	 	<< " \t|"
				<< v_part.at(j).p_p0	 	<< " \t|"
				<< v_part.at(j).p_p1	 	<< " \t|"
				<< v_part.at(j).p_p2	 	<< " \t|"
				<< v_part.at(j).p_ene.size()<< " \t|" << endl;


		total_p_exp += v_part.at(j).p_exp;
		total_partitions = v_part.size();
	}
	cout << " =====================================" << endl;

	cout << " Number of partitions created : " 	<< total_partitions << endl;
	cout << " Total Exposure : " 				<< total_p_exp << endl;

	TH1F* 		h_p[total_partitions];
	// // TCanvas* 	c_p[total_partitions];
	
	ofstream myfile;
	myfile.open (txt_outf);
	myfile << "detector number; partition number; start time; end time; no of events inside; exposure; p0;p1;p2;high z; \n";

	for(int p = 0; p < v_part.size(); p++)
	{
		stringstream 			 	   h_title;			//names for histograms		
		h_title  <<  "h_" <<   v_part.at(p).p_dno << "-" << p;
		string h_tit_str  = 	 h_title.str();

		h_p[p] = hist_style(h_tit_str);

		for(int j = 0; j < v_part.at(p).p_ene.size(); j++)
		{
			h_p[p]->Fill(v_part.at(p).p_ene.at(j));
		}

		myfile  << v_part.at(p).p_dno 			<< "; "
				<< p 							<< "; "
				<< v_part.at(p).p_sta_tim		<< "; "
				<< v_part.at(p).p_end_tim		<< "; "
				<< v_part.at(p).p_ene.size()	<< "; "
				<< v_part.at(p).p_exp			<< "; "
				<< v_part.at(p).p_p0 			<< "; "
				<< v_part.at(p).p_p1 			<< "; "
				<< v_part.at(p).p_p2 			<< "; "
				<< v_part.at(p).p_hiz 			<< "\n";

	}
	// }	
	myfile.close();



	TFile* tf = new TFile(h_p_outf, "RECREATE");
	for (int i = 0; i < total_partitions; i++)
	{
		h_p[i]->Write();
	}
	delete tf;



 	///////////////////////////////////////  T1/2 Optimized Window Counting - Detector Fit ///////// 


	MPFeldman_Cousins* 		mpfc 				= 		new MPFeldman_Cousins(gaus_cutoff, 0.9);


/////////////////////////////////////////////
/////////////PARTITION VARIABLES
/////////////////////////////////////////////

	double C____p[total_partitions];    		///N0 of the partition
	double sigm_p[total_partitions];			///Partition Resolution
	double epsd_p[total_partitions];


cout <<"PARTITION VARIABLES" << endl;


/////////////////////////////////////////////
/////////////Calculate the Quality Factor
/////////////////////////////////////////////
cout << "++++++ Calculating Quality Factor ++++"<< endl;
double q_glob[5];
double q_part 		= 0.0;

double left_m 		= 0.0;
double q__max[5];


double temp_denom 	= 0.0;
double temp_exp 	= 0.0 ; // exponent of beta  --lambda_1 * (Q[2] + lambda_1 * v_part.at(p).Get_p_res( Q[2] ) * 
							//												     v_part.at(p).Get_p_res( Q[2] ) / 2 )
							// 	-> - lambda (Q + lambda*sigma^2 / 2)

for( int q = 0; q < 5; q++)
{
	q_glob[q] = 0.0;
	q__max[q] = 0.0;

	for(unsigned int p = 0; p < total_partitions ; p++)
	{
		C____p[p] = C_d[v_part.at(p).p_dno - 1] * v_part.at(p).p_exp * h_p[p]->GetBinWidth(1);
		sigm_p[p] = v_part.at(p).Get_p_res( Q[q] );
		epsd_p[p] = 1.0;

		temp_exp   	= l_d[v_part.at(p).p_dno - 1] * ( ( l_d[v_part.at(p).p_dno - 1] * sigm_p[p] * sigm_p[p] ) / 2 + Q[q] );

		if(sigm_p[p] != 0)
		{
			temp_denom 	= C____p[p] * sigm_p[p] * sqrt( 2 * TMath::Pi() );
			q_part 		= ( epsd_p[p] / temp_denom ) * TMath::Exp( temp_exp );
	 		q_glob[q]  += q_part / total_partitions;
		}
		if( q_part > q__max[q])
		{
			q__max[q] = q_part;
		}

	}
	cout << "q_glob = " << q_glob[q] << endl;
	cout << "+++++++++++++++++++++++++++++++++++++ "<< endl;
}




/////////////////////////////////////////////
/////////////RHO OVER BETA
/////////////////////////////////////////////

	int    rob_num; 	
	double rob_max[5];									/// rho/beta maximum value
	double rob_stp[5]; 									/// rho/beta step size
	double rob_min[5]; 									/// rho/beta factor

	for(int q = 0; q < 5; q++ )
	{
		rob_num    = 1e3;
		rob_max[q] = 1*q__max[q];
		rob_stp[q] = rob_max[q] / double (rob_num);
		rob_min[q] = rob_stp[q];
	}

	double rob     = 0.00;
cout <<" RHO OVER BETA" << endl;

/////////////////////////////////////////////
/////////////Temporary VARIABLES FOR bounds
/////////////////////////////////////////////

	double gamma 		= 0.0;   // gamma = N0_p*sigm_p*sqrt(2*Pi)*rb/epsd_p 		(inside of ln)
	double sl 	 		= 0.0;	 // sl   = sigm_p*l_d[v_part.at(p).p_dno - 1]                    		(sl - sigma * lambda)
	double lq   		= 0.0;   // lq   = lambda*Q 		   						(sl - lambda * Q )
	double in_sqrt   	= 0.0;   // in_sqrt = (sl)^2 + 2*sl - 2*ln(gamma)/sigma^2  	(inside of the square root term)

cout <<"Temporary VARIABLES FOR bounds" << endl;

/////////////////////////////////////////////
/////////////VARIABLES FOR FINAL T1/2 CALCULATION
/////////////////////////////////////////////


	double*** righ_p = new double**[5]; 				///Right bound Dynamically allocated because memory fails
	double*** left_p = new double**[5]; //[total_partitions][rob_num];		///Left bound

	for (int q = 0; q < 5; q++)
	{
		left_p[q] = new double*[total_partitions];
		righ_p[q] = new double*[total_partitions];
		for(int p = 0 ; p < total_partitions; p++)
		{
			left_p[q][p] = new double[rob_num];
			righ_p[q][p] = new double[rob_num];
		}
	}

	double erf_r		;   // The inside of the error function with r_p
	double erf_l 		;	// The inside of the error function with l_p
	double FC_b 		;   // The inside of the FC sensitivity , number of backgrounds
	
	double numerator 	;	// Numerator of the T1/2 eq				
	double denumerator 	;   // Denumberator of the T1/2 eq
	double constants    ;   // The rest of the T1/2 eq

	double T_half[5]    ;

	double T_half_max = 0.0;
	double bkg_counts = 0.0;

	double mid_val = 0.0;
	double sig_max = 0.0;

	TGraph* Opt_win_Gr[5] ;
	for( int q = 0; q < 5; q++)
	{
		Opt_win_Gr[q] = new TGraph();
	}

cout <<"VARIABLES FOR FINAL T1/2 CALCULATION" << endl;


/////////////////////////////////////////////
/////////////THALF CALCULATING
/////////////////////////////////////////////
	cout <<"++++++++++ THALF START +++++++++" << endl;
	for( int q = 0; q < 5; q++)
	{
		for(int rb = 0; rb < rob_num; rb++)
		{
			rob 			= (rb * rob_stp[q]) + rob_min[q];

			numerator 		= 0.0;
			denumerator 	= 0.0;
			constants    	= 0.0;
			erf_l 			= 0.0;
			erf_r			= 0.0;
			FC_b			= 0.0;


				for(int p = 0; p < total_partitions; p++)
				{
					C____p[p] = C_d[v_part.at(p).p_dno - 1] * v_part.at(p).p_exp * h_p[p]->GetBinWidth(1) ;
					sigm_p[p] = v_part.at(p).Get_p_res( Q[q] );
					epsd_p[p] = 1.0;

					if(sigm_p[p] != 0) //some partitions do not have resolution parameters defined, so they will not be accounted for
					{
						temp_exp   = l_d[v_part.at(p).p_dno - 1] * ( ( l_d[v_part.at(p).p_dno - 1] * sigm_p[p] * sigm_p[p] ) / 2 + Q[q] );
						temp_denom = C____p[p] * sigm_p[p] * sqrt( 2 * TMath::Pi() );

				 		q_part = ( epsd_p[p] / temp_denom ) * TMath::Exp( temp_exp );
					}
					else
					{
				 		q_part = 0.0;
					}

					if( q_part > ( 0.2 * q_glob[q] ) )
					{
						gamma 		= ( C____p[p] * sigm_p[p] * sqrt( 2 * TMath::Pi() ) * rob ) / epsd_p[p] ;
						sl 			= sigm_p[p] * l_d[v_part.at(p).p_dno - 1] ; 
						lq 			= ( l_d[v_part.at(p).p_dno - 1]  * Q[q] ) ;
						in_sqrt 	= ( sl * sl + 2*lq - 2*log(gamma) ) ;

						left_p[q][p][rb] 	= Q[q] + ( sl * sigm_p[p] ) - ( sigm_p[p] * sqrt(in_sqrt) );
						righ_p[q][p][rb] 	= Q[q] + ( sl * sigm_p[p] ) + ( sigm_p[p] * sqrt(in_sqrt) );

						erf_r 		= (righ_p[q][p][rb] - Q[q]) / ( sqrt(2) * sigm_p[p] );
						erf_l 		= (left_p[q][p][rb] - Q[q]) / ( sqrt(2) * sigm_p[p] );

						mid_val 	= Q[q] + ( sl * sigm_p[p] );

						if(in_sqrt > 0)   ///if the inside of the sqrt were negative, the partition is to be ignored
						{
							if(left_p[q][p][rb] > Q[q])
							{
								numerator  += v_part.at(p).p_exp * (epsd_p[p] / 2) * ( erf(erf_r) - erf(erf_l) );
							}
							else
							{
								numerator  += v_part.at(p).p_exp * (epsd_p[p] / 2) * ( erf(erf_r) + erf(erf_l) );
							}

							if( left_p[q][p][rb] < righ_p[q][p][rb] )
							{
								FC_b += ( C____p[p]  /  l_d[v_part.at(p).p_dno - 1]  ) * 
										( v_part.at(p).p_exp      ) *
										( TMath::Exp(-l_d[v_part.at(p).p_dno - 1] * left_p[q][p][rb]) - TMath::Exp(-l_d[v_part.at(p).p_dno - 1] * righ_p[q][p][rb]) );
							}
							else
							{
								FC_b += -1*( C____p[p]  /  l_d[v_part.at(p).p_dno - 1]  ) * 
										( v_part.at(p).p_exp      ) *
										( TMath::Exp(-l_d[v_part.at(p).p_dno - 1] * left_p[q][p][rb]) - TMath::Exp(-l_d[v_part.at(p).p_dno - 1] * righ_p[q][p][rb]) );
							}
							
							if( left_p[q][p][rb] > left_m)
							{
								left_m = left_p[q][p][rb];
								// rb___m = rb;
							}
							if( sigm_p[p] > sig_max )
							{
								sig_max = sigm_p[p];
							}
						}
						else
						{
							left_p[q][p][rb]    = Q[q] + ( sl * sigm_p[p] );
							righ_p[q][p][rb]    = Q[q] + ( sl * sigm_p[p] );
							FC_b       		+= 0;
							numerator  		+= 0;
						}

						gamma 			= 0.0 ;  //reset values
						sl 				= 0.0 ;
						sl 				= 0.0 ;
						in_sqrt 		= 0.0 ;
						
					}
					else
					{

						left_p[q][p][rb]   	= 0.0 ;
						righ_p[q][p][rb]   	= 0.0 ;
						sigm_p[p]   		= 0.0 ;
						C____p[p]   		= 0.0 ;
						epsd_p[p]   		= 0.0 ;
						FC_b       			+= 0;
						numerator  			+= 0;
					}
				}

			constants   = ( N_A * ln2 ) / ( W[q] * a[q] );

			if( FC_b < gaus_cutoff )
			{
				denumerator = mpfc->get_sensitivity( FC_b );
			}
			else
			{
				denumerator = sqrt(FC_b);
			}

			T_half[q] 		= constants * ( numerator / denumerator );

			if( T_half_max < T_half[q])
			{
				T_half_max = T_half[q];
				bkg_counts = FC_b;
			}

			if(rb%100 == 0)
			{
				// cout << "numerator "   << numerator << endl;
				// cout << "denumerator " << denumerator << endl;
				cout << "+++++++++++++++++++" << endl;
				cout << "T_half = " << T_half[q] << endl;
				cout << "+++++++++++++++++++" << endl;
			}

			Opt_win_Gr[q]->SetPoint(rb, rob, T_half[q]);

		}
	}
	cout <<"++++++++++ THALF END +++++++++" << endl;

	TCanvas* Opt_win_Canv[5];
	for(int q = 0; q < 5 ; q++)
	{
		stringstream 			q_name;
		q_name << "Calculated Half-Life for " << Iso[q] << " ;rho/beta [kgy]; half-life [y]" ;
		string s_Q_name = q_name.str();


		Opt_win_Canv[q]= new TCanvas(s_Q_name.c_str(),s_Q_name.c_str() ,1000 ,600 );

		stringstream 			M3_fname;
		M3_fname << "./Figures/Method3/" << Iso[q] << M3_png ;
		string s_M3_fname = M3_fname.str();

		Opt_win_Gr[q]->SetMarkerStyle(21 + q);
		Opt_win_Gr[q]->SetMarkerColor(2 + q);
		Opt_win_Gr[q]->SetTitle(s_Q_name.c_str());

		gROOT->SetBatch(kTRUE);

		Opt_win_Canv[q]->cd();
		Opt_win_Gr[q]->Draw("apl");
		Opt_win_Canv[q]->SaveAs(s_M3_fname.c_str());

		cout << " Best T_half = " 			<< T_half_max << endl
			 << " Expected background = " 	<< bkg_counts << endl;
	}

	

///////////////////////////////////////////////
/////////// GRAPHS OF R, L vs RB
///////////////////////////////////////////////
	TGraph*  		gr_l[5][total_partitions];
	TGraph*  		gr_r[5][total_partitions];
	TGraph*  		gr_q[5][total_partitions];
	TMultiGraph*    mg_b[5][total_partitions];
	TCanvas* 		c_l [5][total_partitions];

	TFile*   tf_bounds = new TFile(g_b_outf, "RECREATE");

	for( int q = 0; q < 5; q++)
	{
		for (int  p = 0; p  < total_partitions; p++ )
		{	
			stringstream  c_l_title;
			c_l_title  << "./Figures/Bounds/" << Iso[q] << p << ".png" ;
			string c_l_tit_str  = 	 c_l_title.str();

			c_l[q][p]  = new TCanvas(c_l_tit_str.c_str(), c_l_tit_str.c_str(), 1000, 600 );
			gr_l[q][p] = new TGraph();
			gr_r[q][p] = new TGraph();
			gr_q[q][p] = new TGraph();
			mg_b[q][p] = new TMultiGraph();

			for (int rb = 0; rb < rob_num; rb++)
			{
				rob  = (rb * rob_stp[q]) + rob_min[q];

				gr_l[q][p]->SetPoint(rb, rob, left_p[q][p][rb]);
				gr_r[q][p]->SetPoint(rb, rob, righ_p[q][p][rb]);
				gr_q[q][p]->SetPoint(rb, rob, Q[q]);
			}

			// gr_l[p]->SetMarkerStyle(22);
			gr_l[q][p]->SetLineColor(3);
			gr_l[q][p]->SetLineWidth(2);
			gr_l[q][p]->SetTitle("Bounds vs rho/beta; rho/beta [kgy]; energy [keV]");
			gr_l[q][p]->GetYaxis()->SetRangeUser(int(Q[q] - 0.5*Q[q]), int(Q[q] + 0.5*Q[q]));
			gr_l[q][p]->GetXaxis()->SetRangeUser(0, 1.1*q__max[q]);
			mg_b[q][p]->Add(gr_l[q][p]);



			// gr_r[p]->SetMarkerStyle(23);
			gr_r[q][p]->SetLineColor(4);
			gr_r[q][p]->SetLineWidth(2);
			gr_r[q][p]->GetYaxis()->SetRangeUser(int(Q[q] - 0.5*Q[q]), int(Q[q] + 0.5*Q[q]));
			gr_r[q][p]->GetXaxis()->SetRangeUser(0, 1.1*q__max[q]);
			mg_b[q][p]->Add(gr_r[q][p]);



			// gr_q[p]->SetMarkerStyle(q);
			gr_q[q][p]->SetLineColor(5);
			gr_q[q][p]->SetLineWidth(2);
			mg_b[q][p]->Add(gr_q[q][p]);
			mg_b[q][p]->SetTitle("Bounds vs rho/beta; rho/beta [kgy]; energy [keV]");


			c_l[q][p]->cd();
			mg_b[q][p]->Draw("al");
			// gr_r[p]->Draw("SAME");
			// gr_q[p]->Draw("SAME");

			if( q == 4  )
			{
				c_l[q][p]->SaveAs(c_l_tit_str.c_str());
			}

			// gr_l[p]->Draw("apl");
			// gr_r[p]->Write();

		}
	}

	delete tf_bounds;

	// gROOT->SetBatch(kTRUE);


///////////////////////////////////////////////
/////////// GRAPHS OF R, L vs RB
///////////////////////////////////////////////

	TH2D* h2d[5]; 
	TCanvas* c2d[5];

	for(int q = 0; q < 5; q++)
	{
		stringstream  c_b_title;
		c_b_title  << "c2d_dfit" << Iso[q] <<".png" ;
		string s_c_b_title  = 	 c_b_title.str();

		stringstream  h_b_title;
		h_b_title  << "c2d_dfit" << Iso[q]  ;
		string s_h_b_title  = 	 h_b_title.str();

		h2d[q] = new TH2D(s_h_b_title.c_str(),s_h_b_title.c_str(), rob_num, 0, 1.1 * q__max[q],
										 int(4 * sig_max), Q[q] - 2 * sig_max, Q[q] + 2 * sig_max);

		c2d[q] = new TCanvas(s_c_b_title.c_str(), s_c_b_title.c_str(), 1000, 600 );

		for (int  p = 0; p  < total_partitions; p++ )
		{
			for (int rb = 0; rb < rob_num; rb++)
			{
				rob  = (rb * rob_stp[q]) + rob_min[q];
				if( (righ_p[q][p][rb] - left_p[q][p][rb] ) > 1e-16)
				{
					h2d[q]->Fill(rob, left_p[q][p][rb]);
					h2d[q]->Fill(rob, righ_p[q][p][rb]);
				}
				
			}
		}

		h2d[q]->SetTitle("Bounds vs rho/beta; rho/beta [kgy]; energy [keV]");
		c2d[q]->cd();
		h2d[q]->Draw("COLZ");
		c2d[q]->SaveAs(s_c_b_title.c_str());
	}


 	///////////////////////////////////////  T1/2 Optimized Window Counting - Detector Fit ///////// 

}
///////////////////////////////////////  T1/2 Optimized Window Counting - Global Fit ///////// 


// 	MPFeldman_Cousins* 		mpfc 				= 		new MPFeldman_Cousins(gaus_cutoff, 0.9);

// /////////////////////////////////////////////
// /////////////PARAMETERS OF FIT Region 1
// /////////////////////////////////////////////

// 	const double 	E_Fit_min_1 =  		400.0;  					///Range of the fit
// 	const double 	E_Fit_max_1 = 	   1100.0; 
// 	const int 		Bin_Width_1 = 		   20;						///Bin width used for the fit
// 	const double 	   C_fit    = 2.35078e+04;   					///norm of the exponential bkg fit 		 [C-fit] = 1
// 	const double 	   lambda_1 = 2.61635e-03;						///lambda of the exponential bkg fit
// 	double C = C_fit / (total_p_exp * Bin_Width_1);    				///Norm of the background shape function [C] = 1/(keV*yr*kg)

// cout <<"PARAMETERS OF FIT Region 1" << endl;

// /////////////////////////////////////////////
// /////////////PARTITION VARIABLES
// /////////////////////////////////////////////

// 	double C____p[total_partitions];    		///N0 of the partition
// 	double sigm_p[total_partitions];			///Partition Resolution
// 	double epsd_p[total_partitions];


// cout <<"PARTITION VARIABLES" << endl;


// /////////////////////////////////////////////
// /////////////Calculate the Quality Factor
// /////////////////////////////////////////////
// cout << "++++++ Calculating Quality Factor ++++"<< endl;
// double q_glob[5];
// double q_part 		= 0.0;

// double left_m 		= 0.0;
// double q__max[5];


// double temp_denom 	= 0.0;
// double temp_exp 	= 0.0 ; // exponent of beta  --lambda_1 * (Q[2] + lambda_1 * v_part.at(p).Get_p_res( Q[2] ) * 
// 							//												     v_part.at(p).Get_p_res( Q[2] ) / 2 )
// 							// 	-> - lambda (Q + lambda*sigma^2 / 2)
// for( int q = 0; q < 5; q++)
// {
// 	q_glob[q] = 0.0;
// 	q__max[q] = 0.0;

// 	for(unsigned int p = 0; p < total_partitions ; p++)
// 	{
// 		C____p[p] = C * v_part.at(p).p_exp * h_p[p]->GetBinWidth(1);
// 		sigm_p[p] = v_part.at(p).Get_p_res( Q[q] );
// 		epsd_p[p] = 1.0;

// 		temp_exp   	= lambda_1 * ( ( lambda_1 * sigm_p[p] * sigm_p[p] ) / 2 + Q[q] );

// 		if(sigm_p[p] != 0)
// 		{
// 			temp_denom 	= C____p[p] * sigm_p[p] * sqrt( 2 * TMath::Pi() );
// 			q_part 		= ( epsd_p[p] / temp_denom ) * TMath::Exp( temp_exp );
// 	 		q_glob[q]  += q_part / total_partitions;
// 		}
// 		if( q_part > q__max[q])
// 		{
// 			q__max[q] = q_part;
// 		}
// 	}
// 	cout << "q_glob = " << q_glob[q] << endl;
// 	cout << "+++++++++++++++++++++++++++++++++++++ "<< endl;
// }




// /////////////////////////////////////////////
// /////////////RHO OVER BETA
// /////////////////////////////////////////////

// 	int    rob_num; 	
// 	double rob_max[5];									/// rho/beta maximum value
// 	double rob_stp[5]; 									/// rho/beta step size
// 	double rob_min[5]; 									/// rho/beta factor

// 	for(int q = 0; q < 5; q++ )
// 	{
// 		rob_num    = 1e3;
// 		rob_max[q] = 0.2*q__max[q];
// 		rob_stp[q] = rob_max[q] / double (rob_num);
// 		rob_min[q] = rob_stp[q];
// 	}

// 	double rob     = 0.00;
// cout <<" RHO OVER BETA" << endl;

// /////////////////////////////////////////////
// /////////////Temporary VARIABLES FOR bounds
// /////////////////////////////////////////////

// 	double gamma 		= 0.0;   // gamma = N0_p*sigm_p*sqrt(2*Pi)*rb/epsd_p 		(inside of ln)
// 	double sl 	 		= 0.0;	 // sl   = sigm_p*lambda_1                    		(sl - sigma * lambda)
// 	double lq   		= 0.0;   // lq   = lambda*Q 		   						(sl - lambda * Q )
// 	double in_sqrt   	= 0.0;   // in_sqrt = (sl)^2 + 2*sl - 2*ln(gamma)/sigma^2  	(inside of the square root term)

// cout <<"Temporary VARIABLES FOR bounds" << endl;

// /////////////////////////////////////////////
// /////////////VARIABLES FOR FINAL T1/2 CALCULATION
// /////////////////////////////////////////////


// 	double*** righ_p = new double**[5]; 				///Right bound Dynamically allocated because memory fails
// 	double*** left_p = new double**[5]; //[total_partitions][rob_num];		///Left bound

// 	for (int q = 0; q < 5; q++)
// 	{
// 		left_p[q] = new double*[total_partitions];
// 		righ_p[q] = new double*[total_partitions];
// 		for(int p = 0 ; p < total_partitions; p++)
// 		{
// 			left_p[q][p] = new double[rob_num];
// 			righ_p[q][p] = new double[rob_num];
// 		}
// 	}

// 	double erf_r		;   // The inside of the error function with r_p
// 	double erf_l 		;	// The inside of the error function with l_p
// 	double FC_b 		;   // The inside of the FC sensitivity , number of backgrounds
	
// 	double numerator 	;	// Numerator of the T1/2 eq				
// 	double denumerator 	;   // Denumberator of the T1/2 eq
// 	double constants    ;   // The rest of the T1/2 eq

// 	double T_half[5]    ;

// 	double T_half_max = 0.0;
// 	double bkg_counts = 0.0;

// 	double mid_val = 0.0;
// 	double sig_max = 0.0;

// 	TGraph* Opt_win_Gr[5] ;
// 	for( int q = 0; q < 5; q++)
// 	{
// 		Opt_win_Gr[q] = new TGraph();
// 	}

// cout <<"VARIABLES FOR FINAL T1/2 CALCULATION" << endl;


// /////////////////////////////////////////////
// /////////////THALF CALCULATING
// /////////////////////////////////////////////
// 	cout <<"++++++++++ THALF START +++++++++" << endl;
// 	for( int q = 0; q < 5; q++)
// 	{
// 		for(int rb = 0; rb < rob_num; rb++)
// 		{
// 			rob 			= (rb * rob_stp[q]) + rob_min[q];

// 			numerator 		= 0.0;
// 			denumerator 	= 0.0;
// 			constants    	= 0.0;
// 			erf_l 			= 0.0;
// 			erf_r			= 0.0;
// 			FC_b			= 0.0;


// 				for(int p = 0; p < total_partitions; p++)
// 				{
// 					C____p[p] = C * v_part.at(p).p_exp * h_p[p]->GetBinWidth(1) ;
// 					sigm_p[p] = v_part.at(p).Get_p_res( Q[q] );
// 					epsd_p[p] = 1.0;

// 					if(sigm_p[p] != 0) //some partitions do not have resolution parameters defined, so they will not be accounted for
// 					{
// 						temp_exp   = lambda_1 * ( ( lambda_1 * sigm_p[p] * sigm_p[p] ) / 2 + Q[q] );
// 						temp_denom = C____p[p] * sigm_p[p] * sqrt( 2 * TMath::Pi() );

// 				 		q_part = ( epsd_p[p] / temp_denom ) * TMath::Exp( temp_exp );
// 					}
// 					else
// 					{
// 				 		q_part = 0.0;
// 					}

// 					if( q_part > ( 0.05 * q_glob[q] ) )
// 					{
// 						gamma 		= ( C____p[p] * sigm_p[p] * sqrt( 2 * TMath::Pi() ) * rob ) / epsd_p[p] ;
// 						sl 			= sigm_p[p] * lambda_1 ; 
// 						lq 			= ( lambda_1  * Q[q] ) ;
// 						in_sqrt 	= ( sl * sl + 2*lq - 2*log(gamma) ) ;

// 						left_p[q][p][rb] 	= Q[q] + ( sl * sigm_p[p] ) - ( sigm_p[p] * sqrt(in_sqrt) );
// 						righ_p[q][p][rb] 	= Q[q] + ( sl * sigm_p[p] ) + ( sigm_p[p] * sqrt(in_sqrt) );

// 						erf_r 		= (righ_p[q][p][rb] - Q[q]) / ( sqrt(2) * sigm_p[p] );
// 						erf_l 		= (left_p[q][p][rb] - Q[q]) / ( sqrt(2) * sigm_p[p] );

// 						mid_val 	= Q[q] + ( sl * sigm_p[p] );

// 						if(in_sqrt > 0)   ///if the inside of the sqrt were negative, the partition is to be ignored
// 						{
// 							if(left_p[q][p][rb] > Q[q])
// 							{
// 								numerator  += v_part.at(p).p_exp * (epsd_p[p] / 2) * ( erf(erf_r) - erf(erf_l) );
// 							}
// 							else
// 							{
// 								numerator  += v_part.at(p).p_exp * (epsd_p[p] / 2) * ( erf(erf_r) + erf(erf_l) );
// 							}

// 							if( left_p[q][p][rb] < righ_p[q][p][rb] )
// 							{
// 								FC_b += ( C____p[p]  /  lambda_1  ) * 
// 										( v_part.at(p).p_exp      ) *
// 										( TMath::Exp(-lambda_1 * left_p[q][p][rb]) - TMath::Exp(-lambda_1 * righ_p[q][p][rb]) );
// 							}
// 							else
// 							{
// 								FC_b += -1*( C____p[p]  /  lambda_1  ) * 
// 										( v_part.at(p).p_exp      ) *
// 										( TMath::Exp(-lambda_1 * left_p[q][p][rb]) - TMath::Exp(-lambda_1 * righ_p[q][p][rb]) );
// 							}
// 							// if(p == 1)
// 							// {
// 							// cout << "Q[q] " 	 << Q[q] 	  << endl
// 							// 	 << "in_sqrt" << in_sqrt << endl
// 							// 	 << "( sqrt(in_sqrt) )" << (  sqrt(in_sqrt) ) << endl
// 							// 	 << "sigm_p[p] " << sigm_p[p] << endl
// 							// 	 // << "C____p[p] " << C____p[p] << endl
// 							// 	 // << "gamma " 	 << C____p[p] << endl
// 							// 	 << "sl * sl " 		 << sl * sl << endl
// 							// 	 << "2*sl " 		 << 2*sl << endl
// 							// 	 << "2*log(gamma) " 	 << 2*log(gamma) << endl
// 							// 	 << "left_p[p][rb] " << left_p[p][rb] << endl
// 							// 	 // << "righ_p[p][rb] " << righ_p[p][rb] << endl
// 							// 	 // << "numerator " << numerator << endl
// 							// 	 << "FC_b      " << FC_b      << endl;
// 							// }

// 							if( left_p[q][p][rb] > left_m)
// 							{
// 								left_m = left_p[q][p][rb];
// 								// rb___m = rb;
// 							}
// 							if( sigm_p[p] > sig_max )
// 							{
// 								sig_max = sigm_p[p];
// 							}
// 						}
// 						else
// 						{
// 							left_p[q][p][rb]    = Q[q] + ( sl * sigm_p[p] );
// 							righ_p[q][p][rb]    = Q[q] + ( sl * sigm_p[p] );
// 							FC_b       		+= 0;
// 							numerator  		+= 0;
// 						}

// 						gamma 			= 0.0 ;  //reset values
// 						sl 				= 0.0 ;
// 						sl 				= 0.0 ;
// 						in_sqrt 		= 0.0 ;
// 						// left_p[p][rb]   = 0.0 ;
// 						// righ_p[p][rb]   = 0.0 ;
// 						// sigm_p[p]   	= 0.0 ;
// 						// C____p[p]   	= 0.0 ;
// 						// epsd_p[p]   	= 0.0 ;
// 					}
// 					else
// 					{

// 						left_p[q][p][rb]   	= 0.0 ;
// 						righ_p[q][p][rb]   	= 0.0 ;
// 						sigm_p[p]   		= 0.0 ;
// 						C____p[p]   		= 0.0 ;
// 						epsd_p[p]   		= 0.0 ;
// 						FC_b       			+= 0;
// 						numerator  			+= 0;
// 					}
// 				}

// 			constants   = ( N_A * ln2 ) / ( W[q] * a[q] );

// 			if( FC_b < gaus_cutoff )
// 			{
// 				denumerator = mpfc->get_sensitivity( FC_b );
// 			}
// 			else
// 			{
// 				denumerator = sqrt(FC_b);
// 				// cout << " USED GAUSSIAN - " << endl;
// 				// cout << " FC_b =  "         << FC_b << endl;
// 			}

// 			T_half[q] 		= constants * ( numerator / denumerator );

// 			if( T_half_max < T_half[q])
// 			{
// 				T_half_max = T_half[q];
// 				bkg_counts = FC_b;
// 			}

// 			if(rb%100 == 0)
// 			{
// 				cout << "numerator "   << numerator << endl;
// 				cout << "denumerator " << denumerator << endl;
// 				cout << "+++++++++++++++++++" << endl;
// 				cout << "T_half = " << T_half[q] << endl;
// 				cout << "+++++++++++++++++++" << endl;
// 			}

// 			Opt_win_Gr[q]->SetPoint(rb, rob, T_half[q]);

// 		}
// 	}
// 	cout <<"++++++++++ THALF END +++++++++" << endl;

// 	TCanvas* Opt_win_Canv[5];
// 	for(int q = 0; q < 5 ; q++)
// 	{
// 		stringstream 			q_name;
// 		q_name << "Calculated Half-Life for " << Iso[q] << " ;rho/beta [kgy]; half-life [y]" ;
// 		string s_Q_name = q_name.str();


// 		Opt_win_Canv[q]= new TCanvas(s_Q_name.c_str(),s_Q_name.c_str() ,1000 ,600 );

// 		stringstream 			M3_fname;
// 		M3_fname << "./Figures/Method3/" << Iso[q] << M3_png ;
// 		string s_M3_fname = M3_fname.str();

// 		Opt_win_Gr[q]->SetMarkerStyle(21 + q);
// 		Opt_win_Gr[q]->SetMarkerColor(2 + q);
// 		Opt_win_Gr[q]->SetTitle(s_Q_name.c_str());

// 		gROOT->SetBatch(kTRUE);

// 		Opt_win_Canv[q]->cd();
// 		Opt_win_Gr[q]->Draw("apl");
// 		Opt_win_Canv[q]->SaveAs(s_M3_fname.c_str());

// 		cout << " Best T_half = " 			<< T_half_max << endl
// 			 << " Expected background = " 	<< bkg_counts << endl;
// 	}

	

// ///////////////////////////////////////////////
// /////////// GRAPHS OF R, L vs RB
// ///////////////////////////////////////////////
// 	TGraph*  		gr_l[5][total_partitions];
// 	TGraph*  		gr_r[5][total_partitions];
// 	TGraph*  		gr_q[5][total_partitions];
// 	TMultiGraph*    mg_b[5][total_partitions];
// 	TCanvas* 		c_l [5][total_partitions];

// 	TFile*   tf_bounds = new TFile(g_b_outf, "RECREATE");

// 	for( int q = 0; q < 5; q++)
// 	{
// 		for (int  p = 0; p  < total_partitions; p++ )
// 		{	
// 			stringstream  c_l_title;
// 			c_l_title  << "./Figures/Bounds/" << Iso[q] << p << ".png" ;
// 			string c_l_tit_str  = 	 c_l_title.str();

// 			c_l[q][p]  = new TCanvas(c_l_tit_str.c_str(), c_l_tit_str.c_str(), 1000, 600 );
// 			gr_l[q][p] = new TGraph();
// 			gr_r[q][p] = new TGraph();
// 			gr_q[q][p] = new TGraph();
// 			mg_b[q][p] = new TMultiGraph();

// 			for (int rb = 0; rb < rob_num; rb++)
// 			{
// 				rob  = (rb * rob_stp[q]) + rob_min[q];

// 				gr_l[q][p]->SetPoint(rb, rob, left_p[q][p][rb]);
// 				gr_r[q][p]->SetPoint(rb, rob, righ_p[q][p][rb]);
// 				gr_q[q][p]->SetPoint(rb, rob, Q[q]);
// 			}

// 			// gr_l[p]->SetMarkerStyle(22);
// 			gr_l[q][p]->SetLineColor(3);
// 			gr_l[q][p]->SetLineWidth(2);
// 			gr_l[q][p]->SetTitle("Bounds vs rho/beta; rho/beta [kgy]; energy [keV]");
// 			gr_l[q][p]->GetYaxis()->SetRangeUser(int(Q[q] - 0.5*Q[q]), int(Q[q] + 0.5*Q[q]));
// 			gr_l[q][p]->GetXaxis()->SetRangeUser(0, 1.1*q__max[q]);
// 			mg_b[q][p]->Add(gr_l[q][p]);



// 			// gr_r[p]->SetMarkerStyle(23);
// 			gr_r[q][p]->SetLineColor(4);
// 			gr_r[q][p]->SetLineWidth(2);
// 			gr_r[q][p]->GetYaxis()->SetRangeUser(int(Q[q] - 0.5*Q[q]), int(Q[q] + 0.5*Q[q]));
// 			gr_r[q][p]->GetXaxis()->SetRangeUser(0, 1.1*q__max[q]);
// 			mg_b[q][p]->Add(gr_r[q][p]);



// 			// gr_q[p]->SetMarkerStyle(q);
// 			gr_q[q][p]->SetLineColor(5);
// 			gr_q[q][p]->SetLineWidth(2);
// 			mg_b[q][p]->Add(gr_q[q][p]);
// 			mg_b[q][p]->SetTitle("Bounds vs rho/beta; rho/beta [kgy]; energy [keV]");


// 			c_l[q][p]->cd();
// 			mg_b[q][p]->Draw("al");
// 			// gr_r[p]->Draw("SAME");
// 			// gr_q[p]->Draw("SAME");

// 			if( q == 4  )
// 			{
// 				c_l[q][p]->SaveAs(c_l_tit_str.c_str());
// 			}

// 			// gr_l[p]->Draw("apl");
// 			// gr_r[p]->Write();

// 		}
// 	}

// 	delete tf_bounds;

// 	// gROOT->SetBatch(kTRUE);


// ///////////////////////////////////////////////
// /////////// GRAPHS OF R, L vs RB
// ///////////////////////////////////////////////

// 	TH2D* h2d[5]; 
// 	TCanvas* c2d[5];

// 	for(int q = 0; q < 5; q++)
// 	{
// 		stringstream  c_b_title;
// 		c_b_title  << "c2d" << Iso[q] <<".png" ;
// 		string s_c_b_title  = 	 c_b_title.str();

// 		stringstream  h_b_title;
// 		h_b_title  << "h2d" << Iso[q]  ;
// 		string s_h_b_title  = 	 h_b_title.str();

// 		h2d[q] = new TH2D(s_h_b_title.c_str(),s_h_b_title.c_str(), rob_num, 0, 1.1 * q__max[q],
// 										 int(4 * sig_max), Q[q] - 2 * sig_max, Q[q] + 2 * sig_max);

// 		c2d[q] = new TCanvas(s_c_b_title.c_str(), s_c_b_title.c_str(), 1000, 600 );

// 		for (int  p = 0; p  < total_partitions; p++ )
// 		{
// 			for (int rb = 0; rb < rob_num; rb++)
// 			{
// 				rob  = (rb * rob_stp[q]) + rob_min[q];
// 				if( (righ_p[q][p][rb] - left_p[q][p][rb] ) > 1e-16)
// 				{
// 					h2d[q]->Fill(rob, left_p[q][p][rb]);
// 					h2d[q]->Fill(rob, righ_p[q][p][rb]);
// 				}

// 				// if(p%100 == 0)
// 				// {
// 				// 	cout << "left_p[4][p][rb] " << left_p[4][p][rb] << endl
// 				// 		 << "righ_p[4][p][rb] " << righ_p[4][p][rb] << endl
// 				// 		 << "dif 			  " << left_p[4][p][rb] - righ_p[4][p][rb] << endl;
// 				// }
// 			}
// 		}

// 		h2d[q]->SetTitle("Bounds vs rho/beta; rho/beta [kgy]; energy [keV]");
// 		c2d[q]->cd();
// 		h2d[q]->Draw("COLZ");
// 		c2d[q]->SaveAs(s_c_b_title.c_str());
// 	}


 	///////////////////////////////////////  T1/2 Optimized Window Counting - Global Fit ///////// 

 		// // CREATION OF PARTITIONS 					 	/////////////////////

	// // ----------------------------------------- 		/////////////////////

	// // CALCULATING T1/2 FROM PARTITION HISTOGRAMS	 		 	/////////////////////

	// TH1F* 		h[v_part.size()];

	// MPFeldman_Cousins* 		obj 			= 		new MPFeldman_Cousins(gaus_cutoff, 0.9);
	// TFile* 					tf_p 			= 		new TFile("20210408_1555-h_p-w_zcut.root");

	// for(int n = 0; n < v_part.size(); n++)
	// {
	// 	stringstream 			 	   h_title;
	// 	h_title  <<  "h_" <<   		v_part.at(n).p_dno<< "-" << n;
	// 	string h_tit_str  = 	 h_title.str();

	// 	TH1F* 	h = (TH1F*) tf_p->Get(h_tit_str.c_str()); // getting all detectors

	// 	for(int i =0; i < 9; i++) 		//Determining ROI. Q_bin_ID is the bin number where Q lies. min_bin_ID is the bin number of -3 sigma from Q. 
	// 	{
	// 		Q_bin_ID[i] 	= ceil(Q[i]/2);  //the assumption is that the binning is 2kev/bin (this might change!)
	// 		resolution[i] 	= Get_Resolution(v_part.at(n).p_dno, Q[i], v_part.at(n));

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
	// 	cout << " =========== Detector number =  " << v_part.at(n).p_dno << " Partition:" << n <<" ==============" << endl;

	// 	// double exposure 	= 	Get_Total_Exposure(v_part.at(n).p_dno, v_det);
	// 	double exposure 	= 	v_part.at(n).p_exp;

	// 	for( int i = 0 ; i < 9 ; i++)
	// 	{
	// 		t_half[i] = ln2 * N_A * a[i] * exposure / ( W[i] * sensitivities[i] );

	// 		cout << "T_half : " << t_half[i] << endl; 

	// 	}

	// 	cout << " ===================================================================" << endl;

	// 	detector_ID += 1;

	// }



	// // CALCULATING T1/2 FROM PARTITION HISTOGRAMS	 		 	/////////////////////

	// // ----------------------------------------- 		/////////////////////

 // 	////////////////////////////////////// Detector Parameters ////////////////////	

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


 // 	////////////////////////////////////// Detector Parameters ////////////////////	


 // 	////////////////////////////////////// T1/2 VS k*sigma     ////////////////////
	// TH1F* 		h[total_partitions];

	// MPFeldman_Cousins* 		mpfc 				= 		new MPFeldman_Cousins(gaus_cutoff, 0.9);
	// // TFile* 					tf_1 				= 		new TFile(h_p_outf);

	// // vector<string> 			root_file_path 	= 		ReadFiles();
	// // TChain* 			 	T_Chain_det		=  		Make_TChain(root_file_path, runs);
	// // vector<CO_detector>  	v_det 			= 		Fill_CO_detector(T_Chain_det);


	// stringstream 			 	   h_title;
	// h_title  <<  "h_2013" 			;
	// string h_tit_str  = 	 h_title.str();


	// for(int p = 0; p < total_partitions; p++)
	// {
	// 	stringstream 			 	   h_title;			//names for histograms		
	// 	h_title  <<  "h_" <<   v_part.at(p).p_dno << "-" << p;
	// 	string h_tit_str  = 	 h_title.str();

	// 	h[p] = (TH1F*) h_p[p];//tf_1->Get(h_tit_str.c_str()); // getting all detectors
	// }
	// double exposure = 0;

	// cout << " +++++ OBTAINING ENTRIES IN ROI and PARTITION RESOLUTIONS ++++++" << endl;
	// for( unsigned int k = 0; k < k_tot; k ++ )
	// {
	// 	k_sig = double(k * k_stp + k_stp);
	// 	// cout << " ++++++++ k = " << k * k_stp + k_stp << " +++++++++++" << endl; 
	// 	for(unsigned int p = 0; p < total_partitions; p++)
	// 	{
	// 		for(unsigned int q = 0; q < 9; q++)
	// 		{
	// 			Q_bin_ID[q] 	= ceil(Q[q]/2.0);  //the assumption is that the binning is 2kev/bin (this might change!)
	// 			resolution[q]   = v_part.at(p).Get_p_res( Q[q] ) ;

	// 			min_E = ( Q[q] - k_sig * resolution[q] ) ; 
	// 			max_E = ( Q[q] + k_sig * resolution[q] ) ; 

	// 			min_bin_ID[q]	= ceil ( ( min_E )/ 2.0);  	
	// 			max_bin_ID[q]   = ceil ( ( max_E )/ 2.0);
			
	// 			entries[q][k] = 0;

	// 			if(min_bin_ID[q] < 0) // || max_bin_ID[q] > 500)
	// 			{
	// 				min_bin_ID[q] = 0;
	// 			}
	// 			else if( max_bin_ID[q] > 5000 )
	// 			{
	// 				max_bin_ID[q] = 5000;
	// 			}
	// 			for( int b = min_bin_ID[q]; b <= max_bin_ID[q]; b++ )
	// 			{
	// 				if(b == min_bin_ID[q]) //proportion of the amount of entries[in[k] relation to position of -k*sigma within bin
	// 				{
	// 					double ratio_factor =  (h[p]->GetXaxis()->GetBinCenter(b) + h[p]->GetXaxis()->GetBinWidth(b) / 2 - min_E) / ( h[p]->GetXaxis()->GetBinWidth(b) );
	// 					double entries_full = h[p]->GetBinContent(b);

	// 					entries[q][k] += ratio_factor * entries_full;

	// 				}
	// 				else if(b == max_bin_ID[q]) //proportion of the amount of entries[in[k] relation to position of +k*sigma within bin
	// 				{
	// 					double ratio_factor =  ( max_E - h[p]->GetXaxis()->GetBinCenter(b) + h[p]->GetXaxis()->GetBinWidth(b) / 2 ) / ( h[p]->GetXaxis()->GetBinWidth(b) );
	// 					double entries_full = h[p]->GetBinContent(b);

	// 					entries[q][k] += ratio_factor * entries_full;
	// 				}
	// 				else
	// 				{
	// 					entries[q][k] += h[p]->GetBinContent(b); 				//sum of all events in ROI
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	// cout << " +++++ OBTAINING ENTRIES IN ROI and PARTITION RESOLUTIONS ++++++" << endl;

	// cout << " START THALF " << endl;

	// for(unsigned int p = 0; p < v_part.size(); p++)
	// {
	// 	exposure 	+= 	v_part.at(p).p_exp;
	// }
	// cout << " TOTAL EXPOSURE IS: " << exposure << endl;

	// 		// double exposure 	= 	Get_Total_Exposure( d_num[d] , v_det );
	// for(unsigned int k = 0; k < k_tot; k ++ )
	// {
	// 	k_sig = double(k * k_stp );
	// 	for( int q = 0 ; q < 9 ; q++)
	// 	{
	// 		cout << "+++++++++++++++++++++++++++" << endl;
			
	// 		if(entries[q][k] > gaus_cutoff)
	// 		{
	// 			sensitivities[q][k] = k_sig * sqrt(entries[q][k]);
	// 		}
	// 		else
	// 		{
	// 			sensitivities[q][k] = mpfc->get_sensitivity(entries[q][k]);
	// 		}
	// 		t_half[q][k] = ln2 * N_A * a[q] * exposure * erf(k_sig/sqrt(2)) / ( W[q] * sensitivities[q][k] );

	// 		cout << "T_half : " << t_half[q][k] << endl; 
	// 		cout << "+++++++++++++++++++++++++++" << endl << endl;

	// 	}
	// }
	
	// cout << " END THALF " << endl;

	// TGraph* 	t_g[9];
	// TCanvas* 	t_c[9];

	// for(int g = 0; g < 9; g ++)
	// {
	// 	stringstream 			 	   t_g_title;
	// 	t_g_title  <<  "Q_" 			<< g;
	// 	string t_g_tit_str  = 	 t_g_title.str();

	// 	t_g[g] = new TGraph();
	// 	t_c[g] = new TCanvas(t_g_tit_str.c_str(), t_g_tit_str.c_str(), 1000, 600);
	// }

	// for(unsigned int k = 0; k < k_tot; k ++ )
	// {
	// 	k_sig = double(k * k_stp );
	// 	for(int q = 0; q < 9 ; q++)
	// 	{
	// 		t_g[q]->SetPoint(k, k_sig , t_half[q][k] );
	// 	}
	// }

	// for(int g = 0; g < 9 ; g++ )
	// {
	// 	stringstream 			 	   t_g_title;
	// 	t_g_title  << " Half-life Limit of " << Iso[g] 			<< " for k*sigma ROI at CL = 90%; k-multiple of sigma; Calculated Half-Life" ;
	// 	string t_g_tit_str  = 	 t_g_title.str();

 // 		t_g[g]->SetMarkerColor(g+1);
	// 	t_g[g]->SetMarkerStyle(g+20);
	// 	t_g[g]->SetTitle(t_g_tit_str.c_str());
	// 	t_c[g]->cd();
	// 	t_g[g]->Draw("apl");
	// }

	// TFile*  	g_tf = new TFile(g_t_outf, "RECREATE");
	// for (int i = 0; i < 9; i++)
	// {
	// 	t_g[i]->Write();
	// }
	// delete g_tf;

	// double TG_k[k_tot];

	// for(unsigned int k = 0; k < k_tot; k ++ )
	// {
	// 	TG_k[k] = double(k * k_stp );
	// }

	// TGraph* TG_ent = new TGraph(k_tot, TG_k, entries[8]);
	// TGraph* TG_sen = new TGraph(k_tot, TG_k, sensitivities[8]);
	// TCanvas* TG_c  = new TCanvas("TG_c","The total Background events and respective Calculates Sensitities as a function of k*sigma; k-multiple",1000,600);
	// TMultiGraph* TG_mg = new TMultiGraph();

	// TG_ent->SetMarkerStyle(21);
	// TG_ent->SetMarkerColor(2);
	// TG_ent->SetName("TG_ent");

	// TG_sen->SetMarkerStyle(22);
	// TG_sen->SetMarkerColor(3);
	// TG_sen->SetName("TG_sen");

	// auto legend = new TLegend(0.1,0.7,0.48,0.9);
	// legend->SetHeader("Legend","C_fit"); 
	// legend->AddEntry("TG_ent","Total Background Events","lp");
	// legend->AddEntry("TG_sen","Calculated Sensitivity","lp");

	// TG_c->cd();
	// TG_mg->Add(TG_ent);
	// TG_mg->Add(TG_sen);
	// TG_mg->Draw("ap");
	// legend->Draw();

 	//////////////////////////////////////// T1/2 VS k*sigma     ////////////////////	

