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


struct paths
{
	string 		directory;
	string 		year;
	string 		folder;
	string 		file;
};
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

void print_files_from_path() 
{ 
	vector<paths>  	 root_file_path = 	ReadFiles();
	for (int i = 0 ; i < root_file_path.size(); i++)
	{
		char* root_file;

		root_file = str_to_char(root_file_path.at(i).directory + "/" + 
								root_file_path.at(i).year      + "/" + 
								root_file_path.at(i).folder    + "/" + 
								root_file_path.at(i).file );

			
		cout << root_file << endl;
	}
	return;
}