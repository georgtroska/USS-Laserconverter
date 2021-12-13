// USS-Laserconverter.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "stdafx.h"
#include <iostream>
#include <TH1F.h>
#include <TCanvas.h>
#include <TSystemDirectory.h>
#include <TFile.h>
#include "LJAnalysis.h"
#include <TStyle.h>
#include "identifyFormat.h"




int main(int argc, const char* argv[]) {
	std::cerr << std::endl;
	std::cerr << "This is the USS Laserconverter converter!" << std::endl;
	std::cerr << "Version 1.0.1, 2021-12-09" << std::endl;
	std::cerr << std::endl;
	std::cerr << "Developed by" << std::endl;
	std::cerr << "\tInfineon Technologies AG, Warstein" << std::endl;
	std::cerr << "\tDr. Georg Troska & Dr. Till Neddermann" << std::endl;
	std::cerr << std::endl;


	gStyle->SetOptStat(0);
	//ArgParser ap;
	//ap.parse(argc, argv);
	//ap.print();

	if (argc > 2) {
		std::cout << "Currently only one or no parameter is possible" << std::endl;
		exit(1);
	}
	std::string dirname;
	if (argc == 1) dirname = std::string("./");
	if (argc == 2) dirname = std::string(argv[2]);

    std::string ext = ".csv";
	TSystemDirectory dir(dirname.c_str(), dirname.c_str());
	TList* files = dir.GetListOfFiles();
	if (files) {

		TSystemFile* file;
		TString fname;
		TIter next(files);
		while ((file = (TSystemFile*)next())) {

			DataKind kind = kUnknown;
			fname = file->GetName();
			std::string myFilename = std::string(fname.Data());
			if (!file->IsDirectory() && fname.EndsWith(ext.c_str())) {

				std::cout << "Analysing: " << fname.Data() << std::endl;
				std::string filename;
				std::string filename_wo_ext = std::string(fname.Data());
				filename_wo_ext.erase(filename_wo_ext.find(".csv", 4));
				if (dirname == "./") { // avoid the ./
					filename = myFilename;
				}
				else {
					filename = dirname + "/" + myFilename;
				}

				//Checking what might be inside??
				kind = identifyFormat(filename);
				
				//Create the root-filename from the csv-filename
				std::string rootfilename = std::string(myFilename); 
				rootfilename.replace(rootfilename.find(".csv"), 4, ".root");
				//Including path
				std::string outfilename;
				if (dirname == "./") { // avoid the ./
					outfilename = std::string(rootfilename);
				}
				else {
					outfilename = std::string(dirname + "/" + rootfilename);
				}
				std::cout << "Outfilename is " << outfilename << std::endl;
				
			
				
				TFile* file = new TFile(outfilename.c_str(), "RECREATE");
				LJAnalysis *ljana = new LJAnalysis(kUnknown);
				ljana->setName(filename_wo_ext);
				TTree *tree = ljana->convert(filename);
				//ljana->analysis();
				file->Write();
				/*
				ljana->print("range");
				//ljana->print("rangesm");
				ljana->print("amplitude");
				ljana->print("cross");
				//ljana->print("crosssm");
				*/

				delete ljana;

			}
		}
	}

	




    std::cout << "Hello World!" << std::endl;
    TH1F* h = new TH1F("h", "Huhu", 10, 0, 10);
    h->Fill(3);
    TCanvas* c1 = new TCanvas("c1", "c1",1200,800);
    h->Draw();
    c1->SaveAs("huhu.gif");
    
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
