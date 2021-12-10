#pragma once

#include <string>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>

#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>




class LJAnalysis {
protected:
	TTree* _tree;

	std::string _name;

	float _frequency = 16000; //Hz;
	float _xResolution = 0.04; //mm/Dots
	float _xOffset = _xResolution/2.;//mm
	static const int _nDots = 400; // number of dots per line (this is the number of colums in file)
	std::string _deliminator;

	//Tilt, offset and pitposition is analysed during conversion phase
	TH1F* _hTiltNom;
	TH1F* _hTiltDenom;
	TH1F* _hTilt;
	float _slope; //This is the tilt of the sensor
	float _offset; //this is the distance in the middle point
	float _pitpos; //This is the position of the pit

	// Variables and methods for TTree analysis
	Int_t GetEntry(Long64_t entry);
	Long64_t LoadTree(Long64_t entry);
	Int_t   fCurrent; //!current Tree number in a TChain

	void rotateBySlope(float &x,  float &z, const float slope);
	void smooth(float* arrIn, float* arrOut, const int size, int smooth = 5);


	//Histograms and stuff like this:

	TH2F * _hRangeNom;
	TH2F * _hRangeDenom;
	TH2F * _hRange;
	TH2F * _hRangeNomSmooth;
	TH2F * _hRangeDenomSmooth;
	TH2F * _hRangeSmooth;

	TH1F* _hRangeCrossNom;
	TH1F* _hRangeCrossDenom;
	TH1F* _hRangeCross;

	TH1F* _hRangeCrossNomSmooth;
	TH1F* _hRangeCrossDenomSmooth;
	TH1F* _hRangeCrossSmooth;

	TGraph * _grAmplPit;
	TGraph * _grNullPos;
	TGraph * _grAmplPos;

	
	




public:
	LJAnalysis();
	void print(std::string opt);//Filename without extension!!!

	virtual ~LJAnalysis();

	TTree * convert(std::string filename);
	void analysis();
	void setName(const std::string name) { _name = name; }

	
};

