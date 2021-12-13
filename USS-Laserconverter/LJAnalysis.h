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
#include "identifyFormat.h"




class LJAnalysis {
protected:
	TTree* _tree;
	DataKind _kind;

	std::string _name;

	float _frequency = 16000; //Hz;
	float _xResolution = 0.04; //mm/Dots
	float _xOffset = _xResolution/2.;//mm
	static const int _nDots = 400; // number of dots per line (this is the number of colums in file)
	std::string _deliminator;

	//Tilt, offset and pitposition is analysed during conversion phase 1 is for sensor A, 2 is for sensor B
	TH1F* _hTiltNom1;
	TH1F* _hTiltDenom1;
	TH1F* _hTilt1;
	float _slope1; //This is the tilt of the sensor
	float _offset1; //this is the distance in the middle point
	float _pitpos1; //This is the position of the pit

	TH1F* _hTiltNom2;
	TH1F* _hTiltDenom2;
	TH1F* _hTilt2;
	float _slope2; //This is the tilt of the sensor
	float _offset2; //this is the distance in the middle point
	float _pitpos2; //This is the position of the pit

	// Variables and methods for TTree analysis
	Int_t GetEntry(Long64_t entry);
	Long64_t LoadTree(Long64_t entry);
	Int_t   fCurrent; //!current Tree number in a TChain

	void rotateBySlope(float &x,  float &z, const float slope);
	void smooth(float* arrIn, float* arrOut, const int size, int smooth = 5);


	//Histograms and stuff like this:
	//SENSOR A
	TH2F * _hRangeNom1;
	TH2F * _hRangeDenom1;
	TH2F * _hRange1;
	TH2F * _hRangeNomSmooth1;
	TH2F * _hRangeDenomSmooth1;
	TH2F * _hRangeSmooth1;

	TH1F* _hRangeCrossNom1;
	TH1F* _hRangeCrossDenom1;
	TH1F* _hRangeCross1;

	TH1F* _hRangeCrossNomSmooth1;
	TH1F* _hRangeCrossDenomSmooth1;
	TH1F* _hRangeCrossSmooth1;

	TGraph * _grAmplPit1;
	TGraph * _grNullPos1;
	TGraph * _grAmplPos1;


	//Sensor B
	TH2F* _hRangeNom2;
	TH2F* _hRangeDenom2;
	TH2F* _hRange2;
	TH2F* _hRangeNomSmooth2;
	TH2F* _hRangeDenomSmooth2;
	TH2F* _hRangeSmooth2;

	TH1F* _hRangeCrossNom2;
	TH1F* _hRangeCrossDenom2;
	TH1F* _hRangeCross2;

	TH1F* _hRangeCrossNomSmooth2;
	TH1F* _hRangeCrossDenomSmooth2;
	TH1F* _hRangeCrossSmooth2;

	TGraph* _grAmplPit2;
	TGraph* _grNullPos2;
	TGraph* _grAmplPos2;

	
	




public:
	LJAnalysis(const DataKind kind);
	void print(std::string opt);//Filename without extension!!!

	virtual ~LJAnalysis();

	TTree * convert(std::string filename);
	void analysis();
	void setName(const std::string name) { _name = name; }

	
};

