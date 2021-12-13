#include "LJAnalysis.h"
#include "global.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <TF1.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TCanvas.h>

LJAnalysis::~LJAnalysis() {
	delete _tree;
}

void LJAnalysis::rotateBySlope(float & x, float & z, const float slope) {
	//This is a simple rotation matrix
	//slope = tan = sin/cos, therefore cos = 1/sqrt(1+tan^2) and sin = tan/sqrt(1+tan^2), simply replace tan by slope
	float denom = sqrt(1 + slope * slope);
	float xout = (x - z * slope) / denom;
	float zout = (x * slope + z) / denom;
	x = xout;
	z = zout;

}

void LJAnalysis::smooth(float* arrIn, float* arrOut, const int size,  int smooth) {
	if (smooth % 2 == 0) smooth++; //Smooth needs to be odd
	if (smooth < 3) smooth = 3; //smaller than 3 makes no sense
	if (smooth > size) smooth = size; //makes no sense
	for (int i = 0; i < size; i++) {
		float sum = 0;
		int n = 0;
		for (int j = i - (smooth - 1) / 2; j <= i + (smooth - 1) / 2; j++) {
			int myJ = j;
			if (j < 0) myJ = -j;
			if (j > size) myJ = size - (j - size + 1);
			if (arrIn[myJ] > -50) {
				sum += arrIn[myJ];
				n++;
			}
			//std::cout << "Calc Smooth: " << i << " " << myJ << " " << n << std::endl;
		}
		if (n > 0) arrOut[i] = sum / n; else arrOut[i] = -999;
	}
}





LJAnalysis::LJAnalysis(const DataKind kind) {
	_kind = kind;

	_tree = new TTree("tree", "tree");
	_deliminator = ";";
	_frequency = 16000; //Hz;
	_xResolution = +0.04; //mm/Dots .. doing so 0 is in the center point of the sensor
	_xOffset = -8.0 + _xResolution / 2.;//mm

	_hTiltNom1 = new TH1F("hTiltNom1", "Tilt Nominator (A)", 400, -8, 8); //This is for raw-data from 8-8
	_hTiltDenom1 = (TH1F*)_hTiltNom1->Clone("hTiltDenom1");
	_hTiltDenom1->SetTitle("Tilt Denominator (A)");

	if (_kind = kTwoLJ7060_800_16) {
		_hTiltNom2 = new TH1F("hTiltNom2", "Tilt Nominator (B)", 400, -8, 8); //This is for raw-data from 8-8
		_hTiltDenom2 = (TH1F*)_hTiltNom2->Clone("hTiltDenom2");
		_hTiltDenom2->SetTitle("Tilt Denominator (B)");
	}
	_hRangeNom1 = new TH2F("hRangeNom1", "Range Nominator (A)", 200, 0, 0.20, 450, -2, 16);
	_hRangeNom1->GetXaxis()->SetTitle("time [s]");
	_hRangeNom1->GetYaxis()->SetTitle("Pos on Sonotrode [mm]");
	
	_hRangeDenom1 = (TH2F*)_hRangeNom1->Clone("hRangeDenom1");
	_hRangeDenom1->SetTitle("Amplitude Denominator (A)");

	_hRangeNomSmooth1 = (TH2F*)_hRangeNom1->Clone("hRangeNomSmooth1");
	_hRangeNomSmooth1->SetTitle("Amplitude Denominator Smooth (A)");

	_hRangeDenomSmooth1 = (TH2F*)_hRangeNom1->Clone("hRangeDenomSmooth1");
	_hRangeDenomSmooth1->SetTitle("Amplitude Denominator Smooth (A)");

	_hRangeCrossNom1 = new TH1F("hRangeCrossNom1", "CrossRange Nominator (A)", 450, -2, 16);
	_hRangeCrossNom1->GetXaxis()->SetTitle("Pos on Sonotrode [mm]");
	_hRangeCrossNom1->GetYaxis()->SetTitle("Mean Range P2P [mm]");

	_hRangeCrossDenom1 = (TH1F*)_hRangeCrossNom1->Clone("hRangeCrossDenom1");
	_hRangeCrossDenom1->SetTitle("Cross Amplitude Denom (A)");

	_hRangeCrossNomSmooth1 = (TH1F*)_hRangeCrossNom1->Clone("hRangeCrossDenomSmooth1");
	_hRangeCrossNomSmooth1->SetTitle("Cross Amplitude Denom (A)");

	_hRangeCrossDenomSmooth1 = (TH1F*)_hRangeCrossNom1->Clone("hRangeCrossDenomSmooth1");
	_hRangeCrossDenomSmooth1->SetTitle("Cross Amplitude Denom (A)");

	_grAmplPit1 = new TGraph();
	//_grAmplPit->GetYaxis()->SetTitle("Amplitude P2P at Max[#mum]");
	//_grAmplPit->GetXaxis()->SetTitle("time [s]");
	_grAmplPit1->SetTitle("Max Amplitude at Pit (A)");
	
	_grNullPos1 = new TGraph();
	_grNullPos1->SetTitle("Position of NP (A)");
	//_grNullPos->GetXaxis()->SetTitle("time [s]");
	//_grNullPos->GetYaxis()->SetTitle("Position on Sonotrode [mm]");
	
	_grAmplPos1 = new TGraph();
	_grAmplPos1->SetTitle("Pos of Max (A)");
	//_grAmplPos->GetXaxis()->SetTitle("time [s]");
	//_grAmplPos->GetYaxis()->SetTitle("Position on Sonotrode [mm]");
	
	if (_kind = kTwoLJ7060_800_16) {
		_hRangeNom2 = new TH2F("hRangeNom2", "Range Nominator (B)", 200, 0, 0.20, 450, -16, 2);
		_hRangeNom2->GetXaxis()->SetTitle("time [s]");
		_hRangeNom2->GetYaxis()->SetTitle("Pos on Sonotrode [mm]");

		_hRangeDenom2 = (TH2F*)_hRangeNom2->Clone("hRangeDenom2");
		_hRangeDenom2->SetTitle("Amplitude Denominator (B)");

		_hRangeNomSmooth2 = (TH2F*)_hRangeNom2->Clone("hRangeNomSmooth2");
		_hRangeNomSmooth2->SetTitle("Amplitude Denominator Smooth (B)");

		_hRangeDenomSmooth2 = (TH2F*)_hRangeNom2->Clone("hRangeDenomSmooth2");
		_hRangeDenomSmooth2->SetTitle("Amplitude Denominator Smooth (B)");

		_hRangeCrossNom2 = new TH2F("hRangeCrossNom2", "CrossRange Nominator (B)", 450, -16, 2);
		_hRangeCrossNom2->GetXaxis()->SetTitle("Pos on Sonotrode [mm]");
		_hRangeCrossNom2->GetYaxis()->SetTitle("Mean Range P2P [mm]");

		_hRangeCrossDenom2 = (TH2F*)_hRangeCrossNom2->Clone("hRangeCrossDenom2");
		_hRangeCrossDenom2->SetTitle("Cross Amplitude Denom (B)");

		_hRangeCrossNomSmooth2 = (TH2F*)_hRangeCrossNom2->Clone("hRangeCrossDenomSmooth2");
		_hRangeCrossNomSmooth2->SetTitle("Cross Amplitude Denom (B)");

		_hRangeCrossDenomSmooth2 = (TH2F*)_hRangeCrossNom2->Clone("hRangeCrossDenomSmooth2");
		_hRangeCrossDenomSmooth2->SetTitle("Cross Amplitude Denom (B)");

		_grAmplPit2 = new TGraph();
		//_grAmplPit->GetYaxis()->SetTitle("Amplitude P2P at Max[#mum]");
		//_grAmplPit->GetXaxis()->SetTitle("time [s]");
		_grAmplPit2->SetTitle("Max Amplitude at Pit (B)");

		_grNullPos2 = new TGraph();
		_grNullPos2->SetTitle("Position of NP (B)");
		//_grNullPos->GetXaxis()->SetTitle("time [s]");
		//_grNullPos->GetYaxis()->SetTitle("Position on Sonotrode [mm]");

		_grAmplPos2 = new TGraph();
		_grAmplPos2->SetTitle("Pos of Max (B)");
		//_grAmplPos->GetXaxis()->SetTitle("time [s]");
		//_grAmplPos->GetYaxis()->SetTitle("Position on Sonotrode [mm]");
	}

	std::cout << "Constructeed LJAnalysis " << std::endl;
}

void LJAnalysis::print(std::string opt) {
	//Add the filenameextension

	

	TCanvas* c = new TCanvas("c", "CANVAS", 1200, 800);
	if (opt == "range") {
		std::string newfilename = std::string(_name) + "-col.gif";
		_hRange->SetTitle(_name.c_str());
		_hRange->GetZaxis()->SetRangeUser(0, 0.2);
		_hRange->Draw("colz");
		_grNullPos->Draw("lp");
		_grAmplPos->Draw("lp");
		c->SaveAs(newfilename.c_str());
	} else if (opt == "amplitude") {
		c->SetGrid(0, 1);
		_grAmplPit->SetTitle(_name.c_str());
		std::string newfilename = std::string(_name) + "-Ampl.gif";
		_grAmplPit->Draw("alp");
		_grAmplPit->GetYaxis()->SetRangeUser(0, 0.2);
		c->SaveAs(newfilename.c_str());
	}
	else if (opt == "cross") {
		c->SetGrid(1, 1);
		_hRangeCross->SetTitle(_name.c_str());
		std::string newfilename = std::string(_name) + "-cross.gif";
		_hRangeCross->GetYaxis()->SetRangeUser(0, 0.2);
		_hRangeCross->Draw("hist");
		c->SaveAs(newfilename.c_str());
	}

	if (opt == "rangesm") {
		std::string newfilename = std::string(_name) + "-colsm.gif";
		_hRangeSmooth->SetTitle(_name.c_str());
		_hRangeSmooth->GetZaxis()->SetRangeUser(0, 0.2);
		_hRangeSmooth->Draw("colz");
		_grNullPos->Draw("lp");
		_grAmplPos->Draw("lp");
		c->SaveAs(newfilename.c_str());
	}
	else if (opt == "crosssm") {
		_hRangeCrossSmooth->SetTitle(_name.c_str());
		std::string newfilename = std::string(_name) + "-crosssm.gif";
		_hRangeCrossSmooth->Draw("");
		c->SaveAs(newfilename.c_str());
	}
	delete c;


}

Int_t LJAnalysis::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!_tree) return 0;
	return _tree->GetEntry(entry);
}
Long64_t LJAnalysis::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!_tree) return -5;
	Long64_t centry = _tree->LoadTree(entry);
	if (centry < 0) return centry;
	if (_tree->GetTreeNumber() != fCurrent) {
		fCurrent = _tree->GetTreeNumber();
		//Notify();
	}
	return centry;
}

void LJAnalysis::analysis() {
	if (!_tree) {
		std::cerr << "No TTree... aborting" << std::endl;
		exit(1);
	}
	/* Analysis runs on TTree level */

	/* Conneting the TTree with the variables */
	
	Float_t         rawX[400];
	Float_t         rawZ1[400];
	Float_t         rawZ2[400];
	Float_t         rawTime;
	Float_t         rawSensor;


	
	TBranch* b_raw_x;   //!
	TBranch* b_raw_z1;   //!
	TBranch* b_raw_z2;   //!
	TBranch* b_raw_time;   //!
	TBranch* b_raw_sensor;   //!

	// Set branch addresses and branch pointers

	
	_tree->SetBranchAddress("raw.x", rawX, &b_raw_x);
	_tree->SetBranchAddress("raw.z1", rawZ1, &b_raw_z1);
	if (_kind = kTwoLJ7060_800_16) {
		_tree->SetBranchAddress("raw.z2", rawZ2, &b_raw_z2);
	}
	_tree->SetBranchAddress("raw.time", &rawTime, &b_raw_time);
	_tree->SetBranchAddress("raw.sensor", &rawSensor, &b_raw_sensor);

	// New TTree to store further analysed variables
	TTree* anatree = new TTree("anatree", "anatree");
	anatree->Branch("raw.x", rawX);
	if (_kind = kTwoLJ7060_800_16) {
		anatree->Branch("raw.z2", rawZ2);
	}
	anatree->Branch("raw.z1", rawZ1);
	anatree->Branch("raw.time", &rawTime);
	anatree->Branch("raw.sensor", &rawSensor);


	static const int nDots = 400;

	float z1[nDots];
	float z2[nDots];
	float x[nDots];
	float deltaZ1[nDots];
	float rangeZ1[nDots];
	float smoothZ1[nDots];
	float smoothDeltaZ1[nDots];
	float smoothRangeZ1[nDots];
	float deltaZ2[nDots];
	float rangeZ2[nDots];
	float smoothZ2[nDots];
	float smoothDeltaZ2[nDots];
	float smoothRangeZ2[nDots];
	anatree->Branch("ana.z1", &z,"ana.z1[400]/F");
	if (_kind = kTwoLJ7060_800_16) {
		anatree->Branch("ana.z2", &z, "ana.z2[400]/F");
	}
	anatree->Branch("ana.x", &x,"ana.x[400]/F");
	anatree->Branch("ana.deltaZ1", &deltaZ1,"ana.deltaZ1[400]/F");
	anatree->Branch("ana.rangeZ1", &rangeZ1,"ana.rangeZ1[400]/F");
	anatree->Branch("ana.smoothZ1", &smoothZ1, "ana.smoothZ1[400]/F");
	anatree->Branch("ana.smoothDeltaZ1", &smoothDeltaZ1, "ana.smoothDeltaZ1[400]/F");
	anatree->Branch("ana.smoothRangeZ1", &smoothRangeZ1, "ana.smoothRangeZ1[400]/F");

	anatree->Branch("ana.deltaZ2", &deltaZ2, "ana.deltaZ2[400]/F");
	anatree->Branch("ana.rangeZ2", &rangeZ2, "ana.rangeZ2[400]/F");
	anatree->Branch("ana.smoothZ2", &smoothZ2, "ana.smoothZ2[400]/F");
	anatree->Branch("ana.smoothDeltaZ2", &smoothDeltaZ2, "ana.smoothDeltaZ2[400]/F");
	anatree->Branch("ana.smoothRangeZ2", &smoothRangeZ2, "ana.smoothRangeZ2[400]/F");


	static const int nElements = 20;
	float lastElements1[nDots][nElements];
	float lastElementsSmooth1[nDots][nElements];
	float lastElements2[nDots][nElements];
	float lastElementsSmooth2[nDots][nElements];




	Long64_t nentries = _tree->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry = 0; jentry < nentries;jentry++) {
		if (jentry % 1000 == 0) {
			std::cout << "Analysing entry " << jentry << std::endl;
		}
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = _tree->GetEntry(jentry);   nbytes += nb;
		for (int col = 0; col < nDots; col++) {
			if (rawZ1[col] < -50) {
				z1[col] = -999;
				rangeZ[col] = -999;
			}
			else {

				//Changing the orientation, by shift and rotation....in the end we will look perendicular to the sonotrode
				//first copy
				x[col] = rawX[col];
				z1[col] = rawZ[col];


				//Shift, so that z(0) = 0!
				z1[col] -= _offset1;

				//Rotate by the negative of the estimated slope
				rotateBySlope(x[col], z1[col], -_slope1);

				//now more the pitpos to 0;
				//Recalculate the pitpos into the coordinate system, too
				float pitPosZ1 = _slope1 * _pitpos1;
				float pitPosX1 = _pitpos1; //Keep the porignal one
				rotateBySlope(pitPosX, pitPosZ, -_slope1);
				//std::cout << "Recalculated pitpos is " << pitPosX << " " << pitPosZ << std::endl;

				x[col] -= pitPosX1;

				//Move in such a way, that the sonotrode's pit is at 0

				//----------------------------------

				lastElements1[col][jentry % nElements] = z[col];
				if (jentry > nElements) { //Calculate the range
					float thisMin = 9999;
					float thisMax = -9999;
					for (int nLast = 0; nLast < nElements; nLast++) {
						if (lastElements1[col][nLast] < thisMin) {
							thisMin = lastElements1[col][nLast];
						}
						if (lastElements1[col][nLast] > thisMax) {
							thisMax = lastElements1[col][nLast];
						}
					}
					if (thisMax < 50 && thisMin > -50 && thisMax > -50 && thisMin < 50) {
						rangeZ1[col] = thisMax - thisMin; //This is the range
						//To calculate the mean filling the sum and the n, x is inverted!!!
						_hRangeNom1->Fill(rawTime, -x[col], rangeZ1[col]);
						_hRangeDenom1->Fill(rawTime, -x[col]);
						_hRangeCrossNom1->Fill(-x[col], rangeZ1[col]);
						_hRangeCrossDenom-1>Fill(-x[col]);
					}
					else rangeZ1[col] = -999;
					
				}
				else {
					rangeZ1[col] = -999;
				}




			}

		}

		smooth(z1, smoothZ1, 400, 20);
		/*
		// Analysis for smoothed data
		//--------------------------------------
		for (int col = 0; col < nDots; col++) {
			lastElementsSmooth[col][jentry % nElements] = z[col];
			if (jentry > nElements) { //Calculate the range
				float thisMin = 9999;
				float thisMax = -9999;
				for (int nLast = 0; nLast < nElements; nLast++) {
					if (lastElementsSmooth[col][nLast] < thisMin) {
						thisMin = lastElementsSmooth[col][nLast];
					}
					if (lastElementsSmooth[col][nLast] > thisMax) {
						thisMax = lastElementsSmooth[col][nLast];
					}
				}
				if (thisMax < 50 && thisMin > -50 && thisMax > -50 && thisMin < 50) {
					smoothRangeZ[col] = thisMax - thisMin; //This is the range
					//To calculate the mean filling the sum and the n, x is inverted!!!
					_hRangeNomSmooth->Fill(rawTime, -x[col], smoothRangeZ[col]);
					_hRangeDenomSmooth->Fill(rawTime, -x[col]);
					_hRangeCrossNomSmooth->Fill(-x[col], rangeZ[col]);
					_hRangeCrossDenomSmooth->Fill(-x[col]);
				}
				else smoothRangeZ[col] = -999;

			}
			else {
				smoothRangeZ[col] = -999;
			}
		}
		*/
		
		//Filling all connected variables into the new anatree
		anatree->Fill();
	} //This is the end of the treelooper

	//Calculation of the mean Range
	_hRange1 = (TH2F*)_hRangeNom1->Clone("hRange1");
	//_hAmpl->Sumw2();
	_hRange1->SetTitle("Mean Range (A)");
	_hRange1->Divide(_hRangeDenom1);
	/*
	//Calculation of the mean Range Smoothed
	_hRangeSmooth = (TH2F*)_hRangeNomSmooth->Clone("hRangeSmooth");
	//_hAmpl->Sumw2();
	_hRangeSmooth->SetTitle("Mean Range Smoothed");
	_hRangeSmooth->Divide(_hRangeDenomSmooth);
	*/
	//Calculation of the mean Range
	_hRangeCross1 = (TH1F*)_hRangeCrossNom1->Clone("hRangeCross1");
	//_hAmpl->Sumw2();
	_hRangeCross1->SetTitle("Mean Range Cross (A)");
	_hRangeCross1->Divide(_hRangeCrossDenom1);
	/*
	//Calculation of the mean Range Smoothed
	_hRangeCrossSmooth = (TH1F*)_hRangeCrossNomSmooth->Clone("hRangeCrossSmoothed");
	//_hAmpl->Sumw2();
	_hRangeCrossSmooth->SetTitle("Mean Range Cross Smoothed");
	_hRangeCrossSmooth->Divide(_hRangeDenomSmooth);
	*/
	


	//Calculation of the NP
	//TF1* fNP = new TF1("fNP", "[0]+[1]*fabs(sin([2]*x +[3]))+ [4]*x", 4, 16);
	TF1* fNP1 = new TF1("fNP1", "[0] + fabs([1]*x+[2])", 4, 16);
	
	fNP1->SetParameter(0, 0);
	fNP1->SetParameter(1, 0);
	fNP1->SetParameter(2, 0.05);
	fNP1->SetParLimits(0, 0, 0.1);
	fNP1->SetParLimits(1, -0.1, 0);
	fNP1->SetParLimits(2, 0, 0.1);
	/*
	fNP->SetParameter(0, 0);
	fNP->SetParameter(1, 1);
	fNP->SetParameter(2, +6.8e-3);
	fNP->SetParameter(3, 10);
	fNP->SetParameter(4, -1.6e-3);

	fNP->SetParLimits(0, 0, 0.5);
	fNP->SetParLimits(1, 0, 10);
	fNP->SetParLimits(2, -1e-2, 1e-2);
	fNP->SetParLimits(3, 5, 15);
	fNP->SetParLimits(4, -0.0001, -0.001);
	*/
	
	int nbins = _hRange1->GetXaxis()->GetNbins();
	for (int i = 1; i <= nbins; i++) { //this is one time slice
		char nameproj[64];
		sprintf(nameproj, "_projy%d", i);
		TH1F* myProj = (TH1F*)_hRange1->ProjectionY(nameproj, i, i);
		float myMaxAmpl = -99;
		float maxPos = -1;
		for (int j = 1; j <= myProj->GetXaxis()->GetNbins(); j++) { //this is one x-pos
			//if (myProj->GetBinContent(j) < 1) continue;
			if (myProj->GetBinContent(j) > myMaxAmpl) {
				myMaxAmpl = myProj->GetBinContent(j);
				maxPos = myProj->GetBinCenter(j);
			}
			if (myProj->GetBinCenter(j) > 5) break;
		}
		if (myMaxAmpl > -99) {
			_grAmplPit->SetPoint(_grAmplPit->GetN(), _hRange->GetXaxis()->GetBinWidth(0) * i, myMaxAmpl);
			_grAmplPos->SetPoint(_grAmplPos->GetN(), _hRange->GetXaxis()->GetBinWidth(0) * i, maxPos);
			//std::cout << "Max Ampl: " << _hAmpl->GetXaxis()->GetBinWidth(0)*i << " " << myMaxAmpl << " " << _grAmplPit->GetN() << std::endl;
		}
		
		myProj->Fit("fNP", "RQ");
	
	
		float NPPos = fNP->GetMinimumX(4, 16);
		//std::cout << "NP: " << NPPos << std::endl;
		if (NPPos > 4 && NPPos < 16) {
			_grNullPos->SetPoint(_grNullPos->GetN(), _hRange->GetXaxis()->GetBinWidth(0) * i, NPPos);
		}
	}
	_grNullPos->GetXaxis()->SetTitle("time [s]");
	_grNullPos->GetYaxis()->SetTitle("Position on Sonotrode [mm]");
	_grNullPos->SetName("NullPos");
	_grNullPos->SetLineWidth(3);
	_grNullPos->Write();

	_grAmplPit->GetXaxis()->SetTitle("time [s]");
	_grAmplPit->GetYaxis()->SetTitle("Amplitude P2P [mm]");
	_grAmplPit->SetLineWidth(3);
	_grAmplPit->SetLineColor(kRed);
	_grAmplPit->SetName("Amplitude");
	_grAmplPit->Write();

	_grAmplPos->GetXaxis()->SetTitle("time [s]");
	_grAmplPos->GetYaxis()->SetTitle("Position on Sonotrode [mm]");
	_grAmplPos->SetName("Max-Position");
	_grAmplPos->SetLineWidth(3);
	_grAmplPos->SetLineColor(kRed);
	_grAmplPos->Write();




}

TTree* LJAnalysis::convert(std::string filename) {
	//std::cout << "convert of " << filename << std::endl;

	/* convert runs on file-io level, but produces the TTree */

	

	static const int nDots = 400;
	float x[nDots];
	float z1[nDots];
	float z2[nDots];
	int sensor;
	float time;
	int nDataSets = 0;


	for (int i = 0; i < nDots; i++) {
		x[i] = -999;
		z1[i] = -999;
		z2[i] = -999;
	
	}


	

	
	_tree->Branch("raw.x", &x,"raw.x[400]/F");
	_tree->Branch("raw.z1", &z1,"raw.z1[400]/F");
	if (_kind == kTwoLJ7060_800_16) {
		_tree->Branch("raw.z2", &z2, "raw.z2[400]/F");
	}
	_tree->Branch("raw.time", &time);
	_tree->Branch("raw.sensor", &sensor);
	
	std::ifstream infile(filename.c_str());
	int linenum = -1;
	std::string sline;
	std::vector <std::string> dataVec;
	bool isFirstTime = true;
	if (!infile.is_open()) {
		std::cerr << "ERROR: Unable to open input file '" << filename << "'!" << std::endl;
		exit(1);
	} else {
		//std::cout << "Starting conversion!" << std::endl;
		while (!infile.eof()) {
			linenum++;
			if (linenum % 1000 == 0) {
				std::cout << "Converting line " << linenum << " of " << filename << std::endl;
			}
			getline(infile, sline);
			if (linenum < 5) continue;
			if (sline.length() < 1) continue;
			if (linenum > 4000) {
				std::cout << "Stopping conversion after " << linenum << " lines!" << std::endl;
				break;
			}
			
			//Replaceing all "," by "."
			std::replace(sline.begin(), sline.end(), ',', '.');
			dataVec = split(sline, _deliminator.at(0));
				
			if (isFirstTime) {
				isFirstTime = false;
				std::cout << "Found line with " << dataVec.size() << " elements!" << std::endl;
			}
			time = linenum / _frequency;
			for (int i = 0; i < dataVec.size(); i++) {
				if (i < 400) { //These data belong to sensor A

					float myZ = atof(dataVec.at(i).c_str());
					float myX = i * _xResolution + _xOffset;
					x[i] = myX;
					sensor = 1;
					if (myZ > -50 && myZ < 50) { //Convert only data in valid range
	
						z1[i] = myZ;

						_hTiltNom1->Fill(myX, myZ);
						_hTiltDenom1->Fill(myX);
					}
					else {
						//x[i] = -999;
						z1[i] = -999;
					}
					
				} else { //These data belong to sensor B
					float myZ = atof(dataVec.at(i).c_str());
					float myX = (i-400) * _xResolution + _xOffset;
					x[i-400] = myX;
					sensor = 2;
					if (myZ > -50 && myZ < 50) { //Convert only data in valid range
						z2[i-400] = myZ;

						_hTiltNom2->Fill(myX, myZ);
						_hTiltDenom2->Fill(myX);
					}
					else {
						z2[i-400] = -999;
					}
					
				}
			}
			
			
			_tree->Fill();
			nDataSets++;
			
		}

		infile.close();
		std::cout << nDataSets << " datasets from file " << filename << " analysed" << std::endl;
		// SENSOR A
		//------------
		//Mean of the z-values (per x) is the sum of the z-values (per x) devided by n
		_hTilt1 = (TH1F*)_hTiltNom1->Clone("hTilt1");
		//_hTilt->Sumw2();
		_hTilt1->SetTitle("Mean Tilt Sensor A");
		_hTilt1->Divide(_hTiltDenom1);

		TF1* fTilt1 = new TF1("fTilt1", "pol1", -8, 8);
		_hTilt1->Fit("fTilt1");
		_offset1 = fTilt1->GetParameter(0);
		_slope1 = fTilt1->GetParameter(1);
		std::cout << "Tilt (A) was determined to " << _slope1 << " (+/- " << fTilt1->GetParError(1) << ")" << std::endl;
		std::cout << "Distance (A) was determined to " << _offset1 << " (+/- " << fTilt1->GetParError(0) << ")" << std::endl;

		//determine the pit position by looking for the last entries of the histogram (small entries < 1% will be ignored)
		Double_t q[1];
		Double_t p[1] = { 1-1e-4 };
		_hTiltDenom1->GetQuantiles(1, q, p);
		std::cout << "Pit (A) position located at " << q[0] << std::endl;
		_pitpos1 = q[0];



		// SENSOR B
		//------------
		if (_kind == kTwoLJ7060_800_16) {
			//Mean of the z-values (per x) is the sum of the z-values (per x) devided by n
			_hTilt2 = (TH1F*)_hTiltNom2->Clone("hTilt2");
			//_hTilt->Sumw2();
			_hTilt2->SetTitle("Mean Tilt Sensor A");
			_hTilt2->Divide(_hTiltDenom2);

			TF1* fTilt2 = new TF1("fTilt2", "pol1", -8, 8);
			_hTilt2->Fit("fTilt2");
			_offset2 = fTilt2->GetParameter(0);
			_slope2 = fTilt2->GetParameter(1);
			std::cout << "Tilt (B) was determined to " << _slope2 << " (+/- " << fTilt2->GetParError(1) << ")" << std::endl;
			std::cout << "Distance (B) was determined to " << _offset2 << " (+/- " << fTilt2->GetParError(0) << ")" << std::endl;

			//determine the pit position by looking for the last entries of the histogram (small entries < 1% will be ignored)

			p[0] = { 1e-4 };
			_hTiltDenom2->GetQuantiles(1, q, p);
			std::cout << "Pit (B) position located at " << q[0] << std::endl;
			_pitpos2 = q[0];
		}
	}



	return _tree;


}