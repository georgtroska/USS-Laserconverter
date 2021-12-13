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





LJAnalysis::LJAnalysis() {
	_tree = new TTree("tree", "tree");
	_deliminator = ";";
	_frequency = 16000; //Hz;
	_xResolution = +0.04; //mm/Dots .. doing so 0 is in the center point of the sensor
	_xOffset = -8.0 + _xResolution / 2.;//mm

	_hTiltNom = new TH1F("hTiltNom", "Tilt Nominator", 400, -8, 8); //This is for raw-data from 8-8
	_hTiltDenom = (TH1F*)_hTiltNom->Clone("hTiltDenom");
	_hTiltDenom->SetTitle("Tilt Denominator");

	_hRangeNom = new TH2F("hRangeNom", "Range Nominator", 200, 0, 0.20, 450, -2, 16);
	_hRangeNom->GetXaxis()->SetTitle("time [s]");
	_hRangeNom->GetYaxis()->SetTitle("Pos on Sonotrode [mm]");
	
	_hRangeDenom = (TH2F*)_hRangeNom->Clone("hRangeDenom");
	_hRangeDenom->SetTitle("Amplitude Denominator");

	_hRangeNomSmooth = (TH2F*)_hRangeNom->Clone("hRangeNomSmooth");
	_hRangeNomSmooth->SetTitle("Amplitude Denominator Smooth");

	_hRangeDenomSmooth = (TH2F*)_hRangeNom->Clone("hRangeDenomSmooth");
	_hRangeDenomSmooth->SetTitle("Amplitude Denominator Smooth");

	_hRangeCrossNom = new TH1F("hRangeCrossNom", "CrossRange Nominator", 450, -2, 16);
	_hRangeCrossNom->GetXaxis()->SetTitle("Pos on Sonotrode [mm]");
	_hRangeCrossNom->GetYaxis()->SetTitle("Mean Range P2P [mm]");

	_hRangeCrossDenom = (TH1F*)_hRangeCrossNom->Clone("hRangeCrossDenom");
	_hRangeCrossDenom->SetTitle("Cross Amplitude Denom");

	_hRangeCrossNomSmooth = (TH1F*)_hRangeCrossNom->Clone("hRangeCrossDenomSmooth");
	_hRangeCrossNomSmooth->SetTitle("Cross Amplitude Denom");

	_hRangeCrossDenomSmooth = (TH1F*)_hRangeCrossNom->Clone("hRangeCrossDenomSmooth");
	_hRangeCrossDenomSmooth->SetTitle("Cross Amplitude Denom");

	_grAmplPit = new TGraph();
	//_grAmplPit->GetYaxis()->SetTitle("Amplitude P2P at Max[#mum]");
	//_grAmplPit->GetXaxis()->SetTitle("time [s]");
	_grAmplPit->SetTitle("Max Amplitude at Pit");
	
	_grNullPos = new TGraph();
	_grNullPos->SetTitle("Position of NP");
	//_grNullPos->GetXaxis()->SetTitle("time [s]");
	//_grNullPos->GetYaxis()->SetTitle("Position on Sonotrode [mm]");
	
	_grAmplPos = new TGraph();
	_grAmplPos->SetTitle("Pos of Max");
	//_grAmplPos->GetXaxis()->SetTitle("time [s]");
	//_grAmplPos->GetYaxis()->SetTitle("Position on Sonotrode [mm]");
	

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
	Float_t         rawZ[400];
	Float_t         rawTime;
	Float_t         rawSensor;


	
	TBranch* b_raw_x;   //!
	TBranch* b_raw_z;   //!
	TBranch* b_raw_time;   //!
	TBranch* b_raw_sensor;   //!

	// Set branch addresses and branch pointers

	
	_tree->SetBranchAddress("raw.x", rawX, &b_raw_x);
	_tree->SetBranchAddress("raw.z", rawZ, &b_raw_z);
	_tree->SetBranchAddress("raw.time", &rawTime, &b_raw_time);
	_tree->SetBranchAddress("raw.sensor", &rawSensor, &b_raw_sensor);

	// New TTree to store further analysed variables
	TTree* anatree = new TTree("anatree", "anatree");
	anatree->Branch("raw.x", rawX);
	anatree->Branch("raw.z", rawZ);
	anatree->Branch("raw.time", &rawTime);
	anatree->Branch("raw.sensor", &rawSensor);


	static const int nDots = 400;

	float z[nDots];
	float x[nDots];
	float deltaZ[nDots];
	float rangeZ[nDots];
	float smoothZ[nDots];
	float smoothDeltaZ[nDots];
	float smoothRangeZ[nDots];
	anatree->Branch("ana.z", &z,"ana.z[400]/F");
	anatree->Branch("ana.x", &x,"ana.x[400]/F");
	anatree->Branch("ana.deltaZ", &deltaZ,"ana.deltaZ[400]/F");
	anatree->Branch("ana.rangeZ", &rangeZ,"ana.rangeZ[400]/F");
	anatree->Branch("ana.smoothZ", &smoothZ, "ana.smoothZ[400]/F");
	anatree->Branch("ana.smoothDeltaZ", &smoothDeltaZ, "ana.smoothDeltaZ[400]/F");
	anatree->Branch("ana.smoothRangeZ", &smoothRangeZ, "ana.smoothRangeZ[400]/F");


	static const int nElements = 20;
	float lastElements[nDots][nElements];
	float lastElementsSmooth[nDots][nElements];




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
			if (rawZ[col] < -50) {
				x[col] = -999;
				z[col] = -999;
				rangeZ[col] = -999;
			}
			else {

				//Changing the orientation, by shift and rotation....in the end we will look perendicular to the sonotrode
				//first copy
				x[col] = rawX[col];
				z[col] = rawZ[col];


				//Shift, so that z(0) = 0!
				z[col] -= _offset;

				//Rotate by the negative of the estimated slope
				rotateBySlope(x[col], z[col], -_slope);

				//now more the pitpos to 0;
				//Recalculate the pitpos into the coordinate system, too
				float pitPosZ = _slope * _pitpos;
				float pitPosX = _pitpos; //Keep the porignal one
				rotateBySlope(pitPosX, pitPosZ, -_slope);
				//std::cout << "Recalculated pitpos is " << pitPosX << " " << pitPosZ << std::endl;

				x[col] -= pitPosX;

				//Move in such a way, that the sonotrode's pit is at 0

				//----------------------------------

				lastElements[col][jentry % nElements] = z[col];
				if (jentry > nElements) { //Calculate the range
					float thisMin = 9999;
					float thisMax = -9999;
					for (int nLast = 0; nLast < nElements; nLast++) {
						if (lastElements[col][nLast] < thisMin) {
							thisMin = lastElements[col][nLast];
						}
						if (lastElements[col][nLast] > thisMax) {
							thisMax = lastElements[col][nLast];
						}
					}
					if (thisMax < 50 && thisMin > -50 && thisMax > -50 && thisMin < 50) {
						rangeZ[col] = thisMax - thisMin; //This is the range
						//To calculate the mean filling the sum and the n, x is inverted!!!
						_hRangeNom->Fill(rawTime, -x[col], rangeZ[col]);
						_hRangeDenom->Fill(rawTime, -x[col]);
						_hRangeCrossNom->Fill(-x[col], rangeZ[col]);
						_hRangeCrossDenom->Fill(-x[col]);
					}
					else rangeZ[col] = -999;
					
				}
				else {
					rangeZ[col] = -999;
				}




			}

		}

		smooth(z, smoothZ, 400, 20);
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
	_hRange = (TH2F*)_hRangeNom->Clone("hRange");
	//_hAmpl->Sumw2();
	_hRange->SetTitle("Mean Range");
	_hRange->Divide(_hRangeDenom);
	/*
	//Calculation of the mean Range Smoothed
	_hRangeSmooth = (TH2F*)_hRangeNomSmooth->Clone("hRangeSmooth");
	//_hAmpl->Sumw2();
	_hRangeSmooth->SetTitle("Mean Range Smoothed");
	_hRangeSmooth->Divide(_hRangeDenomSmooth);
	*/
	//Calculation of the mean Range
	_hRangeCross = (TH1F*)_hRangeCrossNom->Clone("hRangeCross");
	//_hAmpl->Sumw2();
	_hRangeCross->SetTitle("Mean Range Cross");
	_hRangeCross->Divide(_hRangeCrossDenom);
	/*
	//Calculation of the mean Range Smoothed
	_hRangeCrossSmooth = (TH1F*)_hRangeCrossNomSmooth->Clone("hRangeCrossSmoothed");
	//_hAmpl->Sumw2();
	_hRangeCrossSmooth->SetTitle("Mean Range Cross Smoothed");
	_hRangeCrossSmooth->Divide(_hRangeDenomSmooth);
	*/
	


	//Calculation of the NP
	//TF1* fNP = new TF1("fNP", "[0]+[1]*fabs(sin([2]*x +[3]))+ [4]*x", 4, 16);
	TF1* fNP = new TF1("fNP", "[0] + fabs([1]*x+[2])", 4, 16);
	
	fNP->SetParameter(0, 0);
	fNP->SetParameter(1, 0);
	fNP->SetParameter(2, 0.05);
	fNP->SetParLimits(0, 0, 0.1);
	fNP->SetParLimits(1, -0.1, 0);
	fNP->SetParLimits(2, 0, 0.1);
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
	
	int nbins = _hRange->GetXaxis()->GetNbins();
	for (int i = 1; i <= nbins; i++) { //this is one time slice
		char nameproj[64];
		sprintf(nameproj, "_projy%d", i);
		TH1F* myProj = (TH1F*)_hRange->ProjectionY(nameproj, i, i);
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
	float z[nDots];
	float sensor;
	float time;
	int nDataSets = 0;


	

	
	_tree->Branch("raw.x", &x,"raw.x[400]/F");
	_tree->Branch("raw.z", &z,"raw.z[400]/F");
	_tree->Branch("raw.time", &time);
	_tree->Branch("raw.sensor", &sensor);
	
	std::ifstream infile(filename.c_str());
	int linenum = -1;
	std::string sline;
	std::vector <std::string> dataVec;

	if (!infile.is_open()) {
		std::cerr << "ERROR: Unable to open input file '" << filename << "'!" << std::endl;
		exit(1);
	}
	else {
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
			dataVec = split(sline,_deliminator.at(0));
			for (int i = 0; i < dataVec.size(); i++) {
				float myZ = atof(dataVec.at(i).c_str());
				float myX = i * _xResolution + _xOffset;
				if (myZ > -50 && myZ < 50) { //Convert only data in valid range
					x[i] = myX;
					z[i] = myZ;
					
					_hTiltNom->Fill(myX, myZ);
					_hTiltDenom->Fill(myX);
				}
				else {
					x[i] = -999;
					z[i] = -999;
				}
			}
			
			
			time = linenum / _frequency;
			_tree->Fill();
			nDataSets++;
			
		}

		infile.close();
		std::cout << nDataSets << " datasets from file " << filename << " analysed" << std::endl;

		//Mean of the z-values (per x) is the sum of the z-values (per x) devided by n
		_hTilt = (TH1F*)_hTiltNom->Clone("hTilt");
		//_hTilt->Sumw2();
		_hTilt->SetTitle("Mean Tilt");
		_hTilt->Divide(_hTiltDenom);

		TF1* fTilt = new TF1("fTilt", "pol1", -8, 8);
		_hTilt->Fit("fTilt");
		_offset = fTilt->GetParameter(0);
		_slope = fTilt->GetParameter(1);
		std::cout << "Tilt was determined to " << _slope << " (+/- " << fTilt->GetParError(1) << ")" << std::endl;
		std::cout << "Distance was determined to " << _offset << " (+/- " << fTilt->GetParError(0) << ")" << std::endl;


		//determine the pit position by looking for the last entries of the histogram (small entries < 1% will be ignored)
		Double_t q[1];
		Double_t p[1] = { 1-1e-5 };
		_hTiltDenom->GetQuantiles(1, q, p);
		std::cout << "Pit position located at " << q[0] << std::endl;
		_pitpos = q[0];
	}



	return _tree;


}