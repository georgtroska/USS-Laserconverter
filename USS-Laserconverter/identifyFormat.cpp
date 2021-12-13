#include "identifyFormat.h"

DataKind identifyFormat(std::string filename) {
	std::ifstream infile(filename.c_str());
	std::vector <std::string> dataVec;
	int linenum = -1;

	std::string sline;

	if (!infile.is_open()) {
		std::cerr << "ERROR: Unable to open input file '" << filename << "'!" << std::endl;
		exit(1);
	}
	else {
		//std::cout << "Starting conversion!" << std::endl;
		while (!infile.eof()) {
			linenum++;
			getline(infile, sline);
			if (linenum < 10) continue;
			dataVec = split(sline, ';');
			break;
		}
		if (dataVec.size() == 400) {
			std::cout << "Assuming Data from one LJ7060, with 16KHz, and 0.04mm pitch" << std::endl;
			infile.close();
			return kOneLJ7060_400_16;
		}
		else if (dataVec.size() == 800) {
			std::cout << "Assuming Data from two LJ7060, with 16KHz, and 0.04mm pitch" << std::endl;
			infile.close();
			return kTwoLJ7060_800_16;
		}
	}
	return kUnknown;


}