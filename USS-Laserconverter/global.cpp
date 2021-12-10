#include "global.h"
#include "stdafx.h"
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

void split(const std::string& s, char delim, std::vector<std::string>& elems) {
	// Split a string 's' at delimiter 'delim' into a vector 'elems' 
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
}

std::vector<std::string> split(const std::string& s, char delim) {
	// Split a string 's' at delimiter 'delim' and return a vector containing the elements
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}


///Keep name conventions for ROOT TObject names
std::string RootNameConditioning(const std::string in) {
	std::string allowedCharsFirst = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_";
	std::string allowedChars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_";
	// +-/*()#[]<>%.:'"´` ... are no good idea
	std::string out;
	for (unsigned int i = 0; i < in.length(); ++i) {
		if ((i == 0 && allowedCharsFirst.find(in.at(i)) != std::string::npos) || (i > 0 && allowedChars.find(in.at(i)) != std::string::npos)) {
			out.push_back(in.at(i));
		}
	}
	return out;
}

std::string reformatTypes(const std::string in, std::string &options) {
	std::string out = in;
	//separation of format and option
	std::size_t f1 = out.find("(");
	std::size_t f2 = out.find("(");
	if (f1 != std::string::npos && f2 != std::string::npos) {
		//Copy the string between the brackets into the options
		options = std::string(out, f1 + 1, f2 - f1 - 1);
		//Erase the options string
		out.erase(out.begin() + f1, out.begin() + f1 + 1);
	}
	else if (f1 != std::string::npos) { //only "(" no ")"
		//Erase the option string
		out.erase(out.begin() + f1, out.end());
	}
	std::transform(out.begin(), out.end(), out.begin(), ::toupper);
	return out;
}
