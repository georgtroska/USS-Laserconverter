#pragma once
#ifndef ARGS_H
#define ARGS_H

#include <string>
#include <vector>
class ArgParser {
	protected:
	std::string _delim1;
	std::string _delim2;
	std::string _outFilename;
	std::vector<std::string> _inputFilenameList;
	std::string _treename;
public:
	ArgParser() {}
	virtual ~ArgParser() {}
	void parse(int argc, const char* argv[]);
	void print() {
		std::cout << "Delim1: " << _delim1 << std::endl;
		std::cout << "Delim2: " << _delim2 << std::endl;
		std::cout << "outFilename: " << _outFilename << std::endl;
		std::cout << "treename: " << _treename << std::endl;
		std::cout << "inputFilenameList: ";
		for (unsigned int i = 0; i < _inputFilenameList.size(); i++) std::cout << std::string(_inputFilenameList.at(i) + " ");
		std::cout << std::endl;
	}
	std::string getDelimitor1() { return _delim1; }
	std::string getDelimitor2() { return _delim2; }
	std::string getOutFileName() { return _outFilename; }
	std::vector<std::string> getInputFiles() { return _inputFilenameList; }
	std::string getTreeName() { return _treename;  }
};




	



#endif //ARGS_H