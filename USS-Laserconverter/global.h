#pragma once
#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include <vector>
#include <sstream>
void split(const std::string& s, char delim, std::vector<std::string>& elems);
std::vector<std::string> split(const std::string& s, char delim);
///Keep name conventions for ROOT TObject names
std::string RootNameConditioning(const std::string in);
std::string reformatTypes(const std::string in, std::string& options);



#endif