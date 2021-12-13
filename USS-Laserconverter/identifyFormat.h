#pragma once

#ifndef IDENTIFYFORMAT_H
#define IDENTIFYFORMAT_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "global.h"


enum DataKind {
	kUnknown,
	kOneLJ7060_400_16,
	kTwoLJ7060_800_16
		
};
	

DataKind identifyFormat(std::string filename);

#endif
