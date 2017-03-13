#include "ConfigParams.h"
#include <iostream>

namespace usac
{

bool ConfigParams::initUSACParamsFromConfigFile()
{
	// read in parameters 
}

bool ConfigParams::initParamsFromConfigFile(std::string& configFilePath) {
	std::cout << "Implement this in the derived class" << std::endl;
	return false;
}

}
