//#include <fstream>
//#include <iostream>
#include <string>
#include <vector>

void displayValues(std::vector<std::vector<float>> values);

std::string extractValue(std::string stringOfValues, size_t &beginIndex);

std::vector<std::vector<float>> readValuesFromString(std::string arrayValues);
  
// read values from histo file filename
// search for data corresponding to histogram named histoname
std::string readValuesFromFile(std::string filename, std::string histoname);
