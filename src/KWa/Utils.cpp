#include <stdio.h>  // defines FILENAME_MAX
#include <unistd.h> // for getcwd()
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void displayValues(std::vector<std::vector<float>> values) {
  for (auto row : values) {
    for (auto element : row) {
      std::cout << element << '\t';
    }
    std::cout << '\n';
  }
}

std::string extractValue(std::string stringOfValues, size_t &beginIndex) {
  const char delimiter = ' ';
  while (stringOfValues[beginIndex] == delimiter) {
    beginIndex++;
  }

  // extract token between two delimiters
  size_t endIndex = stringOfValues.find(delimiter, beginIndex);
  std::string token = stringOfValues.substr(beginIndex, endIndex - beginIndex);

  // set beginIndex to the last delimiter, for the next extraction
  beginIndex = endIndex;

  return token;
}

std::vector<std::vector<float>> readValuesFromString(std::string arrayValues) {
  // extract number of columns
  size_t beginIndex = 0;
  std::string token;

  const unsigned int numberOfColumns =
      std::stoi(extractValue(arrayValues, beginIndex));

  std::vector<std::vector<float>> values;
  int i = 0;

  // extract array values
  while (beginIndex < arrayValues.length()) {
    std::vector<float> row;
    for (unsigned int i = 0; i < numberOfColumns; i++) {
      row.push_back(std::stof(extractValue(arrayValues, beginIndex)));
    }
    values.push_back(row);
  }

  return values;
}

// read values from histo file filename
// search for data corresponding to histogram named histoname
std::string readValuesFromFile(std::string filename, std::string histoname) {
  // Open the input file
  std::string absolutePathFilename;
  const char* variableName = "EPO";
  if (const char* path = std::getenv(variableName)) {
    absolutePathFilename = path + filename;
    std::cout << "Import array from file : " << absolutePathFilename << '\n';
  } else {
    std::cerr << "The environment variable " << variableName << " is not defined..." << std::endl;
    exit(1);
  }
  
  std::ifstream inputFile(absolutePathFilename);

  // Check if the file is successfully opened
  if (!inputFile.is_open()) {
    std::cerr << "Error opening the file for importarray command: " << absolutePathFilename << std::endl;
    exit(1);
  }

  std::string line; // Declare a string variable to store each line of the file
  std::string beginArray = "array";
  std::string endArray = "endarray";
  std::string arrayString, arrayValues;
  bool histonameIsFound = false;
  bool arrayIsFound = false;

  // Read each line of the file
  // search histoname
  // extract data between "openhisto name 'histoname'" and "endarray"
  while ((getline(inputFile, line)) && (!arrayIsFound)) {
    if (line.length() > 0) { // skip empty lines
      if (!histonameIsFound) {
        size_t histonameIndex = line.find(histoname);
	if (histonameIndex != std::string::npos) {
	    size_t nextCharacterIndex = line.find(histoname) + histoname.length();
	    
	    // histoname was found in line
	    // the found string matchs exactly histoname
	    bool noCharacterBefore =
	      (histonameIndex == 0) || (line[histonameIndex - 1] == ' ');
	    bool noCharacterAfter = (line[nextCharacterIndex] == ' ') ||
	      (nextCharacterIndex == line.length());
	    if (noCharacterBefore && noCharacterAfter) {
	      histonameIsFound = true;
	    }
	  }
        // if histoname is not found, the line is skipped
      } else {
        size_t endArrayIndex = line.find(endArray);
        if (endArrayIndex != std::string::npos) {
          arrayIsFound = true;
          // stop reading the file, line by line
        } else {
          // if histoname was found but endArray is not yet found
          // the line is added to arrayString
          arrayString += line;
        }
      }
    }
  }
  // Close the file
  inputFile.close();

  // extract the values from beginArray
  // values are stored in a string beginning by the number of columns
  size_t beginArrayIndex = arrayString.find(beginArray) + beginArray.length();
  if ((beginArrayIndex != std::string::npos)) {
    arrayValues = arrayString.substr(beginArrayIndex);
  }

  return arrayValues;
}

