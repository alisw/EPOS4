//
// This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
// (See COPYING file for the text of the licence)
//

#include <iostream>

void checkTime();
void defineStorageSettings();
void initializeElectronProtonPart();                              
void eposEnd();
void eposStart();
void finalizeSimulation();
void generateEposEvent(int&);
void initializeEpos();
void initializeEventCounters();
void initializeHeavyQuarkPart();
int numberOfEnergyValues();
int numberOfEvents();
void readInputFile();
void rewindInputFile();
int setEnergyIndex(int&);
void showMemoryAtStart();
void showMemoryAtEnd();
void writeStatistics();


int main( int argc, char *argv[]){
  std::cout << "start main in main.cpp" << std::endl;

  showMemoryAtStart();
  checkTime();

  eposStart();
  readInputFile();
  initializeHeavyQuarkPart();
  initializeElectronProtonPart();                              

  for (int k=1; k<=numberOfEnergyValues(); k++){
    setEnergyIndex(k);
    initializeEpos();
    initializeEventCounters();
    defineStorageSettings();
    for (int n=1; n<=numberOfEvents(); n++){
      generateEposEvent(n);
    }   
    writeStatistics();
  }
  finalizeSimulation();
  rewindInputFile();
  readInputFile();
  eposEnd();

  showMemoryAtEnd();
}

