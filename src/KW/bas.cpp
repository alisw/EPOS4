//
// This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
// (See COPYING file for the text of the licence)
//

#include <iostream>
#include <sstream>
#include "Analysis.hpp"
#include "Ico.hpp"
#include "Outlist.hpp"
using std::string;

int numberOfEvents();
void normalizeAnalysis(int numberOfEvents);
void listParticlesCheck(string,int,int,int,string);

extern "C" void aaalist_(char *text, int *n, int *m, char *file);

extern "C" {
  void memo_(int *, const char*, int);
  void checktime_(const char*, int);
  void eposstart_();
  void aread_();
  void ainit_();                   
  void eventsini_();
  void eposevent_(int *);
  void bstora_();
  void bfinal_();                              
  void swopen_();                              
  void eposend_();
  void astati_();
  void closecccptl_();
  int ifevent_();
  void setebin_(int *);
  int ifebin_();
}

extern Outlist *outlist;

/**
 * Get OutputListScreen flag:
 *
 *     1 - To activate printout of particle list to the screen
 *     0 - No action
 *
 */
void getOutputListScreen(int *flag){
  *flag=outlist->getToscreen();
}

/**
 * Get OutputListFile flag:
 *
 *     1 - To activate printout of particle list into some file
 *     0 - No action
 *
 * Get name of the file
 *
 */
void getOutputListFile(int *flag, string& file){
  int le=outlist->getTofile();
  if(le>0){
    *flag=1;
    //Array<char> * name = outlist->getThename();
    file="";
    for(int i = 1; i < le+1; i++){
      char ch = outlist->getThenameElem(i);
      file+=ch;
    }
  }else{
    *flag=0;
  }
}

/**
 * Destroy cccptl object (EPOS particle list).
 *
 */
void closeCccptl(){
  closecccptl_();
}

/**
Print EPOS particle list into check file.
\param text Some text explaining particle list
\param n First index to be printed
\param m Last index to be printed
\param file Name of the check file
*/
void aaalist_(char *text, int *n, int *m, char *file){
  string s = text;
  int i=s.find('&');
  string su = s.substr(0,i); 
  string fi=file;
  //std::cout<<"AAALIST"<< su <<" "<<i<<" "<< *n<<" "<<*m<<" "<< fi <<std::endl;
  listParticlesCheck(su,i,*n,*m,fi);
}

/**
 * Show memory usage at the start of the process.
 *
 */
void showMemoryAtStart(){
  const char* start_str = (const char *)"start program;";
  int indice = 0;
  memo_(&indice, &start_str[0], (int)std::char_traits<char>::length(start_str));
} 


/**
 * Show memory usage at the end of the process.
 *
 */
void showMemoryAtEnd(){
  const char* stop_str = (const char *)("stop program;");
  int indice = 0;
  memo_(&indice, &stop_str[0], (int)std::char_traits<char>::length(stop_str));
} 


/**
 * Initialize execution time.
 *
 */
void checkTime(){
  const char* start_str = (const char *)"start program;";
  checktime_(&start_str[0], (int)std::char_traits<char>::length(start_str));
}


/**
 * Initialize EPOS - first round.
 *
 */
void eposStart(){
  eposstart_();
}


/**
 * Read the input file (optns file), containing the configuration
 * concerning the reaction and the output format. 
 *
 */
void readInputFile(){
  aread_();
}


/**
 * Initialize EPOS - second round.
 *
 */
void initializeEpos(){
  ainit_();
}

    
/**
 * Initialize the event counters.
 *
 */
void initializeEventCounters(){
  eventsini_();
}


/**
 * Define the storage settings.
 *
 */
void defineStorageSettings(){
  bstora_();
}


/**
 * Genarate EPOS event.
 *
 * @param n event number
 *
 */
void generateEposEvent(int& n){
  eposevent_(&n);
  // the map analysis_function is defined in src/KWa/AnalysisFunctions.hpp
  // the following loop iterates the map and executes all analysis functions.
  std::vector<Analysis*> analysesVector = getAnalysesVector();
  for (Analysis* analysis : analysesVector) {
    analysis->analyze();
  }
  closeCccptl();
}


/**
 * Finalize simulation.
 *
 */
void finalizeSimulation(){
  bfinal_();
  // normalize dataset before writing to the file
  //  int analysisIdx = 0;
  //  normalizeAnalysis("PtSpectra", numberOfEvents());  
  // the map analysis_function is defined in src/KWa/AnalysisFunctions.hpp
  // the following loop iterates the map and executes all analysis functions.
  normalizeAnalysis(numberOfEvents());  
}                         


/**
 * Rewind the input file.
 *
 */
void rewindInputFile(){
  swopen_();  
}    
     
                  
/**
 * Delete all C++ objects.
 *
 */
void eposEnd(){
  eposend_();                             
}

    
/**
 * Write statistics
 *
 */
void writeStatistics(){
  astati_();
}


/**
 * Get the total number of events
 *
 * @return Total number of Events
 *
 */
int numberOfEvents(){
  return ifevent_();
}


/**
 * Get the total number of energy values
 *
 * @return Total number of energy values
 *
 */
int numberOfEnergyValues(){
  return ifebin_();
}


/**
 * Set the energy index
 *
 * @param n Energy index
 */
void setEnergyIndex(int& n){
  setebin_(&n);
}


