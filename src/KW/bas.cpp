//
// This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
// (See COPYING file for the text of the licence)
//

#include <iostream>


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
  int ifevent_();
  void setebin_(int *);
  int ifebin_();
}


/**
 * Show memory usage at the start of the process
 *
 */
void showMemoryAtStart(){
  const char* start_str = (const char *)"start program;";
  int indice = 0;
  memo_(&indice, &start_str[0], (int)std::char_traits<char>::length(start_str));
} 


/**
 * Show memory usage at the end of the process
 *
 */
void showMemoryAtEnd(){
  const char* stop_str = (const char *)("stop program;");
  int indice = 0;
  memo_(&indice, &stop_str[0], (int)std::char_traits<char>::length(stop_str));
} 


/**
 * Execution time
 *
 */
void checkTime(){
  const char* start_str = (const char *)"start program;";
  checktime_(&start_str[0], (int)std::char_traits<char>::length(start_str));
}


/**
 * Start EPOS
 *
 */
void eposStart(){
  eposstart_();
}


/**
 * Read the input file
 *
 */
void readInputFile(){
  aread_();
}


/**
 * Initialize EPOS
 *
 */
void initializeEpos(){
  ainit_();
}

    
/**
 * Initialize the event counters
 *
 */
void initializeEventCounters(){
  eventsini_();
}


/**
 * Define the storage settings
 *
 */
void defineStorageSettings(){
  bstora_();
}


/**
 * Write statistics
 *
 * @param n event number
 *
 */
void generateEposEvent(int& n){
  eposevent_(&n);
}


/**
 * Finalize simulation
 *
 */
void finalizeSimulation(){
  bfinal_();     
}                         


/**
 * Rewind the input file
 *
 */
void rewindInputFile(){
  swopen_();  
}    
     
                  
/**
 * End EPOS
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

