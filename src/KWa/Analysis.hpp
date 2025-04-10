#ifndef ANALYSIS
#define ANALYSIS

#include "Subanalysis.hpp"
#include <cmath>
#include <iostream>

/** Analysis class.
 *  The Analysis class represents an analysis of simulation data, in general
 * composed of several subanalyses, characterized by a set of "common
 * parameters", and several sets of "additional parameters defining the
 * different subanalyses.\n There are at least four common parameters, and the
 * first four have a well defined meaning (minvalue, maxvalue, nrbins,
 * binning, normalization).\n The parameters which define the analyses
 * are provided via the standard EPOS input (optns file), using the syntax:
 *  <pre>
 *  beginanalysis
 *    name \<name of analysis\>
 *    observable \<name of observable\>
 *    minvalue \<max value\>
 *    maxvalue \<min value\>
 *    nrbins \<number of bins\>
 *    binning \<lin|log\>
 *    normalization 1
 *    moreparameters \<number of more common parameters\>
 *      \<liste of more common parameters\>
 *    subanalyses \<number of subanalyses\>
 *    subparameters \<number of additional parameters\>
 *      \<liste of subanalyses\>
 *  endanalysis
 *  </pre>
 *
 * for example:
 *  <pre>
 *  beginanalysis
 *   name PtSpectra
 *   minvalue 0.
 *   maxvalue 10.
 *   nrbins   10
 *   binning lin
 *   normalization 1
 *   moreparameters 1
 *     3     !another common parameter
 *   subanalyses 2 !which allows to make similar analyses in parallel
 *   subparameters 3
 *      !First set
 *      120  !particle id
 *      -1.0 !min rapidity
 *       1.0 !max rapidity
 *      !Second set
 *      2130 !particle id
 *      -1.0 !min rapidity
 *       1.0 !max rapidity
 *  endanalysis
 *  </pre>
 *
 *  Each occurence  of such a "named analysis structure" will result in the
 * creation of an Analysis instance with the corresponding name, containing the
 * parameter information.\n Each analysis amounts to updating tables with
 * nr_of_bins rows and nr_of_columns columns (one table per subanalysis) at the
 * end of each event simulation, by calling Analysis::analyse.\n The latter is
 * defined in the child class of name <name of analysis>.\n
 */

class Analysis {
private:
  /**
   * @brief analysis name
   *        This name corresponds to the child class type.
   *
   */
  const std::string name;

  /**
   * @brief observable
   *        This name corresponds to the observable name.
   *
   */
  const std::string observable;

  /**
   * @brief minvalue
   *        Minimum X value.
   *
   */
  const float minvalue;

  /**
   * @brief maxvalue
   *        Maximum X value.
   *
   */
  const float maxvalue;

  /**
   * @brief numberOfBins
   *        Number of bins in histogram.
   *
   */
  const int numberOfBins;

  /**
   * @brief binning
   *        Binning type : linear or logarithmic
   *
   */
  const std::string binning;

  /**
   * @brief normalization
   *        This name corresponds to the child class type.
   *
   */
  const int normalization;

  /**
   * @brief vector of vector of imported data
   *        Data are imported from a histogram defined in a histofile
   *        These imported data are common to all subanalyses.
   */
  std::vector<std::vector<float>> importarray;

  /**
   * @brief vector of common parameters
   *        These parameters are common to all subanalyses.
   */
  std::vector<float> moreparameters;

  /**
   * @brief vector of subanalyses
   *        These sets of parameters allow to add more filters and triggers.
   * conditions.
   *
   */
  std::vector<Subanalysis> subanalyses;

  /**
   * @brief histoweight
   *
   */
  double histoweight;

public:
  /**
   * @brief Construct a new Analysis object
   *
   */
  Analysis();

  /**
   * @brief Construct a new Analysis object
   *
   * @param name : child class type
   * @param minvalue : min value
   * @param maxvalue : max value
   * @param nrBins : number of bins
   * @param binning : binning type - linear (lin) or logarithmic (log)
   * @param normalization : normalization
   * @param nbMoreparameters : number of additionnal common parameters
   * @param nbSubanalyses : number of subanalyses equals to the number of
   * parameter sets
   */
  Analysis(std::string name, float minvalue, float maxvalue, int nrBins,
           std::string binning, int normalization, int nbMoreparameters,
           int nbSubanalyses);

  /**
   * @brief Construct a new Analysis object
   *
   * @param name : child class type
   * @param observable : observale name
   * @param minvalue : min value
   * @param maxvalue : max value
   * @param nrBins : number of bins
   * @param binning : binning type - linear (lin) or logarithmic (log)
   * @param normalization : normalization
   * @param nbMoreparameters : number of additional common parameters at least 4
   * @param nbSubanalyses : number of subanalyses equals to the number of
   * parameter sets
   */
  Analysis(std::string name, std::string observable, float minvalue,
           float maxvalue, int nrBins, std::string binning, int normalization,
           int nbMoreparameters, int nbSubanalyses);

  /**
   * @brief Construct a new Analysis object
   *
   * @param name : child class type
   * @param observable : observale name
   * @param minvalue : min value
   * @param maxvalue : max value
   * @param nbBins : number of bins
   * @param binning : binning type - linear (lin) or logarithmic (log)
   * @param normalization : normalization
   * @param histofilename : histo file name
   * @param histoname : histo name
   * @param nbMoreparameters : number of additional common parameters at least 4
   * @param nbSubanalyses : number of subanalyses equals to the number of
   * parameter sets
   */
  Analysis(std::string name, std::string observable, float minvalue,
	   float maxvalue, int nbBins, std::string binning,
	   int normalization,
	   std::string histofilename, std::string histoname, 
	   int nbMoreparameters, int nbSubanalyses);
  
  /**
   * @brief Destroy the Analysis object
   *
   */
  virtual ~Analysis() = default;

  /**
   * @brief Get common parameter value
   *
   * @param index
   * @return float value
   */
  float getMoreparametersValue(unsigned int index) const;

  /**
   * @brief Get the name object
   *
   * @return std::string
   */
  std::string getName() const;

  /**
   * @brief Get the observable object
   *
   * @return std::string
   */
  std::string getObservable() const;

  /**
   * @brief Get the Xmin object
   *
   * @return float
   */
  float getMinvalue() const;

  /**
   * @brief Get the Xmax object
   *
   * @return float
   */
  float getMaxvalue() const;

  /**
   * @brief Get the Number Of Bins object - third common parameter
   *
   * @return unsigned int
   */
  unsigned int getNumberOfBins() const;

  /**
   * @brief Get the Normalization object - fouth common parameter
   * <pre>
   * // c     the normalization factor is a number jk which normalizes to:
   * // c
   * // c  k  0:  * 1
   * // c     1:  / number of events
   * // c     2:  / number of triggered events
   * // c     4:  / bin-counts
   * // c
   * // c  j  0:  * 1
   * // c     1:  / bin-width
   * // c
   * </pre>
   *
   * @return unsigned int
   */
  unsigned int getNormalization() const;

  /**
   * @brief Get the Histoweight object
   *
   * @return double
   */
  double getHistoweight() const;

  /**
   * @brief Set the Histoweight object
   *
   * @param histoweight
   */
  void setHistoweight(double histoweight);

  /**
   * @brief Get the Number Of Parameter Sets object
   *
   * @return unsigned int
   */
  unsigned int getNumberOfSubanalyses() const;

  /**
   * @brief Get the importarray object
   *
   * @return std::vector<std::vector<float>>
   */
  std::vector<std::vector<float>> getImportarray() const;

  /**
   * @brief Set the moreparameters object
   *
   * @param moreparametersParam
   */
  void setMoreparameters(std::vector<float> moreparametersParam);

  /**
   * @brief Get the moreparameters object
   *
   * @return std::vector<float>
   */
  std::vector<float> getMoreparameters() const;

  /**
   * @brief Set the Subanalyses object
   *
   * @param subanalysesParam
   */
  void setSubanalyses(std::vector<Subanalysis> subanalysesParam);

  /**
   * @brief Get the Parameter Sets object
   *
   * @return std::vector<Subanalysis>
   */
  std::vector<Subanalysis> &getSubanalyses();

  /**
   * @brief
   *
   * @param subanalysisIndex
   * @param meaning
   * @param subparameters
   */
  void addSubparameters(unsigned int subanalysisIndex, std::string meaning,
                        std::vector<float> subparameters);

  /**
   * @brief Initialize importarray
   * Extract data from histoname defined in histofilename
   */  
  void initializeImportarray(std::string histofilename, std::string histoname);
     
  /**
   * @brief Initialize binCounts for each subanalysis setting bin values in
   * first column
   *
   */
  void initializeSubanalysesBinCounts();

  /**
   * @brief Normalize data in binCounts
   * <pre>
   * // c       from src/KW/xan.f
   * // c.......here normalization.......................................
   * // c           see also   "..........fill histogram"
   * // c.................................................................
   * // c     the norm ( inorm(n) ) is a number hijk which normalizes to:
   * // c
   * // c  k  0:  * 1
   * // c     1:  / number of events
   * // c     2:  / number of triggered events
   * // c     4:  / bin-counts
   * // c
   * // c  j  0:  * 1
   * // c     1:  / bin-width
   * // c
   * // c  h and i are not yet implemented
   * // c.................................................................
   * </pre>
   *
   * @param numberOfEvents
   */
  void normalize(unsigned int numberOfEvents);

  /**
   * @brief Print binCounts values according to histo format
   *
   * @param subanalysisIndex : index of the parameter set
   * @param os                : output stream (file, str::out...)
   */
  void printBinCounts(unsigned int subanalysisIndex, std::ostream &os);

  /**
   * @brief Get the Bin Nb object
   *        Index of the bin corrsponding to the float value
   *
   * @param value
   * @return unsigned int
   */
  unsigned int getBinNb(float value);

  /**
   * @brief Build the histograms according to output particules properties
   *        This method is implemented only in child class.
   */
  virtual void analyze() = 0;
};

/**
 * @brief
 *
 * @param os
 * @param analysis
 * @return std::ostream&
 */
std::ostream &operator<<(std::ostream &os, Analysis &analysis);

/**
 * @brief Get the Analysis object
 *
 * @param analysisName
 * @return Analysis*
 */
Analysis *getAnalysis(std::string analysisName);

/**
 * @brief Get the Analyses Vector object
 *
 * @return std::vector<Analysis*>
 */
std::vector<Analysis *> getAnalysesVector();

#endif
