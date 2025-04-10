#ifndef SUB_ANALYSIS
#define SUB_ANALYSIS

#include <array>
#include <iostream>
#include <vector>
#define NR_OF_COLUMNS 3

/** Subanalysis class
 *  A part of analysis configuration corresponds to
 *  <pre>
 *   subparameters 3
 *      !First set
 *      120  !particle id
 *      -1.0 !min rapidity
 *       1.0 !max rapidity
 *      !Second set
 *      2130 !particle id
 *      -1.0 !min rapidity
 *       1.0 !max rapidity
 *  </pre>
 *  This class contains informations about subparameters used in the
 * Analysis class. Each set of subparameters corresponds to an instance
 * of the Subanalysis class. Each set of subparameters is stored in a
 * c++ vector of a size given by subparameters value. The subparameters
 * are used to select the EPOS output data to build the histograms.
 *  The corresponding histogram values are stored in an array, a vector composed
 * of NR_OF_COLUMNS float rows.
 */

class Subanalysis {
private:
  std::string meaning;

  /**
   * @brief set of subparameters of type float
   *
   */
  std::vector<float> subparameters;

  /**
   * @brief array of data : array of number of bins rows and NR_OF_COLUMNS
   * columns
   *
   */
  std::vector<std::array<float, NR_OF_COLUMNS>> binCounts;

  unsigned int numberOfTriggeredEvents;

public:
  /**
   * @brief Construct a new Parameter Set object
   *
   */
  Subanalysis();

  /**
   * @brief Construct a new Parameter Set object
   *
   * @param nbOfSubparameters : number of subparameters
   */
  Subanalysis(const int nbOfSubparameters);

  /**
   * @brief Construct a new Parameter Set object
   *
   * @param nbOfBins : number of bins
   * @param nbOfSubparameters : number of subparameters
   */
  Subanalysis(const int nbOfBins, const int nbOfSubparameters);

  /**
   * @brief Construct a new Parameter Set object
   *
   * @param meaning : meaning of the subanalysis
   * @param nbOfBins : number of bins
   * @param nbOfSubparameters : number of sub parameters
   */
  Subanalysis(const std::string meaning, const int nbOfBins,
              const int nbOfSubparameters);

  /**
   * @brief Destroy the Parameter Set object
   *
   */
  ~Subanalysis();

  /**
   * @brief Get subparameter value
   *
   * @param index
   * @return float value
   */
  float getSubparametersValue(unsigned int index) const;

  /**
   * @brief Set the subparameters object
   *
   * @param subparametersParam
   */
  void setSubparameters(const std::vector<float> subparametersParam);

  /**
   * @brief Get the subparameters object
   *
   * @return std::vector<float>
   */
  std::vector<float> getSubparameters() const;

  /**
   * @brief Create a vector of nbOfBins arrays of NR_OF_COLUMNS floats
   *
   */
  void createBinCounts(unsigned int nbOfBins);

  /**
   * @brief Set the BinCounts object
   *
   * @param binCounts
   */
  void setBinCounts(std::vector<std::array<float, NR_OF_COLUMNS>> &binCounts);

  /**
   * @brief Return a vector of nbOfBins array of NR_OF_COLUMNS floats
   *
   */
  std::vector<std::array<float, NR_OF_COLUMNS>> getBinCounts() const;

  /**
   * @brief Return meaning of subanalysis
   *
   */
  std::string getMeaning() const;

  /**
   * @brief Set meaning of subanalysis parameters
   *
   * @param meaning
   *
   */
  void setMeaning(std::string meaning);

  /**
   * @brief Initialize a vector of number of bins binCounts of NR_OF_COLUMNS
   * floats to bin values for the first column and 0 for the second and third
   * one
   *
   * @param binvalues
   */
  void initializeBinCounts(std::vector<float> binvalues);

  /**
   * @brief Set the Number Of Triggered Events object
   *
   * @param numberOfTriggeredEvents
   */
  void setNumberOfTriggeredEvents(unsigned int numberOfTriggeredEvents);

  /**
   * @brief Get the Number Of Triggered Events object
   *
   * @return unsigned int
   */
  unsigned int getNumberOfTriggeredEvents();

  /**
   * @brief increment the number of triggered events
   *
   */
  void incrementNumberOfTriggeredEvents();

  /**
   * @brief Update value and error in the binCounts of Subanalysis
   *
   * @param binNb bin number
   * @param value
   * @param error
   */
  void update(unsigned int binNb, float value, float error = 1.);

  /**
   * @brief Reset the binCounts values of Subanalysis to 0.
   *
   */
  void reset();

  /**
   * @brief Print binCounts values according to histo format
   *
   * @param os
   */
  void printBinCounts(std::ostream &os);

  /**
   * @brief normalize the binCounts data
   *
   * @param factor
   * @param normalizeByBinCounts
   */
  void normalize(float factor, bool normalizeByBinCounts);

  /**
   * @brief Get the Row Number object
   *
   * @return unsigned int
   */
  unsigned int getRowNumber();

  /**
   * @brief Get the BinCounts Values object
   *
   * @param rowNumber
   * @return float* : array of NR_OF_COLUMNS values corresponding to row
   * rowNumber
   */
  float *getBinCountsValues(int rowNumber);
};

/**
 * @brief overload the operator <<
 *
 * @param os output stream
 * @param subanalysis parameter set to display
 * @return std::ostream&
 */
std::ostream &operator<<(std::ostream &os, const Subanalysis &subanalysis);

#endif
