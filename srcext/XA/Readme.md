(A nice and simple way to visualize this file is "grip -b Readme.md")

### New open source analysis tools

Analysis tools exist in the EPOS framework since a long time, but essentially for "internal use" by the developers, not well documented, and almost impossible to extend by external collaborators due to the complexity of the Fortran source code. It was therefore decided to develop a new tool, call "Named analysis", based on a very transparent structure, programmed in C++, with some basis classes defining the general structure, and child classes for particular analyses. Some of these child classes are already prepared, as examples, with the possibility to easily add new analyses by simply creating a new file containing a new child class analysis method. 

### Philosophy of the EPOS analysis tools

The main idea is to have a program in the form of a C++ child class, which represents not just a particular analysis, but a quite general analysis method, which allows to accomodate many individual analyses. The specification of an individual and concrete analysis should be done with the help of some text with a well-defined syntax in the EPOS input file (optns file). So creating a new particular analysis does not require any new code, but some more text lines in the the optns file - provided of course that a sufficiently general analysis method exists already.

### The analysis configurations in the optns file

The above-mentioned text lines in the EPOS input file (the optns file) are referred to as "analysis configuration", characterized always by the structure

       beginanalysis
         name <name_of_analysis>
         <commands>
       endanalysis

with <name_of_analysis> being the name of the analysis and  with \<commands\> characterizing in a simple and transparent way the analysis of simulation data, usually a distribution of the number of particles (or events) as a function of some variable. 
This analysis structure is followed by some write statements (write, writearry). 

The input file is read twice: 
1) The first iteration (in the very beginning) stores all the information provided by \<commands\>, to be used later (at the end of each event simulation) to do the analysis of the current event. The write statements are ignored.
2) In the second iteration, the analysis part (beginanalysis ... endanalysis) is ignored, here the write staments become important, to write some text (write) or to write the result of the analysis in the form of a table (writearray).

Any number of analyses and subsequent write statement may occur, based on one or several analysis methods. 

A very important element is the name of the analysis: it must correspond to the name of the above-mentioned child class, which defines an analysis method, the latter defined in srcext/XA/\<name_of_analysis\>.cpp

So eventually we have a limited number of (well tested) analysis methods defined in C++ files in src/XA/ (one class per file), but a huge number of possible analyses via the analysis definition in the optns file.

### An example of an analysis configuration

To explain the precise structure of an analysis configuration, we present a typical example:

       beginanalysis
         name ParticleSpectraMultiplicityTrigger
         observable ptr
         minvalue 0.    
         maxvalue 10.  
         nrbins   10   
         binning lin              !lin or log
         normalization 11         !divide by number of events and bin width   
         moreparameters 3
            -1.0  1.0  1
         subanalyses 4          !allows to make several similar analyses in parallel
            meaning TrgRapTrgMult  !allows different choices for parameter interpretation
            subparameters 7 !id,yMin,yMax,multMin,multMax,etaMinMult,etaMaxMult
               120   -1.0 1.0   0 50 -2 2   !subanalysis 1
               2130  -1.0 1.0   0 50 -2 2   !subanalysis 2 
               120   -1.0 1.0   0 10 -2 2   !subanalysis 3
               2130  -1.0 1.0   0 10 -2 2   !subanalysis 4 
       endanalysis
       write "ParticleSpectraMultiplicityTrigger distribution, subanalyses 1,2,3,4:" !just text
       writearray -subanalysis 1 3      !result, subanalysis 1, three columns (ptr, dn/dptr, error)
       writearray -subanalysis 2 3      !result, subanalysis 2
       writearray -subanalysis 3 3      !result, subanalysis 3
       writearray -subanalysis 4 3      !result, subanalysis 4

Most important is the name, here ParticleSpectraMultiplicityTrigger, which makes the link to the analysis child class of the same name, which contains all the details of the analysis method. Here we consider a transverse momentum distribution (ptr) in the range 0 to 10, with 10 bins, using a linear binning (lin). 

The keyword moreparameters allows to define more parameters with no predefined meaning, they are defined in the class ParticleSpectraMultiplicityTrigger. 

We have four subanalyses corresponding to four cases: In the first two (subanalyses 1,2), we consider pi+ (120) and lambdas (2130) with rapidity between -1 and 1, by requiring the charged multiplicity to be between 0 and 50, with the multiplicity defined for pseudorapidities between -2 and 2. The subanalyses 3 and 4 refer to the same situation, but with multiplicities between 0 and 10. As results we will se in the output (histo) file the text "ParticleSpectraMultiplicityTrigger distribution, subanalyses 1,2,3,4:" followed by four tables representing the four ptr distributions. The parameters have no predefined meaning, they are defined in the class according to the desired analysis method. 

The keyword "meaning" allows to define a string (any string is allowed) which characterizes the meaning of the parameter sets of the subanalyses. This string will be known inside the analysis child class, which allows to define several options concerning the interpretation of the parameters. So the string may be chosen, but should correspond to what is used in the analysis child class. 

## Analysis configuration and analysis child class

The procedure to create a new analysis method amounts to first write an **analysis configuration** (beginanalysis ... endanalysis), like the one above, and then create the corresponding **child class**. Several parameters are well defined, via keywords (observable, minvalue, maxvalue, etc), but others are not, in order to allow a large flexibility of the method. The keyword **moreparameters** allows to add more parameters, with no predefined meaning. The same is true for the keyword **subparameters**. It is up to the author of the child class to define the meaning of these parameters, and interpret them correspondingly in the code. 

## More or less general analysis child classes

The analysis structure allows creating child classes which represent a large number of possible analyses, like ParticleSpectraMultiplicityTrigger, which accommodates different types of particle distributions with different trigger conditions, based on one single class. This should be the 'standard' approach, in order to avoid a huge amount of class definitions (and the corresponding files). Nevertheless, it is possible to create simple analyses, for a very specific application, like

       beginanalysis
         name PtSpectra
         minvalue 0.    
         maxvalue 10.  
         nrbins   10   
         normalization 11  
       endanalysis
       write "PtSpectra" 
       writearray -subanalysis 1 3 

where we consider a transverse momentum distribution in the range 0 to 10, with 10 bins. But the same result may be obtained by using ParticleSpectraMultiplicityTrigger and choosing 'observable ptr'.

## The Analysis class

The above-mentioned analyses are based on the **class Analysis**, defined in src/KWa/Analysis.hpp, where the private members contain parameters which characterize the analysis (like observable,  minvalue,  maxvalue, etc), as well as a vector of subanalyses of type Subanalysis, with the **class Subanalysis** defined in Subanalysis.hpp. As already mentioned, defining several subanalyses allows  the possibility to have a single program in a single file being able to cover several analyses, being similar but using different parameter sets. The public members of the class Analysis are essentially get and set methods to access private attributes.

The **class Subanalysis** contains as private members the subanalysis parameters, and binCounts, which is of the type  
 
        std::vector<type>   

with  "type" being 
 
        std::array<float, NR_OF_COLUMNS>. 
        
So  binCounts represents a table with NR_OF_COLUMNS columns and a number of rows corresponding eventually to the number of bins. The table will be updated after each event, and contains at the end the final distribution. Another  private member is the numberOfTriggeredEvents. The public members of the class Subanalysis are essentially get and set methods to access private attributes.

## Child classes

There is another public member of the class Analysis: the method "analyze()", which is only defined in **child classes**. Whereas the classes Analysis and Subanalysis represent the general structure, we create for each analysis method, like ParticleSpectraMultiplicityTrigger, a child class, derived from Analysis, and the particular analysis for this case is encoded in srcext/XA/ParticleSpectraMultiplicityTrigger.cpp, where one defines the method ParticleSpectraMultiplicityTrigger::analyze().  

Although the analysis code is open sorce, and new features may be added to Analysis and Subanalysis (and AnalysisInterface.cpp), the main objective is the creation of new child classes, as explained in the following. In that case, it is not necessary to understand all the details of the basic code in src/KWa.

## Create new analysis via modifying MyAnalysis

The simplest way to create a new analysis is to modify MyAnalysis.cpp and correspondingly the beginanalysis ... endanalysis input in the optns file.

## Create new analysis child class

To create a really new analysis, say NewAnalysis, one has to create a file NewAnalysis.cpp, which contains the definition of the method NewAnalysis::analyze(), containing the complete analysis. In addition one needs to create include/NewAnalysis.hpp, which defines the class NewAnalysis as child class derived from Analysis. And one has to add few lines in some interface file, as explained in the following. Probably the best way is to take one of the examples in this directory, create new files via copy/paste, and then  apply modifications. So here is the procedure:


1) Copy one of the existing analyses 
     (like ParticleSpectraMultiplicityTrigger.cpp or a simple one like PtSpectra.cpp) 
      to NewAnalysis.cpp

2) Copy the corresponding include file in the include/ directory
     to include/NewAnalysis.hpp.

3) Change the first line of NewAnalysis.cpp  to: #include "NewAnalysis.hpp".

4) In NewAnalysis.cpp, look for ::analyze() and change this line into:  

           void NewAnalysis::analyze() {  
           
   Then make all the changes you want in NewAnalysis.hpp to create the desired analysis method.

5) In include/NewAnalysis.hpp, change the first two lines to  

           #ifndef NEW_ANALYSIS  
           #define NEW_ANALYSIS
           
6) Still in include/NewAnalysis.hpp, change the class name to:  

           class NewAnalysis  

7) Modify CMakeLists.txt: Look for the expression "set(XA_SRC" 
     and below this line, add to the list of files the new filename:  
     
         NewAnalysis.cpp  

8) Go to src/KWa/ and modify AnalysisInterface.cpp as follows:  
   add a line   
   
        #include "NewAnalysis.hpp"
       
   in the beginning. Then  look for the expression
      
        //add a new Analysis instance ...  
      
      and just before that line, add the following:

        } else if (analysisName == "NewAnalysis") {  
          analysis = new NewAnalysis(  
              analysisName, observableName, *minvalue, *maxvalue, *nrbins,  
              binningType, *normalization, *nbOfMoreparameters, *nbOfSubanalyses);

## Examples

In the current directory (srcext/XA/), one finds a couple of "typical" analysis methods, meant to be templates for the creation of new ones, being (hopefully)
sufficiently well documented:

          MultiplicityDistributions.cpp 
          PtSpectraMultTrigger.cpp
          ParticleSpectraCentralityTrigger.cpp
          ParticleSpectraMultiplicityTrigger.cpp
          ParticleSpectraPercentileTrigger.cpp

Examples of the corresponding analysis configurations (to be put in the optns files) can be found in the directory srcext/XA/optns/

