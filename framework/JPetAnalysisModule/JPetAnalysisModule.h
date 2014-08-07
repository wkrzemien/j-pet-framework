/**
 *  @copyright Copyright (c) 2013, Wojciech Krzemien
 *  @file JPetAnalysisModule.h 
 *  @author Wojciech Krzemien, wojciech.krzemien@if.uj.edu.pl
 *  @brief
 */ 

#ifndef JPETANALYSISMODULE_H 
#define JPETANALYSISMODULE_H

#include <vector>
#include <TNamed.h>
#include <TTree.h>
#include <TList.h>

class JPetAnalysisModule: public TNamed {
 public:
  JPetAnalysisModule();
  JPetAnalysisModule(const char* name, const char* title, TTree * shared_tree = NULL);
  virtual ~JPetAnalysisModule(); 
  virtual void CreateInputObjects(const char* inputFilename=0)=0; //
  virtual void CreateOutputObjects(const char* outputFilename=0)=0; //
  virtual void Exec()=0; // called for every event
  virtual long long GetEventNb()=0;
  virtual void RunSubmodules();
  virtual void Terminate()=0; // called once when analysis terminates

  int AddStatsObject(TObject * statObj);
  const TList * GetStatsObjects() const; 

  ClassDef(JPetAnalysisModule,1);

protected:
  virtual void AddSubmodule( JPetAnalysisModule* new_submodule );
  TTree fSubmoduleSharedTree;
  TTree* fSuperSharedTree;
  std::vector< JPetAnalysisModule* > fSubmodules;
  TList fStats; ///< a list to store all objects for statistics of the processing, i.e. histograms
};
#endif /*  !JPETANALYSISMODULE_H */
