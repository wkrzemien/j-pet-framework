// JPetReader.h - Reader
#ifndef JPETREADER_H 
#define JPETREADER_H 

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>

#ifndef __CINT__
#include <boost/noncopyable.hpp>
#else
namespace boost;
class boost::noncopyable;
#endif /* __CINT __ */

/*
#include "../JPetScin/JPetScin.h"
#include "../JPetPM/JPetPM.h"
#include "../JPetFEB/JPetFEB.h"
#include "../JPetTRB/JPetTRB.h"
*/
#include "../JPetTreeHeader/JPetTreeHeader.h"

#include "../../JPetLoggerInclude.h"

/**
 * @brief A class responsible for reading any data from ROOT trees.
 *
 * All objects inheriting from JPetAnalysisModule should use this class in order to access and read data from ROOT files.
 */
class JPetReader : private boost::noncopyable
{
public:
  JPetReader(void);
  explicit JPetReader(const char* p_filename);
  virtual ~JPetReader(void);
  
  virtual void CloseFile();
  virtual long long GetEntries () const { return fTree->GetEntries(); }
  virtual int GetEntry (int entryNo) {return fTree->GetEntry(entryNo); } /// the name of the function is bad but it mimics ROOT function 
  virtual bool OpenFile(const char* filename);
  virtual void ReadData(const char* objname = "tree");
  virtual TNamed& GetData () {return *fObject;}
  JPetTreeHeader * GetHeaderClone()const;
  virtual TObject * GetObject(const char * name) const {return fFile->Get(name);}
  virtual bool isOpen() const {if (fFile) return fFile->IsOpen(); 
                               else return false;}
  
  template <class T>
  void fillContainer(std::vector<T> &p_container, const std::string &p_objectName);
  virtual void closeTFile(void){};
  
protected:
  TBranch* fBranch;
  TNamed* fObject;
  TTree* fTree;
  TFile* fFile;
  
};


template <class T>
void JPetReader::fillContainer(std::vector<T> &p_container, const std::string &p_objectName)
{
  TList *l_TList = (TList*)fFile->Get(p_objectName.c_str());
  TObject *l_obj;
  
  TIter next(l_TList);
  while(l_obj = next())
  {
    T *l_item = static_cast<T*>(l_obj);
    p_container.push_back(*l_item);
  }
}

#endif	// JPETREADER_H
