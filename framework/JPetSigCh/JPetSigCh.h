#ifndef _JPETSIGCH_H_
#define _JPETSIGCH_H_

#include <cassert>
#include <vector>
#include <map>
#include <TClass.h>

/* #include "../JPetPM/JPetPM.h" */
/* #include "../JPetBarrelSlot/JPetBarrelSlot.h" */
/* #include "../JPetScin/JPetScin.h" */
/* #include "../JPetTRB/JPetTRB.h" */
#include "../../JPetLoggerInclude.h"


class JPetSigCh: public TNamed
{
 public:
  enum EdgeType { kRising, kFalling, kSlow };
  /* typedef std::pair < float, float > Channels; */
  const static float kTimeUnset;

  /// @warning here I dont know who should be the owner of JPetTRB etc elements
  friend void my_swap(JPetSigCh& first, JPetSigCh& second) {
    using std::swap;
    swap(first.fValue, second.fValue);
    swap(first.fType, second.fType);
    swap(first.fPMID, second.fPMID);
    swap(first.fThreshold, second.fThreshold);

    /* swap(first.fChannels, second.fChannels); */
    /* swap(first.fIsSlow, second.fIsSlow); */
    /* swap(first.fIsComplete, second.fIsComplete); */
    /* swap(first.fPM, second.fPM); */
    /* swap(first.fScin, second.fScin); */
    /* swap(first.fBarrelSlot, second.fBarrelSlot); */
  }

  JPetSigCh() { init(); }
  JPetSigCh(const JPetSigCh& obj);
  JPetSigCh& operator= (const JPetSigCh obj);
  /* JPetSigCh(float EdgeTime, float FallEdgeTime); */
  JPetSigCh(EdgeType Edge, float EdgeTime);
  ~JPetSigCh() {}
  /* inline bool isComplete() const { return fIsComplete; } */
  inline float getValue() const { return fValue; }
  /**
   * @warning This method may cause seg fault, when is called with kFalling as first argument and object is of type "slow".
   */
  /* inline JPetPM getPM() const { return fPM; } */
  /* inline JPetTRB getTRB() const {return fTRB; } */
  /* inline JPetScin getScin() const { return fScin; } */
  /* inline JPetBarrelSlot getBarrelSlot() const { return fBarrelSlot; } */
  inline EdgeType getType() const { return fType; }

  bool isSlow() const ;
  /* float getTime(EdgeType type) const ; */
 
  /* 
  inline Channels getChannels() const { return fChannels; }
  void addCh(float rise_edge_time, float fall_edge_time);
  */
  /* inline void setPM(const JPetPM& pm) { fPM = pm; } */
  /* inline void setTRB(const JPetTRB& trb) { fTRB = trb; } */
  /* inline void setScin(const JPetScin& scin) { fScin = scin; } */
  /* inline void setBarrelSlot(const JPetBarrelSlot& barrel_slot) { fBarrelSlot = barrel_sl ot; } */
  inline void setValue( float val ) { fValue = val; }
  inline void setType( EdgeType type ) { fType = type; }

  inline void setPMID( Int_t pmid ) { fPMID = pmid; }
  inline void setThreshold( float thr ) { fThreshold = thr; }
  
  inline Int_t getPMID() const { return fPMID; }
  inline float getThreshold() const { return fThreshold; }

/*
  inline void setSlow( bool isSlow ) { fIsSlow = isSlow; }
  inline void setComplete( bool isComplete ) { fIsComplete = isComplete; }
*/	
  ClassDef(JPetSigCh, 1);

 protected:
  Int_t fPMID;
  float fThreshold;
  EdgeType fType;
  float fValue;
  /*
    Channels fChannels; /// fChannels.first is rising edge and fChannels.second is falling edge
    bool fIsSlow;
    bool fIsComplete;
  */
  /* JPetPM fPM; */
  /* JPetTRB fTRB; */
  /* JPetScin fScin; */
  /* JPetBarrelSlot fBarrelSlot; */

  /* template <class T> void set(T** dest, const T& source) throw(std::bad_alloc); */
  void init();
};

#endif
