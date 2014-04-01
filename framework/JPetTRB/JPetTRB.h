#ifndef _JPETTRB_H_
#define _JPETTRB_H_

#include "TNamed.h"

class JPetTRB: public TNamed
{
 public:
  JPetTRB();
  JPetTRB(int id, int type, int channel);
  inline int getID() const { return fID; }
  inline int getType() const { return fType; }
  inline int getChannel() const { return fChannel; }
  inline void setID(int id) { fID = id; }
  inline void setType(int type) { fType = type; }
  inline void setChannel(int ch) { fChannel = ch; }

 private:
  int fID;
  int fType;
  int fChannel;
  /// @todo do implementacji
  //JPetKB* KBId;
  //KBType;
  //KBChan;
  //
  ClassDef(JPetTRB, 1);
};

#endif
