// JPetTRB.cpp - TRB
#include "JPetTRB.h"


ClassImp(JPetTRB);

JPetTRB::JPetTRB() : 
		    fId(0),
		    fIsActive(false),
		    fStatus(std::string("")),
		    fMaxch(0),
		    fName(std::string("")),
		    fDescription(std::string("")),
		    fVersion(0)
{}

JPetTRB::JPetTRB(int p_id,
		 bool p_isActive, 
		 std::string p_status, 
		 int p_maxch, 
		 std::string p_name, 
		 std::string p_description, 
		 int p_version,
		 const JPetUser &p_user) :
					  fId(p_id),
					  fIsActive(p_isActive),
					  fStatus(p_status),
					  fMaxch(p_maxch),
					  fName(p_name),
					  fDescription(p_description),
					  fVersion(p_version),
					  fUser(p_user)
{}

JPetTRB::~JPetTRB()
{}
