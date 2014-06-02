// JPetScin.cpp - Scintillator
#include "JPetScin.h"
#include <cassert>


ClassImp(JPetScin);

JPetScin::JPetScin() :
		      fId(0),
		      fIsActive(false),
		      fStatus(std::string("")),
		      fScinDimensions(0.f, 0.f, 0.f),
		      fDescription(std::string("")),
		      fScinType(0, std::string(""), std::string("")),
		      fScinCalibration(0, std::string(""), 0.f),
		      fUser()
{}

JPetScin::JPetScin(int p_id, 
		   bool p_isActive,
		   std::string p_status, 
		   float p_length, 
		   float p_height, 
		   float p_width,
		   std::string p_description, 
		   int p_typeId, 
		   std::string p_typeName, 
		   std::string p_typeDescription,
		   int p_calibrationId, 
		   std::string p_calibrationName, 
		   float p_calibrationAttlength,
		   const JPetUser &p_user) :
					    fId(p_id),
					    fIsActive(p_isActive),
					    fStatus(p_status),
					    fScinDimensions(p_length, p_height, p_width),
					    fDescription(p_description),
					    fScinType(p_typeId, p_typeName, p_typeDescription),
					    fScinCalibration(p_calibrationId, p_calibrationName, p_calibrationAttlength),
					    fUser(p_user)
{}

JPetScin::~JPetScin()
{}

float JPetScin::getScinDimension(const JPetScin::Dimension &p_dimension) const
{
  switch (p_dimension) 
  {
  case kLength:
    return fScinDimensions.fLength;
  case kHeight:
    return fScinDimensions.fHeight;
  case kWidth:
    return fScinDimensions.fWidth;
  default:
    assert(1 == 0);	//do przemyslenia
  }
}

void JPetScin::setScinDimension(const JPetScin::Dimension &p_dimension, const float p_value)
{
  switch (p_dimension) 
  {
  case kLength:
    fScinDimensions.fLength = p_value;
    break;
  case kHeight:
    fScinDimensions.fHeight = p_value;
    break;
  case kWidth:
    fScinDimensions.fWidth = p_value;
    break;
  default:
    assert(1 == 0);	//do przemyslenia
  }
}
