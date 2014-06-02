// JPetPM.cpp - Photomultiplier
#include "JPetPM.h"


ClassImp(JPetPM);

JPetPM::JPetPM() :
		  fSide(kLeft),
		  fId(0),
		  fIsActive(false),
		  fStatus(std::string("")),
		  fName(std::string("")),
		  fMaxHV(0.0),
		  fDescription(std::string("")),
		  fProducer(std::string("")),
		  fBoughtDate(std::string("")),
		  fSerialNumber(std::string("")),
		  fTakeSpositiveVoltage(false),
		  fPMModel(0, std::string("")),
		  fPMCalibration(0, 
				 std::string(""), 
				 0.f,
				 0.f, 
				 0.f, 
				 0.f, 
				 0.f),
		  fUser()
{}

JPetPM::JPetPM(JPetPM::Side p_side,
	       int p_id,
	       bool p_isActive, 
	       std::string p_status, 
	       std::string p_name, 
	       double p_maxHV, 
	       std::string p_description, 
	       std::string p_producer, 
	       std::string p_boughtDate,
	       std::string p_serialNumber,
	       bool p_takeSpositiveVoltage,
	       int p_modelId,
	       std::string p_modelName,
	       int p_calibrationId,
	       std::string p_calibrationName,
	       float p_calibrationOptHV,
	       float p_calibrationC2e_1,
	       float p_calibrationC2e_2,
	       float p_calibrationGainAlpha,
	       float p_calibrationGainBeta,
	       JPetUser p_user) :
				  fSide(p_side),
				  fId(p_id),
				  fIsActive(p_isActive),
				  fStatus(p_status),
				  fName(p_name),
				  fMaxHV(p_maxHV),
				  fDescription(p_description),
				  fProducer(p_producer),
				  fBoughtDate(p_boughtDate),
				  fSerialNumber(p_serialNumber),
				  fTakeSpositiveVoltage(p_takeSpositiveVoltage),
				  fPMModel(p_modelId, p_modelName),
				  fPMCalibration(p_calibrationId, 
						 p_calibrationName, 
						 p_calibrationOptHV,
						 p_calibrationC2e_1,
						 p_calibrationC2e_2, 
						 p_calibrationGainAlpha, 
						 p_calibrationGainBeta),
				  fUser(p_user)
{}

JPetPM::~JPetPM()
{}
