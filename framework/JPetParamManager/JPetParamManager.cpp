// JPetParamManager.cpp - JPetParamManager
#include "JPetParamManager.h"
//#include "../DBHandler/HeaderFiles/ParamServer.h"
#include "../DBHandler/HeaderFiles/DBHandler.h"
#include "../CommonTools/CommonTools.h"
#include "../JPetReader/JPetReader.h"
#include "../JPetWriter/JPetWriter.h"
#include <boost/lexical_cast.hpp>
#include <sstream>


JPetParamManager::JPetParamManager() :
				      fScintilatorsObjectName(std::string("")),
				      fPMsObjectName(std::string("")),
				      fKBsObjectName(std::string("")),
				      fTRBsObjectName(std::string("")),
				      fTOMBsObjectName(std::string(""))
{}

JPetParamManager::JPetParamManager(const std::string &p_scintilatorsObjectName,
				    const std::string &p_PMsObjectName,
				    const std::string &p_KBsObjectName,
				    const std::string &p_TRBsObjectName,
				    const std::string &p_TOMBsObjectName) :
									    fScintilatorsObjectName(p_scintilatorsObjectName),
									    fPMsObjectName(p_PMsObjectName),
									    fKBsObjectName(p_KBsObjectName),
									    fTRBsObjectName(p_TRBsObjectName),
									    fTOMBsObjectName(p_TOMBsObjectName)
{}

JPetParamManager::~JPetParamManager()
{}

void JPetParamManager::readFile(const char * file_name){
	assert(file_name != NULL);
	
	std::string line;
	
	ifstream file(file_name);
	if( ! file.is_open() ){
		ERROR("No such file!");
		return;
	}
	
	while ( !file.eof() ){
		int trb = 0, scin = 0;
		
		getline(file, line);
		std::istringstream iss(line);
		
		iss >> trb;
		iss >> scin;
		//fTRBNumbers.push_back(trb);
		//fScinNumbers.push_back(scin);
	}
}


void JPetParamManager::fillScintillators(const int p_run_id)
{
  DB::SERVICES::DBHandler &l_dbHandlerInstance = DB::SERVICES::DBHandler::getInstance();
  
  std::string l_run_id = boost::lexical_cast<std::string>(p_run_id);
  
  std::string l_sqlQuerry = "SELECT * FROM getScintillatorsData(" + l_run_id + ");";
  pqxx::result l_runDbResults = l_dbHandlerInstance.querry(l_sqlQuerry);

  size_t l_sizeResultQuerry = l_runDbResults.size();

  if(l_sizeResultQuerry)
  {
    for(pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row)
    {
      int l_scintillator_id = row["scintillator_id"].as<int>();
      bool l_scintillator_isactive = row["scintillator_isactive"].as<bool>();
      std::string l_scintillator_status = row["scintillator_status"].as<std::string>();
      
      double l_scintillator_length = row["scintillator_length"].as<double>();
      double l_scintillator_width = row["scintillator_width"].as<double>();
      double l_scintillator_height = row["scintillator_height"].as<double>();
      std::string l_scintillator_description = row["scintillator_description"].as<std::string>();

      int l_scintillatortype_id = row["scintillatortype_id"].as<int>();
      std::string l_scintillatortype_name = row["scintillatortype_name"].as<std::string>();
      std::string l_scintillatortype_description = row["scintillatortype_description"].as<std::string>();
      
      int l_scintillatorcalibration_id = row["scintillatorcalibration_id"].as<int>();
      std::string l_scintillatorcalibration_name = row["scintillatorcalibration_name"].as<std::string>();
      float l_scintillatorcalibration_attlength = row["scintillatorcalibration_attlength"].as<float>();
      
      int l_user_id = row["user_id"].as<int>();
      std::string l_user_name = row["user_name"].as<std::string>();
      std::string l_user_lastName = row["user_lastName"].as<std::string>();
      std::string l_user_login = row["user_login"].as<std::string>();
      std::string l_user_password = row["user_password"].as<std::string>();
      bool l_user_isroot = row["user_isroot"].as<bool>();
      std::string l_user_creationdate = row["user_creationdate"].as<std::string>();
      std::string l_user_lastlogindate = row["user_lastlogindate"].as<std::string>();
      
      int l_setup_id = row["setup_id"].as<int>();
      int l_run_id = row["run_id"].as<int>();
      
      JPetUser l_user(
		      l_user_id,
		      l_user_name,
		      l_user_lastName,
		      l_user_login,
		      l_user_password,
		      l_user_isroot,
		      l_user_creationdate,
		      l_user_lastlogindate);
      
      JPetScin l_scin(
		      l_scintillator_id,
		      l_scintillator_isactive,
		      l_scintillator_status,

		      l_scintillator_length,
		      l_scintillator_width,
		      l_scintillator_height,
		      l_scintillator_description,

		      l_scintillatortype_id,
		      l_scintillatortype_name,
		      l_scintillatortype_description,

		      l_scintillatorcalibration_id,
		      l_scintillatorcalibration_name,
		      l_scintillatorcalibration_attlength,

		      l_user);
      
      fScintillators.push_back(l_scin);
    }
  }
}

void JPetParamManager::fillPMs(const int p_run_id)
{
  DB::SERVICES::DBHandler &l_dbHandlerInstance = DB::SERVICES::DBHandler::getInstance();
  
  std::string l_run_id = boost::lexical_cast<std::string>(p_run_id);
  
  std::string l_sqlQuerry = "SELECT * FROM getPhotoMultipliersData(" + l_run_id + ");";
  pqxx::result l_runDbResults = l_dbHandlerInstance.querry(l_sqlQuerry);

  size_t l_sizeResultQuerry = l_runDbResults.size();

  if(l_sizeResultQuerry)
  {
    for(pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row)
    {
      int l_hvpmconnection_id = row["hvpmconnection_id"].as<int>();
      bool l_hvpmconnection_isrightside = row["hvpmconnection_isrightside"].as<bool>();
	    
      int l_photomultiplier_id = row["photomultiplier_id"].as<int>();
      bool l_photomultiplier_isactive = row["photomultiplier_isactive"].as<bool>();
      std::string l_photomultiplier_status = row["photomultiplier_status"].as<std::string>();   
      std::string l_photomultiplier_name = row["photomultiplier_name"].as<std::string>();
      double l_photomultiplier_maxhv = row["photomultiplier_maxhv"].as<double>();
      std::string l_photomultiplier_description = row["photomultiplier_description"].as<std::string>();
      std::string l_photomultiplier_producer = row["photomultiplier_producer"].as<std::string>();
      std::string l_photomultiplier_boughtdate = row["photomultiplier_boughtdate"].as<std::string>();
      std::string l_photomultiplier_serialnumber = row["photomultiplier_serialnumber"].as<std::string>();
      bool l_photomultiplier_takespositivevoltage = row["photomultiplier_takespositivevoltage"].as<bool>();
	    
      int l_pmmodel_id = row["pmmodel_id"].as<int>();
      std::string l_pmmodel_name = row["pmmodel_name"].as<std::string>();

      int l_pmcalibration_id = row["pmcalibration_id"].as<int>();
      std::string l_pmcalibration_name = row["pmcalibration_name"].as<std::string>();
      double l_pmcalibration_opthv = row["pmcalibration_opthv"].as<double>();
      double l_pmcalibration_c2e_1 = row["pmcalibration_c2e_1"].as<double>();
      double l_pmcalibration_c2e_2 = row["pmcalibration_c2e_2"].as<double>();
      double l_pmcalibration_gainalpha = row["pmcalibration_gainalpha"].as<double>();
      double l_pmcalibration_gainbeta = row["pmcalibration_gainbeta"].as<double>();
      
      int l_user_id = row["user_id"].as<int>();
      std::string l_user_name = row["user_name"].as<std::string>();
      std::string l_user_lastName = row["user_lastName"].as<std::string>();
      std::string l_user_login = row["user_login"].as<std::string>();
      std::string l_user_password = row["user_password"].as<std::string>();
      bool l_user_isroot = row["user_isroot"].as<bool>();
      std::string l_user_creationdate = row["user_creationdate"].as<std::string>();
      std::string l_user_lastlogindate = row["user_lastlogindate"].as<std::string>();
      
      int l_setup_id = row["setup_id"].as<int>();
      int l_run_id = row["run_id"].as<int>();
      
      JPetUser l_user(
		      l_user_id,
		      l_user_name,
		      l_user_lastName,
		      l_user_login,
		      l_user_password,
		      l_user_isroot,
		      l_user_creationdate,
		      l_user_lastlogindate);
      
      JPetPM::Side l_side = (l_hvpmconnection_isrightside) ? JPetPM::Side::kRight : JPetPM::Side::kLeft;
      
      JPetPM l_pm(
		  l_side,
		  l_photomultiplier_id,
		  l_photomultiplier_isactive,
		  l_photomultiplier_status, 
		  l_photomultiplier_name,
		  l_photomultiplier_maxhv,
		  l_photomultiplier_description,
		  l_photomultiplier_producer,
		  l_photomultiplier_boughtdate,
		  l_photomultiplier_serialnumber,
		  l_photomultiplier_takespositivevoltage,
			
		  l_pmmodel_id,
		  l_pmmodel_name,

		  l_pmcalibration_id,
		  l_pmcalibration_name,
		  l_pmcalibration_opthv,
		  l_pmcalibration_c2e_1,
		  l_pmcalibration_c2e_2,
		  l_pmcalibration_gainalpha,
		  l_pmcalibration_gainbeta,
		  
		  l_user);
      
      fPMs.push_back(l_pm);
    }
  }
}

void JPetParamManager::fillKBs(const int p_run_id)
{
  DB::SERVICES::DBHandler &l_dbHandlerInstance = DB::SERVICES::DBHandler::getInstance();
  
  std::string l_run_id = boost::lexical_cast<std::string>(p_run_id);
  
  std::string l_sqlQuerry = "SELECT * FROM getKonradBoardsData(" + l_run_id + ");";
  pqxx::result l_runDbResults = l_dbHandlerInstance.querry(l_sqlQuerry);

  size_t l_sizeResultQuerry = l_runDbResults.size();

  if(l_sizeResultQuerry)
  {
    for(pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row)
    {
      int l_konradboard_id = row["konradboard_id"].as<int>();
      bool l_konradboard_isactive = row["konradboard_isactive"].as<bool>();
      std::string l_konradboard_status = row["konradboard_status"].as<std::string>();
      std::string l_konradboard_description = row["konradboard_description"].as<std::string>();
      int l_konradboard_version = row["konradboard_version"].as<int>();
      int l_konradboard_creator_id = row["konradboard_creator_id"].as<int>();
      
      int l_user_id = row["user_id"].as<int>();
      std::string l_user_name = row["user_name"].as<std::string>();
      std::string l_user_lastName = row["user_lastName"].as<std::string>();
      std::string l_user_login = row["user_login"].as<std::string>();
      std::string l_user_password = row["user_password"].as<std::string>();
      bool l_user_isroot = row["user_isroot"].as<bool>();
      std::string l_user_creationdate = row["user_creationdate"].as<std::string>();
      std::string l_user_lastlogindate = row["user_lastlogindate"].as<std::string>();
      
      int l_setup_id = row["setup_id"].as<int>();
      int l_run_id = row["run_id"].as<int>();
      
      JPetUser l_user(
		      l_user_id,
		      l_user_name,
		      l_user_lastName,
		      l_user_login,
		      l_user_password,
		      l_user_isroot,
		      l_user_creationdate,
		      l_user_lastlogindate);

      JPetKB l_KB(
		  l_konradboard_id, 
		  l_konradboard_isactive, 
		  l_konradboard_status, 
		  l_konradboard_description, 
		  l_konradboard_version, 
		  l_user
		 );
      
      fKBs.push_back(l_KB);
    }
  }
}

void JPetParamManager::fillTRBs(const int p_run_id)
{
  DB::SERVICES::DBHandler &l_dbHandlerInstance = DB::SERVICES::DBHandler::getInstance();
  
  std::string l_run_id = boost::lexical_cast<std::string>(p_run_id);
  
  std::string l_sqlQuerry = "SELECT * FROM getTRBsData(" + l_run_id + ");";
  pqxx::result l_runDbResults = l_dbHandlerInstance.querry(l_sqlQuerry);

  size_t l_sizeResultQuerry = l_runDbResults.size();
  
  if(l_sizeResultQuerry)
  {
    for(pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row)
    {
      int l_TRB_id = row["TRB_id"].as<int>();
      bool l_TRB_isactive = row["TRB_isactive"].as<bool>();
      std::string l_TRB_status = row["TRB_status"].as<std::string>();
      int l_TRB_maxch = row["TRB_maxch"].as<int>();
      std::string l_TRB_name = row["TRB_name"].as<std::string>();
      std::string l_TRB_description = row["TRB_description"].as<std::string>();
      int l_TRB_version = row["TRB_version"].as<int>();
      
      int l_user_id = row["user_id"].as<int>();
      std::string l_user_name = row["user_name"].as<std::string>();
      std::string l_user_lastName = row["user_lastName"].as<std::string>();
      std::string l_user_login = row["user_login"].as<std::string>();
      std::string l_user_password = row["user_password"].as<std::string>();
      bool l_user_isroot = row["user_isroot"].as<bool>();
      std::string l_user_creationdate = row["user_creationdate"].as<std::string>();
      std::string l_user_lastlogindate = row["user_lastlogindate"].as<std::string>();
      
      int l_setup_id = row["setup_id"].as<int>();
      int l_run_id = row["run_id"].as<int>();
      
      JPetUser l_user(
		      l_user_id,
		      l_user_name,
		      l_user_lastName,
		      l_user_login,
		      l_user_password,
		      l_user_isroot,
		      l_user_creationdate,
		      l_user_lastlogindate);

      JPetTRB l_trb(
		    l_TRB_id,
		    l_TRB_isactive,
		    l_TRB_status,
		    l_TRB_maxch,
		    l_TRB_name,
		    l_TRB_description,
		    l_TRB_version,
		    l_user);
      
      fTRBs.push_back(l_trb);
    }
  }
}

void JPetParamManager::fillTOMBs(const int p_run_id)
{
  DB::SERVICES::DBHandler &l_dbHandlerInstance = DB::SERVICES::DBHandler::getInstance();
  
  std::string l_run_id = boost::lexical_cast<std::string>(p_run_id);
  
  std::string l_sqlQuerry = "SELECT * FROM getTOMBData(" + l_run_id + ");";
  pqxx::result l_runDbResults = l_dbHandlerInstance.querry(l_sqlQuerry);

  size_t l_sizeResultQuerry = l_runDbResults.size();
  
  if(l_sizeResultQuerry)
  {
    for(pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row)
    {
      int l_TOMB_id = row["TOMB_id"].as<int>();
      std::string l_TOMB_description = row["TOMB_description"].as<std::string>();
      int l_TOMB_setup_id = row["TOMB_setup_id"].as<int>();
      
      int l_setup_id = row["setup_id"].as<int>();
      int l_run_id = row["run_id"].as<int>();

      JPetTOMB l_tomb(
		      l_TOMB_id,
		      l_TOMB_description,
		      l_TOMB_setup_id);
      
      fTOMBs.push_back(l_tomb);
    }
  }
}

void JPetParamManager::fillContainers(const int p_run_id)
{ 
  fillScintillators(p_run_id);
  fillPMs(p_run_id);
  fillKBs(p_run_id);
  fillTRBs(p_run_id);
  fillTOMBs(p_run_id);
}

void JPetParamManager::fillContainers(const JPetParamManager::ContainerType &p_containerType, const char *p_fileName)
{
  switch(p_containerType) 
  {
    case kScintillator:
    {
      fillScintillators(p_fileName);
      break;
    }
    case kPM:
    {
      fillPMs(p_fileName);
      break;
    }
    case kKB:
    {
      fillKBs(p_fileName);
      break;
    }
    case kTRB:
    {
      fillTRBs(p_fileName);
      break;
    }
    case kTOMB:
    {
      fillTOMBs(p_fileName);
      break;
    }
  }
}

void JPetParamManager::generateRootFile(const JPetParamManager::ContainerType &p_containerType, const char *p_fileName)
{
  switch(p_containerType) 
  {
    case kScintillator:
    {
      generateRootFileForScintillators(p_fileName);
      break;
    }
    case kPM:
    {
      generateRootFileForPMs(p_fileName);
      break;
    }
    case kKB:
    {
      generateRootFileForKBs(p_fileName);
      break;
    }
    case kTRB:
    {
      generateRootFileForTRBs(p_fileName);
      break;
    }
    case kTOMB:
    {
      generateRootFileForTOMBs(p_fileName);
      break;
    }
  }
}

void JPetParamManager::fillScintillators(const char *p_fileName)
{
  std::string l_fileToRead = p_fileName;
  
  if(CommonTools::findSubstring(l_fileToRead, std::string(".root")) == std::string::npos)
  {
    l_fileToRead.append(".root");
  }
  
  JPetReader l_reader(l_fileToRead.c_str());
  l_reader.fillContainer(fScintillators, fScintilatorsObjectName);
  l_reader.closeTFile();
}

void JPetParamManager::fillPMs(const char *p_fileName)
{
  std::string l_fileToRead = p_fileName;
  
  if(CommonTools::findSubstring(l_fileToRead, std::string(".root")) == std::string::npos)
  {
    l_fileToRead.append(".root");
  }
  
  JPetReader l_reader(l_fileToRead.c_str());
  l_reader.fillContainer(fPMs, fPMsObjectName);
  l_reader.closeTFile();
}

void JPetParamManager::fillKBs(const char *p_fileName)
{
  std::string l_fileToRead = p_fileName;
  
  if(CommonTools::findSubstring(l_fileToRead, std::string(".root")) == std::string::npos)
  {
    l_fileToRead.append(".root");
  }
  
  JPetReader l_reader(l_fileToRead.c_str());
  l_reader.fillContainer(fKBs, fTRBsObjectName);
  //std::cout << "l_kb->isActive() = " << fKBs[0].isActive() << std::endl;
  l_reader.closeTFile();
}

void JPetParamManager::fillTRBs(const char *p_fileName)
{
  std::string l_fileToRead = p_fileName;
  
  if(CommonTools::findSubstring(l_fileToRead, std::string(".root")) == std::string::npos)
  {
    l_fileToRead.append(".root");
  }
  
  JPetReader l_reader(l_fileToRead.c_str());
  l_reader.fillContainer(fTRBs, fTRBsObjectName);
  l_reader.closeTFile();
}

void JPetParamManager::fillTOMBs(const char *p_fileName)
{
  std::string l_fileToRead = p_fileName;
  
  if(CommonTools::findSubstring(l_fileToRead, std::string(".root")) == std::string::npos)
  {
    l_fileToRead.append(".root");
  }
  
  JPetReader l_reader(l_fileToRead.c_str());
  l_reader.fillContainer(fTOMBs, fTOMBsObjectName);
  l_reader.closeTFile();
}

void JPetParamManager::generateRootFileForScintillators(const char *p_fileName)
{
  if(!fScintillators.empty())
  {
    std::string l_fileToSave = p_fileName;
    if(CommonTools::findSubstring(l_fileToSave, std::string(".root")) == std::string::npos)
    {
      l_fileToSave.append(".root");
    }
  
    JPetWriter l_writer(l_fileToSave.c_str());
    l_writer.write(fScintillators, fScintilatorsObjectName);
    l_writer.closeTFile();
  }
}

void JPetParamManager::generateRootFileForPMs(const char *p_fileName)
{
  if(!fPMs.empty())
  {
    std::string l_fileToSave = p_fileName;
    if(CommonTools::findSubstring(l_fileToSave, std::string(".root")) == std::string::npos)
    {
      l_fileToSave.append(".root");
    }
  
    JPetWriter l_writer(l_fileToSave.c_str());
    l_writer.write(fPMs, fPMsObjectName);
    l_writer.closeTFile();
  }
}

void JPetParamManager::generateRootFileForKBs(const char *p_fileName)
{
  if(!fKBs.empty())
  {
    std::string l_fileToSave = p_fileName;
    if(CommonTools::findSubstring(l_fileToSave, std::string(".root")) == std::string::npos)
    {
      l_fileToSave.append(".root");
    }
  
    JPetWriter l_writer(l_fileToSave.c_str());
    l_writer.write(fKBs, fKBsObjectName);
    l_writer.closeTFile();
  }
}

void JPetParamManager::generateRootFileForTRBs(const char *p_fileName)
{
  if(!fTRBs.empty())
  {
    std::string l_fileToSave = p_fileName;
    if(CommonTools::findSubstring(l_fileToSave, std::string(".root")) == std::string::npos)
    {
      l_fileToSave.append(".root");
    }
  
    JPetWriter l_writer(l_fileToSave.c_str());
    l_writer.write(fTRBs, fTRBsObjectName);
    l_writer.closeTFile();
  }
}

void JPetParamManager::generateRootFileForTOMBs(const char *p_fileName)
{
  if(!fTOMBs.empty())
  {
    std::string l_fileToSave = p_fileName;
    if(CommonTools::findSubstring(l_fileToSave, std::string(".root")) == std::string::npos)
    {
      l_fileToSave.append(".root");
    }
  
    JPetWriter l_writer(l_fileToSave.c_str());
    l_writer.write(fTOMBs, fTOMBsObjectName);
    l_writer.closeTFile();
  }
}

template <typename T>
std::vector<T>& JPetParamManager::getContainer(const JPetParamManager::ContainerType &p_containerType)
{
  switch(p_containerType) 
  {
  case kScintillator:
    return fScintillators;
  case kPM:
    return fPMs;
  case kKB:
    return fKBs;
  case kTRB:
    return fTRBs;
  case kTOMB:
    return fTOMBs;
    
  default:
    assert(1 == 0);	//do przemyslenia
  }
}

template <typename T>
bool JPetParamManager::getData(const JPetParamManager::ContainerType &p_containerType, unsigned int p_index, T &p_data)
{
  switch(p_containerType) 
  {
    case kScintillator:
    {
      if(p_index < fScintillators.size())
      {
	p_data = fScintillators[p_index];
	return true;
      }
    }
    case kPM:
    {
      if(p_index < fPMs.size())
      {
	p_data = fPMs[p_index];
	return true;
      }
    }
    case kKB:
    {
      if(p_index < fKBs.size())
      {
	p_data = fKBs[p_index];
	return true;
      }
    }
    case kTRB:
    {
      if(p_index < fTRBs.size())
      {
	p_data = fTRBs[p_index];
	return true;
      }
    }
    case kTOMB:
    {
      if(p_index < fTOMBs.size())
      {
	p_data = fTOMBs[p_index];
	return true;
      }
    }
  }
  
  return false;
}

int JPetParamManager::getDataSize(const JPetParamManager::ContainerType &p_containerType) const
{
  switch(p_containerType) 
  {
  case kScintillator:
    return fScintillators.size();
  case kPM:
    return fPMs.size();
  case kKB:
    return fKBs.size();
  case kTRB:
    return fTRBs.size();
  case kTOMB:
    return fTOMBs.size();
    
  default:
    assert(1 == 0);	//do przemyslenia
  }
}
