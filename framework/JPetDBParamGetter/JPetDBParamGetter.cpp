/**
  *  @copyright Copyright (c) 2014, Wojciech Krzemien
  *  @file JPetDBParamGetter.cpp
  *  @author Wojciech Krzemien, wojciech.krzemien@if.uj.edu.pl
  */ 

#include "./JPetDBParamGetter.h"
#include <boost/lexical_cast.hpp>
#include "../DBHandler/HeaderFiles/DBHandler.h"

JPetDBParamGetter::JPetDBParamGetter()
{
}

/// @param DBConfigFile configuration file with the database connection settings
JPetDBParamGetter::JPetDBParamGetter(const char* dBConfigFile)
{
  DB::SERVICES::DBHandler::getInstance(dBConfigFile); /// this command aims to load the configuration file, the return value is irrelevant

}

JPetParamBank JPetDBParamGetter::generateParamBank(const int p_run_id) 
{
  JPetParamBank paramBank;
  fillScintillators(p_run_id, paramBank);
  fillPMs(p_run_id, paramBank);
  fillFEBs(p_run_id, paramBank);
  fillTRBs(p_run_id, paramBank);
  fillTOMB(p_run_id, paramBank);
  fillAllTRefs(p_run_id, paramBank);
  return paramBank;
}

//void JPetDBParamGetter::fillParamContainer(ParamObjectType type, const int p_run_id, JPetParamBank& paramBank)
//{
//
//  // enum ParamObjectType {kScintillator, kPM, kFEB, kTRB, kTOMB, SIZE};
//  static const char* sqlFunctions[] = {"getDataFromScintillators", "getDataFromPhotoMultipliers","getDataFromKonradBoards", "getDataFromTRBs" ,"getDataFromTOMB"};
//
//  INFO("Start filling Scintillators container.");
//    //"SELECT * FROM getDataFromScintillators(" + l_run_id + ");";
//  pqxx::result l_runDbResults = getDataFromDB("getDataFromScintillators",p_run_id);
//
//
//  size_t l_sizeResultQuerry = l_runDbResults.size();
//
//  if (l_sizeResultQuerry) {
//    for (pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row) {
//      JPetScin l_scin = generateScintillator(row);
//      paramBank.addScintillator(l_scin);
//    }
//  } else {
//    printErrorMessageDB(sqlFunctions[type], p_run_id);
//  }
//
//}


void JPetDBParamGetter::fillScintillators(const int p_run_id, JPetParamBank& paramBank)
{
  INFO("Start filling Scintillators container.");
    //"SELECT * FROM getDataFromScintillators(" + l_run_id + ");";
  std::string l_run_id = boost::lexical_cast<std::string>(p_run_id);
  pqxx::result l_runDbResults = getDataFromDB("getDataFromScintillators",l_run_id);


  size_t l_sizeResultQuerry = l_runDbResults.size();

  if (l_sizeResultQuerry) {
    for (pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row) {
      JPetScin l_scin = generateScintillator(row);
      paramBank.addScintillator(l_scin);
    }
  } else {
    printErrorMessageDB("getDataFromScintillators", p_run_id);
  }
}

std::string JPetDBParamGetter::generateSelectQuery(const std::string& sqlFun, const std::string& arguments) 
{
  std::string sqlQuerry = "SELECT * FROM ";
  sqlQuerry += sqlFun;
  sqlQuerry += "(" +arguments + ");";
  return sqlQuerry;
}

/// @brief method calls the remote PostgreSQL function sqlfunction with the id argument and returns results from database
pqxx::result JPetDBParamGetter::getDataFromDB(const std::string& sqlfunction,const  std::string& arguments) 
{
  //std::string l_run_id = boost::lexical_cast<std::string>(id);
  std::string l_sqlQuerry = generateSelectQuery(sqlfunction,arguments);
  DB::SERVICES::DBHandler& l_dbHandlerInstance = DB::SERVICES::DBHandler::getInstance();
  return  l_dbHandlerInstance.querry(l_sqlQuerry);
}


void JPetDBParamGetter::printErrorMessageDB(std::string sqlFunction, int p_run_id) {
    std::string l_error(sqlFunction);
    l_error += "() querry for run_id = ";
    l_error += p_run_id + " return 0 records.";
    ERROR(l_error.c_str());
}



void JPetDBParamGetter::fillPMs(const int p_run_id, JPetParamBank& paramBank)
{
  INFO("Start filling PMs container.");
//  std::string l_sqlQuerry = "SELECT * FROM getDataFromPhotoMultipliers(" + l_run_id + ");";

  std::string l_run_id = boost::lexical_cast<std::string>(p_run_id);
  pqxx::result l_runDbResults = getDataFromDB("getDataFromPhotoMultipliers",l_run_id);

  size_t l_sizeResultQuerry = l_runDbResults.size();

  if (l_sizeResultQuerry) {
    for (pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row) {
      JPetPM l_pm = generatePM(row);
      paramBank.addPM(l_pm);
    }
  } else {
    printErrorMessageDB("getDataFromPhotoMultipliers", p_run_id);
  }
}

void JPetDBParamGetter::fillFEBs(const int p_run_id, JPetParamBank& paramBank)
{
  INFO("Start filling FEBs container.");

  //std::string l_sqlQuerry = "SELECT * FROM getDataFromKonradBoards(" + l_run_id + ");";
  std::string l_run_id = boost::lexical_cast<std::string>(p_run_id);
  pqxx::result l_runDbResults = getDataFromDB("getDataFromKonradBoards",l_run_id);

  size_t l_sizeResultQuerry = l_runDbResults.size();

  if (l_sizeResultQuerry) {
    for (pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row) {
      JPetFEB l_FEB = generateFEB(row);
      paramBank.addFEB(l_FEB);
    }
  } else {
    printErrorMessageDB("getDataFromKonradBoards", p_run_id);
  }
}

void JPetDBParamGetter::fillTRBs(const int p_run_id, JPetParamBank& paramBank)
{
  INFO("Start filling TRBs container.");


//  std::string l_sqlQuerry = "SELECT * FROM getDataFromTRBs(" + l_run_id + ");";
  std::string l_run_id = boost::lexical_cast<std::string>(p_run_id);
  pqxx::result l_runDbResults = getDataFromDB("getDataFromTRBs",l_run_id);

  size_t l_sizeResultQuerry = l_runDbResults.size();

  if (l_sizeResultQuerry) {
    for (pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row) {
      JPetTRB l_TRB = generateTRB(row);
      paramBank.addTRB(l_TRB);
    }
  } else {
    printErrorMessageDB("getDataFromTRBs", p_run_id);
  }
}

void JPetDBParamGetter::fillTOMB(const int p_run_id, JPetParamBank& paramBank)
{
  INFO("Start filling TOMBs container.");

//  std::string l_sqlQuerry = "SELECT * FROM getDataFromTOMB(" + l_run_id + ");";
  std::string l_run_id = boost::lexical_cast<std::string>(p_run_id);
  pqxx::result l_runDbResults = getDataFromDB("getDataFromTOMB",l_run_id);

  size_t l_sizeResultQuerry = l_runDbResults.size();

  if (l_sizeResultQuerry) {
    //for(pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row)
    {
      pqxx::result::const_iterator row = l_runDbResults.begin();
      JPetTOMB l_TOMB = generateTOMB(row);
      paramBank.setTOMB(l_TOMB);
    }
  } else {
    printErrorMessageDB("getDataFromTOMB", p_run_id);
  }
}

JPetScin JPetDBParamGetter::generateScintillator(pqxx::result::const_iterator row) {
      int l_scintillator_id = row["scintillator_id"].as<int>();

      double l_scintillator_length = row["scintillator_length"].as<double>();
      double l_scintillator_width = row["scintillator_width"].as<double>();
      double l_scintillator_height = row["scintillator_height"].as<double>();

      int l_setup_id = row["setup_id"].as<int>();
      int l_run_id = row["run_id"].as<int>();
      

      JPetScin l_scin(l_scintillator_id,
                      0.f,			/// @todo what is attenuation length in database?
                      l_scintillator_length,
                      l_scintillator_height,
                      l_scintillator_width);

      return l_scin;
}


JPetPM JPetDBParamGetter::generatePM(pqxx::result::const_iterator row) {
      int l_hvpmconnection_id = row["hvpmconnection_id"].as<int>();
      bool l_hvpmconnection_isrightside = row["hvpmconnection_isrightside"].as<bool>();

      int l_photomultiplier_id = row["photomultiplier_id"].as<int>();

      int l_setup_id = row["setup_id"].as<int>();
      int l_run_id = row["run_id"].as<int>();

      JPetPM::Side l_side = (l_hvpmconnection_isrightside) ? JPetPM::Side::kRight : JPetPM::Side::kLeft;

      JPetPM l_pm;
      l_pm.setID(l_photomultiplier_id);
      l_pm.setSide(l_side);

      return  l_pm;
}

JPetFEB JPetDBParamGetter::generateFEB(pqxx::result::const_iterator row) {
      int l_konradboard_id = row["konradboard_id"].as<int>();
      bool l_konradboard_isactive = row["konradboard_isactive"].as<bool>();
      std::string l_konradboard_status = row["konradboard_status"].as<std::string>();
      std::string l_konradboard_description = row["konradboard_description"].as<std::string>();
      int l_konradboard_version = row["konradboard_version"].as<int>();
      int l_konradboard_creator_id = row["konradboard_creator_id"].as<int>();

      int l_setup_id = row["setup_id"].as<int>();
      int l_run_id = row["run_id"].as<int>();

      JPetFEB l_FEB(l_konradboard_id,
                    l_konradboard_isactive,
                    l_konradboard_status,
                    l_konradboard_description,
                    l_konradboard_version,
                    l_konradboard_creator_id);
    return l_FEB;
}

JPetTRB JPetDBParamGetter::generateTRB(pqxx::result::const_iterator row) {
      int l_TRB_id = row["TRB_id"].as<int>();

      int l_setup_id = row["setup_id"].as<int>();
      int l_run_id = row["run_id"].as<int>();

      JPetTRB l_TRB(l_TRB_id,
                    0,		/// @todo what is type in database
                    0);		/// @todo what is channel in database
  return l_TRB;
}


JPetTOMB JPetDBParamGetter::generateTOMB(pqxx::result::const_iterator row) {
      int l_TOMB_id = row["TOMB_id"].as<int>();
      std::string l_TOMB_description = row["TOMB_description"].as<std::string>();

      int l_setup_id = row["setup_id"].as<int>();
      int l_run_id = row["run_id"].as<int>();


      JPetTOMB l_TOMB(l_TOMB_id, l_TOMB_description);
    return l_TOMB;
}


void JPetDBParamGetter::fillScintillatorsTRefs(const int p_run_id, JPetParamBank& paramBank)
{
  INFO("Start filling Scintillators TRefs.");

  int l_scintillatorsSize = paramBank.getScintillatorsSize();
  int l_PMsSize = paramBank.getPMsSize();

  if (l_scintillatorsSize > 0 && l_PMsSize > 0) {

    for (unsigned int l_scintillator_index = 0u; l_scintillator_index < l_scintillatorsSize; ++l_scintillator_index) {
///wk!!      paramBank.getScintillator(l_scintillator_index).clearTRefPMs();
//      ((JPetScin*)fScintillators[l_scintillator_index])->clearTRefPMs();

//      std::string l_scitillator_id = boost::lexical_cast<std::string>(((JPetScin*)fScintillators[l_scintillator_index])->getID());
  //    std::string l_sqlQuerry = "SELECT * FROM getPhotoMultipliersForScintillator(" + l_scitillator_id + ");";
    std::string scin_id = boost::lexical_cast<std::string>(paramBank.getScintillator(l_scintillator_index).getID());
  std::string l_run_id = boost::lexical_cast<std::string>(p_run_id);
    std::string args = scin_id + "," + l_run_id;
      pqxx::result l_runDbResults = getDataFromDB("getPhotoMultipliersForScintillator", args);

      size_t l_sizeResultQuerry = l_runDbResults.size();

      if (l_sizeResultQuerry) {
        for (pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row) {
          int l_SLSCConnection_id = row["SLSCConnection_id"].as<int>();
          int l_TOMB_id = row["Slot_id"].as<int>();
          int l_HVPMConnection_id = row["HVPMConnection_id"].as<int>();
          int l_PhotoMultiplier_id = row["PhotoMultiplier_id"].as<int>();

          for (unsigned int l_PM_index = 0u; l_PM_index < l_PMsSize; ++l_PM_index) {
            //        int l_PM_id = ((JPetPM*)fPMs[l_PM_index])->getID();
            int l_PM_id = paramBank.getPM(l_PM_index).getID();

            if (l_PM_id == l_PhotoMultiplier_id) {
              //            JPetPM::Side l_PM_side = ((JPetPM*)fPMs[l_PM_index])->getSide();
              JPetPM::Side l_PM_side = paramBank.getPM(l_PM_index).getSide();

              if (l_PM_side == JPetPM::Side::kLeft) {
                //      ((JPetScin*)fScintillators[l_scintillator_index])->setLeftTRefPM( *((JPetPM*)fPMs[l_PM_index]) );
                paramBank.getScintillator(l_scintillator_index).setLeftTRefPM(paramBank.getPM(l_PM_index));
              } else if (l_PM_side == JPetPM::Side::kRight) {
                //         ((JPetScin*)fScintillators[l_scintillator_index])->setRightTRefPM( *((JPetPM*)fPMs[l_PM_index]) );
                paramBank.getScintillator(l_scintillator_index).setRightTRefPM( paramBank.getPM(l_PM_index) );
              }
            }
          }
        }
      }
    }
  } else {
    if (l_scintillatorsSize == 0)
      ERROR("Scintillators container is empty.");
    if (l_PMsSize == 0)
      ERROR("PMs container is empty.");
  }
}

void JPetDBParamGetter::fillPMsTRefs(const int p_run_id, JPetParamBank& paramBank)
{
  INFO("Start filling PMs TRefs.");

  int l_PMsSize = paramBank.getPMsSize();
  int l_FEBsSize = paramBank.getFEBsSize();

  if (l_PMsSize > 0 && l_FEBsSize > 0) {

    for (unsigned int l_PM_index = 0u; l_PM_index < l_PMsSize; ++l_PM_index) {
//      ((JPetPM*)fPMs[l_PM_index])->clearTRefFEBs();
    ///wk!!  paramBank.getPM(l_PM_index).clearTRefKBs();
    std::string pm_id = boost::lexical_cast<std::string>(paramBank.getPM(l_PM_index).getID());
  std::string l_run_id = boost::lexical_cast<std::string>(p_run_id);
    std::string args = pm_id + "," + l_run_id;

      pqxx::result l_runDbResults = getDataFromDB("getKonradBoardsForPhotoMultiplier", args);

      size_t l_sizeResultQuerry = l_runDbResults.size();

      if (l_sizeResultQuerry) {
        for (pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row) {
          int l_PMFEBConnection_id = row["PMKBConnection_id"].as<int>();
          int l_KonradBoardInput_id = row["KonradBoardInput_id"].as<int>();
          int l_KonradBoard_id = row["KonradBoard_id"].as<int>();

          for (unsigned int l_FEB_index = 0u; l_FEB_index < l_FEBsSize; ++l_FEB_index) {
//            int l_FEB_id = ((JPetFEB*)fFEBs[l_FEB_index])->id();
            int l_FEB_id = paramBank.getFEB(l_FEB_index).id();

            if (l_FEB_id == l_KonradBoard_id) {
//              ((JPetPM*)fPMs[l_PM_index])->setTRefFEB( *((JPetFEB*)fFEBs[l_FEB_index]) );
              paramBank.getPM(l_PM_index).setTRefKB(paramBank.getFEB(l_FEB_index) );
            }
          }
        }
      }
    }
  } else {
    if (l_PMsSize == 0)
      ERROR("PMs container is empty.");
    if (l_FEBsSize == 0)
      ERROR("FEBs container is empty.");
  }
}

void JPetDBParamGetter::fillFEBsTRefs(const int p_run_id, JPetParamBank& paramBank)
{
  INFO("Start filling FEBs TRefs.");

  int l_FEBsSize = paramBank.getFEBsSize();
  int l_TRBsSize = paramBank.getTRBsSize();

  if (l_FEBsSize > 0 && l_TRBsSize > 0) {

    for (unsigned int l_FEB_index = 0u; l_FEB_index < l_FEBsSize; ++l_FEB_index) {

//     ((JPetFEB*)fFEBs[l_FEB_index])->clearTRefTRBs();
 ///wk!!!     paramBank.getFEB(l_FEB_index).clearTRefTRBs();

//      std::string l_FEB_id = boost::lexical_cast<std::string>(((JPetFEB*)fFEBs[l_FEB_index])->id());
    std::string feb_id = boost::lexical_cast<std::string>(paramBank.getFEB(l_FEB_index).id());
  std::string l_run_id = boost::lexical_cast<std::string>(p_run_id);
    std::string args = feb_id + "," + l_run_id;
      pqxx::result l_runDbResults = getDataFromDB("getTRBsForKonradBoard",args);

      size_t l_sizeResultQuerry = l_runDbResults.size();

      if (l_sizeResultQuerry) {
        for (pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row) {
          int l_KonradBoardOutput_id = row["KonradBoardOutput_id"].as<int>();
          int l_FEBTRBConnection_id = row["KBTRBConnection_id"].as<int>();
          int l_TRBInput_id = row["TRBInput_id"].as<int>();
          int l_TRB_id = row["TRB_id"].as<int>();

          for (unsigned int l_TRB_index = 0u; l_TRB_index < l_TRBsSize; ++l_TRB_index) {
//            int l_TRBId = ((JPetTRB*)fTRBs[l_TRB_index])->getID();
            int l_TRBId = paramBank.getTRB(l_TRB_index).getID();

            if (l_TRBId == l_TRB_id) {
              // ((JPetFEB*)fFEBs[l_FEB_index])->setTRefTRB( *((JPetTRB*)fTRBs[l_TRB_index]) );
              paramBank.getFEB(l_FEB_index).setTRefTRB( paramBank.getTRB(l_TRB_index));
            }
          }
        }
      }
    }
  } else {
    if (l_FEBsSize == 0)
      ERROR("FEBs container is empty.");
    if (l_TRBsSize == 0)
      ERROR("TRBs container is empty.");
  }
}

void JPetDBParamGetter::fillTRBsTRefs(const int p_run_id, JPetParamBank& paramBank)
{
  INFO("Start filling TRBs TRefs.");

  int l_TRBsSize = paramBank.getTRBsSize();
  int l_TOMBSize = 1;

  if (l_TRBsSize > 0 && l_TOMBSize > 0) {

    for (unsigned int l_TRB_index = 0u; l_TRB_index < l_TRBsSize; ++l_TRB_index) {
//      ((JPetTRB*)fTRBs[l_TRB_index])->clearTRefTOMB();
 ///wk!!!     paramBank.getTRB(l_TRB_index).clearTRefTOMB();

      //std::string l_TRB_id = boost::lexical_cast<std::string>(((JPetTRB*)fTRBs[l_TRB_index])->getID());
    std::string trb_id = boost::lexical_cast<std::string>(paramBank.getTRB(l_TRB_index).getID());
  std::string l_run_id = boost::lexical_cast<std::string>(p_run_id);
    std::string args = trb_id + "," + l_run_id;
      pqxx::result l_runDbResults = getDataFromDB("getTOMBForTRB",args);

      size_t l_sizeResultQuerry = l_runDbResults.size();

      if (l_sizeResultQuerry) {
        for (pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row) {
          int l_TRBOutput_id = row["TRBOutput_id"].as<int>();
          int l_TRBTOMBConnection_id = row["TRBTOMBConnection_id"].as<int>();
          int l_TOMBInput_id = row["TOMBInput_id"].as<int>();
          int l_TOMB_id = row["TOMB_id"].as<int>();

          for (unsigned int l_TOMB_index = 0u; l_TOMB_index < l_TOMBSize; ++l_TOMB_index) {
//            int l_TOMBId = ((JPetTOMB*)fTOMB[l_TOMB_index])->id();
            int l_TOMBId = paramBank.getTOMB().id();

            if (l_TOMBId == l_TOMB_id) {
              // ((JPetTRB*)fTRBs[l_TRB_index])->setTRefTOMB( *((JPetTOMB*)fTOMB[l_TOMB_index]) );
              paramBank.getTRB(l_TRB_index).setTRefTOMB(paramBank.getTOMBAddress());
            }
          }
        }
      }
    }
  } else {
    if (l_TRBsSize == 0)
      ERROR("TRBs container is empty.");
    if (l_TOMBSize == 0)
      ERROR("TOMBs container is empty.");
  }
}

void JPetDBParamGetter::fillTRefs(JPetDBParamGetter::ParamObjectType type) 
{
//  std::string typeText = getParamObjectTypeAsText(type);
//  std::string entryInfo = "Start filling " + typeText + "s TRefs";
//  INFO(entryInfo);
// 
//
//  int l_TRBsSize = paramBank.getTRBsSize();
//  int l_TOMBSize = 1;
//
//  if (l_TRBsSize > 0 && l_TOMBSize > 0) {
//
//    for (unsigned int l_TRB_index = 0u; l_TRB_index < l_TRBsSize; ++l_TRB_index) {
////      ((JPetTRB*)fTRBs[l_TRB_index])->clearTRefTOMB();
//      paramBank.getTRB(l_TRB_index).clearTRefTOMB();
//
//      //std::string l_TRB_id = boost::lexical_cast<std::string>(((JPetTRB*)fTRBs[l_TRB_index])->getID());
//      pqxx::result l_runDbResults = getDataFromDB("getTOMBForTRB",paramBank.getTRB(l_TRB_index).getID());
//
//      size_t l_sizeResultQuerry = l_runDbResults.size();
//
//      if (l_sizeResultQuerry) {
//        for (pqxx::result::const_iterator row = l_runDbResults.begin(); row != l_runDbResults.end(); ++row) {
//          int l_TRBOutput_id = row["TRBOutput_id"].as<int>();
//          int l_TRBTOMBConnection_id = row["TRBTOMBConnection_id"].as<int>();
//          int l_TOMBInput_id = row["TOMBInput_id"].as<int>();
//          int l_TOMB_id = row["TOMB_id"].as<int>();
//
//          for (unsigned int l_TOMB_index = 0u; l_TOMB_index < l_TOMBSize; ++l_TOMB_index) {
////            int l_TOMBId = ((JPetTOMB*)fTOMB[l_TOMB_index])->id();
//            int l_TOMBId = paramBank.getTOMB().id();
//
//            if (l_TOMBId == l_TOMB_id) {
//              // ((JPetTRB*)fTRBs[l_TRB_index])->setTRefTOMB( *((JPetTOMB*)fTOMB[l_TOMB_index]) );
//              paramBank.getTRB(l_TRB_index).setTRefTOMB(paramBank.getTOMBAddress());
//            }
//          }
//        }
//      }
//    }
//  } else {
//    if (l_TRBsSize == 0)
//      ERROR("TRBs container is empty.");
//    if (l_TOMBSize == 0)
//      ERROR("TOMBs container is empty.");
//  }
//
}


void JPetDBParamGetter::fillAllTRefs(const int p_run_id, JPetParamBank& paramBank)
{
  if (paramBank.getScintillatorsSize() > 0
      && paramBank.getPMsSize() > 0
      && paramBank.getFEBsSize() > 0
      && paramBank.getTRBsSize() > 0
     ) {
    fillScintillatorsTRefs(p_run_id,  paramBank);
    fillPMsTRefs(p_run_id, paramBank);
    fillFEBsTRefs(p_run_id, paramBank);
    fillTRBsTRefs(p_run_id, paramBank);
  } else {
    ERROR("Containers are empty.");
  }
}

