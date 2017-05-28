/**
 *  @copyright Copyright 2017 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  @file JPetOptionValidator.cpp
 */

#include "./JPetOptionValidator.h"
#include "../JPetCommonTools/JPetCommonTools.h"
#include "../JPetLoggerInclude.h"
#include "./JPetAdditionalValidators.h"

using boost::any_cast;

JPetOptionValidator::JPetOptionValidator(){
  fValidatorMap = generateValidationMap();
}

bool JPetOptionValidator::areCorrectOptions(const std::map<std::string, boost::any>& optionsMap, std::vector<std::string>& isOption)
{
  //auto validationMap = generateValidationMap();
  for (auto & checkGroup : fValidatorMap) {
    if (std::find(isOption.begin(), isOption.end(), checkGroup.first ) != isOption.end()) {
      for (auto & checkFunc : checkGroup.second) {
  //       std::cout<<"areCorrectOptions: "<<checkGroup.first<<std::endl;
        if (( !checkFunc(std::make_pair(checkGroup.first, optionsMap.at(checkGroup.first))) )) {
          ERROR("ERROR VALIDATON FOR " + checkGroup.first);
          return false;
        }
      }
    }
  }
  return true;
}

std::map<std::string, std::vector<bool(*)(std::pair <std::string, boost::any>)> > JPetOptionValidator::generateValidationMap()
{
  std::map<std::string, std::vector<bool(*)(std::pair <std::string, boost::any>)> > validationMap;
  validationMap["range_std::vector<int>"].push_back(&isNumberBoundsInRangeValid);
  validationMap["range_std::vector<int>"].push_back(&isRangeOfEventsValid);
  validationMap["type_std::string"].push_back(&isCorrectFileType);
  validationMap["file_std::vector<std::string>"].push_back(&areFilesValid);
  validationMap["runId_int"].push_back(&isRunIdValid);
  validationMap["localDB_std::string"].push_back(&isLocalDBValid);
  validationMap["outputPath_std::string"].push_back(&isOutputDirectoryValid);
  return validationMap;
}

void JPetOptionValidator::addValidatorFunction(const std::string& name, bool(*validatorFunction)(std::pair <std::string, boost::any>) )
{
  fValidatorMap[name].push_back(validatorFunction);
}

bool JPetOptionValidator::isNumberBoundsInRangeValid(std::pair <std::string, boost::any> option)
{
  if ( any_cast<std::vector<int>>(option.second).size() != 2) {
    ERROR("Wrong number of bounds in range.");
    std::cerr << "Wrong number of bounds in range: " << (any_cast<std::vector<int>>(option.second).size()) << std::endl;
    return false;
  }
  return true;
}

bool JPetOptionValidator::isRangeOfEventsValid(std::pair <std::string, boost::any> option)
{
  if ( any_cast<std::vector<int>>(option.second).at(0) > any_cast<std::vector<int>>(option.second).at(1)) {
    ERROR("Wrong number of bounds in range.");
    std::cerr << "Wrong number of bounds in range: " << (any_cast<std::vector<int>>(option.second).size()) << std::endl;
    return false;
  }
  return true;
}

bool JPetOptionValidator::isCorrectFileType(std::pair <std::string, boost::any> option)
{
  std::string type = any_cast<std::string>(option.second);
  if (type == "hld" || type == "root" || type == "scope" || type == "zip") {
    return true;
  } else {
    ERROR("Wrong type of file.");
    std::cerr << "Possible options: hld, zip, root or scope" << std::endl;
    return false;
  }
}

bool JPetOptionValidator::isRunIdValid(std::pair <std::string, boost::any> option)
{
  if ( any_cast<int>(option.second) <= 0) {
    ERROR("Run id must be a number larger than 0.");
    return false;
  }
  return true;
}

bool JPetOptionValidator::isLocalDBValid(std::pair <std::string, boost::any> option)
{
  if ( !JPetCommonTools::ifFileExisting(any_cast<std::string>(option.second)) ) {
    ERROR("File doed not exist.");
    return false;
  }
  return true;
}

bool JPetOptionValidator::areFilesValid(std::pair <std::string, boost::any> option)
{
  std::vector<std::string> fileNames = any_cast<std::vector<std::string>>(option.second);
  for ( unsigned int i = 0; i < fileNames.size(); i++ ) {
    if ( !JPetCommonTools::ifFileExisting(fileNames[i])) {
      ERROR("File : " + fileNames[i] + " does not exist.");
      return false;
    }
  }
  return true;
}

bool JPetOptionValidator::isOutputDirectoryValid(std::pair <std::string, boost::any> option)
{
  if (!JPetCommonTools::isDirectory(any_cast<std::string>(option.second))) {
    ERROR("Output directory does not exist.");
    return false;
  }
  return true;
}
