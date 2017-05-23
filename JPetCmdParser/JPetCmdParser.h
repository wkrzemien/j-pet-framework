/**
 *  @copyright Copyright 2016 The J-PET Framework Authors. All rights reserved.
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
 *  @file JPetCmdParser.h
 */

#ifndef _JPET_CMD_PARSER_H_
#define _JPET_CMD_PARSER_H_

class JPetCmdParser;

#include "boost/program_options.hpp" // Library parsing command line arguments
#include <string>
#include "../JPetOptions/JPetOptions.h"
#include "../JPetOptionsGenerator/JPetOptionsGenerator.h"


namespace po = boost::program_options;

class JPetCmdParser
{
public:
  JPetCmdParser();
  ~JPetCmdParser();
  std::vector<JPetOptions> parseAndGenerateOptions(int argc, const char** argv);

  inline const po::options_description getOptionsDescription() const {
    return fOptionsDescriptions;
  }
  JPetOptionsGenerator& getGenerator(){
    return fGenerator;
  }
protected:
  po::options_description fOptionsDescriptions;

private:
  JPetOptionsGenerator& fGenerator;
  JPetCmdParser(const JPetCmdParser& cmdParser);
  JPetCmdParser& operator=(const JPetCmdParser& cmdParser);
};

#endif /* _JPET_CMD_PARSER_H_ */
