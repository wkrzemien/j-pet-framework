/**
 *  @copyright Copyright 2018 The J-PET Framework Authors. All rights reserved.
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
 *  @file JPetLoggerInclude.h
 *  @brief Configuration file for the Logger class
 *  Four independent level of logging are defined: INFO, WARNING, ERROR and DEBUG.
 *  Also general macro is defined that allows to create user own level of logging.
 */

/**
  * @brief Configuration file for the Logger class
  * Four independent level of logging are defined: INFO, WARNING, ERROR and DEBUG.
  * Also general macro is defined that allows to create user own level of logging.
  * See
  * http://www.cs.technion.ac.il/users/yechiel/c++-faq/macros-with-multi-stmts.html
  */

#ifndef JPETLOGGER_INCLUDE_H
#define JPETLOGGER_INCLUDE_H

#include "./JPetLogger/JPetLogger.h"

#ifndef __CINT__
#include <boost/log/attributes/scoped_attribute.hpp>
#include <boost/log/utility/manipulators/add_value.hpp>
#endif

#define CUSTOM_LOG(logger, sev, X)                                                                                                                   \
  if (true) {                                                                                                                                        \
    BOOST_LOG_SEV(logger, sev) << boost::log::add_value("Line", __LINE__) << boost::log::add_value("File", __FILE__)                                 \
                               << boost::log::add_value("Function", __func__) << X;                                                                  \
  } else                                                                                                                                             \
  (void)0

#define INFO(X) CUSTOM_LOG(JPetLogger::getInstance().getSeverity(), boost::log::trivial::info, X)
#define WARNING(X) CUSTOM_LOG(JPetLogger::getInstance().getSeverity(), boost::log::trivial::warning, X)
#define ERROR(X) CUSTOM_LOG(JPetLogger::getInstance().getSeverity(), boost::log::trivial::error, X)
#define DEBUG(X) CUSTOM_LOG(JPetLogger::getInstance().getSeverity(), boost::log::trivial::debug, X)
#define LOG(X, sev) CUSTOM_LOG(JPetLogger::getInstance().getSevarity(), sev, X)

#define SET_MINIMAL_LOG_ERROR JPetLogger::getInstance().setLogLevel(boost::log::trivial::error)     // prints only error messages
#define SET_MINIMAL_LOG_WARNING JPetLogger::getInstance().setLogLevel(boost::log::trivial::warning) // prints error + warning messages
#define SET_MINIMAL_LOG_INFO JPetLogger::getInstance().setLogLevel(boost::log::trivial::info)       // prints error + warning + info messages
#define SET_MINIMAL_LOG_DEBUG JPetLogger::getInstance().setLogLevel(boost::log::trivial::debug)     // prints error + warning + info + debug messages

#define SET_MINIMAL_LOG_LEVEL(X) JPetLogger::getInstance().setLogLevel(X)

// for backward compability
#define ENABLE_DEBUG JPetLogger::getInstance().setLogLevel(boost::log::trivial::debug)
#define DISABLE_DEBUG JPetLogger::getInstance().setLogLevel(boost::log::trivial::info)

#endif /* !JPETLOGGER_INCLUDE_H */
