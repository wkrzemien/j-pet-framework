/**
 *  @copyright Copyright (c) 2015, Wojciech Krzemien
 *  @file JPetTaskIO.h 
 *  @author Wojciech Krzemien, wojciech.krzemien@if.uj.edu.pl
 *  @brief
 */ 

#ifndef JPETTASKIO_H 
#define JPETTASKIO_H 
#include "../JPetTaskInterface/JPetTaskInterface.h"
#include "JPetTaskIO.h"

class JPetWriter;
class JPetReader;
class JPetTreeHeader;
class JPetTask;

/**
 * @brief Class representing computing task with input/output operations.
 *
 */
class JPetTaskIO: public JPetTaskInterface
{
 public:
  JPetTaskIO();
  virtual void init(const JPetTaskInterface::Options& opts);
  virtual void exec();
  virtual void terminate();
  virtual ~JPetTaskIO();
 protected:
  virtual void createInputObjects(const char* inputFilename);
  virtual void createOutputObjects(const char* outputFilename);
  JPetTask* fTask;
  int fEventNb;
  int fEvent;
  JPetWriter* fWriter;
  JPetReader* fReader;
  JPetTreeHeader* fHeader;

};
#endif /*  !JPETTASKIO_H */
