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
 *  @file JPetTaskRunner.h
 */

#ifndef JPETTASKRUNNER_H
#define JPETTASKRUNNER_H

#include "../JPetTaskInterface/JPetTaskInterface.h"
#include "../JPetTaskRunnerInterface/JPetTaskRunnerInterface.h"

class JPetTaskRunner: public JPetTaskRunnerInterface
{
public:
  JPetTaskRunner();
  virtual ~JPetTaskRunner();
  virtual void setTask(std::shared_ptr<JPetTaskInterface> task);
  virtual std::shared_ptr<JPetTaskInterface> getTask() const;
  virtual void runTask() = 0;
protected:
  std::shared_ptr<JPetTaskInterface> fTask; /// maybe as unique_ptr ?
private:
  void operator=(const JPetTaskRunner&);
  JPetTaskRunner(const JPetTaskRunner&);
};
#endif /*  !JPETTASKRUNNER_H */
