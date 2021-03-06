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
 *  @file JPetProgressBarManager.h
 */

#ifndef JPETPROGRESSBARMANAGER_H
#define JPETPROGRESSBARMANAGER_H

/**
 * @brief Class managing the progress bar used in while processing events.
 *
 */
class JPetProgressBarManager
{
public:
  void display(long long currentEventNumber, long long numberOfEvents) const;
  float getCurrentValue(int currentEventNumber, int numberOfEvents) const;
};
#endif /*  !JPETPROGRESSBARMANAGER_H */
