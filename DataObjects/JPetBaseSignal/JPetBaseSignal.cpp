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
 *  @file JPetBaseSignal.cpp
 */

#include "./JPetBaseSignal.h"

ClassImp(JPetBaseSignal);

JPetBaseSignal::JPetBaseSignal() :
  TObject(), fPM(0), fBarrelSlot(0) {
}
JPetBaseSignal::JPetBaseSignal(bool isNull):
  fIsNullObject(isNull)
{}

bool JPetBaseSignal::isNullObject() const
{
   return fIsNullObject;
}

JPetBaseSignal& JPetBaseSignal::getDummyResult()
{
   static JPetBaseSignal DummyResult(true);
   return DummyResult;
}

JPetBaseSignal::~JPetBaseSignal() {
}

void JPetBaseSignal::Clear(Option_t *){

  fPM = NULL;
  fBarrelSlot = NULL;

}
