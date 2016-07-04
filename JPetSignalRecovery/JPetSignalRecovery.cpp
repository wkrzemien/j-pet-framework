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
 *  @file JPetSignalRecovery.cpp
 */

#include "JPetSignalRecovery.h"
#include "HelperMathFunctions.h"

std::vector<double> JPetSignalRecovery::recoverFullSignal(const JPetRawSignal& signal)
{
  using namespace boost::numeric::ublas;
  //those parameters must be read from some other place
  matrix<double> A; 
  matrix<double> P; 
  vector<double> m; 
  vector<int> omega;

  int iNumPointsLead = signal.getNumberOfLeadingEdgePoints();
  int iNumPointsTrai = signal.getNumberOfTrailingEdgePoints();
  int iNumPoints = iNumPointsLead + iNumPointsTrai;

  vector<double> time(iNumPoints), yb(iNumPoints);

  for (unsigned int j = 0; j < iNumPointsLead; j++) 
  {
      time[j] = signal.getPoints(JPetSigCh::Leading, JPetRawSignal::ByThrValue).at(j).getValue();
      yb[j] = signal.getPoints(JPetSigCh::Leading, JPetRawSignal::ByThrValue).at(j).getThreshold();
  }
  for (unsigned int j = 0; j < iNumPointsTrai; j++) 
  {
      // sort elements with respect to time (increasing order)
      int j_inner = iNumPoints - 1 - j;
      time[j_inner] = signal.getPoints(JPetSigCh::Trailing, JPetRawSignal::ByThrValue).at(j).getValue();
      yb[j_inner] = signal.getPoints(JPetSigCh::Trailing, JPetRawSignal::ByThrValue).at(j).getThreshold();
  }
  
  if (yb[0] == -0.06) { 
    // necessary condition to procced signal recovery
 
    // find indexes
    vector<int> omega;
    omega = establish_index(time);
    
    // recover signal - prior data A, P, m are needed !!!!
    vector<double> y_hat_JPET = sig_recovery(A, P, m, omega, yb);
  }
}
