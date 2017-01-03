#include "./JPetRecoImageTools.h"
#include "../../JPetLoggerInclude.h"
#include <vector>
#include <utility>
#include <cmath>
#include <functional>
#include <cassert>
#include <memory>

JPetRecoImageTools::JPetRecoImageTools() {}

JPetRecoImageTools::~JPetRecoImageTools() {}

double JPetRecoImageTools::nearestNeighbour(Matrix2D& emissionMatrix, double a, double b,
    int center, int x, int y, bool sang)
{
  if (sang) {
    y = (int)std::round(a * x + b) + center;
    if (y >= 0 && y < (int)emissionMatrix[0].size())
      return emissionMatrix[x + center][y];
    else return 0;
  } else {
    x = (int)std::round(a * y + b) + center + 1; // not really know why need to +1 to x, but it works
    if (x >= 0 && x < (int)emissionMatrix.size())
      return emissionMatrix[x][y + center];
    else return 0;
  }
}

double JPetRecoImageTools::linear(Matrix2D& emissionMatrix, double a, double b,
                                  int center, int x, int y, bool sang)
{
  if (sang) {
    y = (int)std::round(a * x + b) + center;
    double weight = std::abs((a * x + b) - std::ceil(a * x + b));
    if (y >= 0 && y < (int)emissionMatrix[0].size()) {
      return (1 - weight) * emissionMatrix[x + center][y] + weight * emissionMatrix[x + center][y + 1];
    } else return 0;
  } else {
    x = (int)std::round(a * y + b) + center + 1; // not really know why need to +1 to x, but it works
    double weight = std::abs((a * y + b) - std::ceil(a * y + b));
    if (x >= 0 && x + 1 < (int)emissionMatrix.size()) {
      if (weight == 0)
        return emissionMatrix[x][y + center];
      else
        return (1 - weight) * emissionMatrix[x][y + center] + weight * emissionMatrix[x + 1][y + center];
    } else return 0;
  }
}


JPetRecoImageTools::Matrix2DProj JPetRecoImageTools::sinogram(Matrix2D& emissionMatrix,
    int views, int scans,
    InterpolationFunc interpolationFunction,
    float ang1, float ang2, bool scaleResult, int min, int max)
{
  assert(emissionMatrix.size() > 0);
  assert(emissionMatrix.size() == emissionMatrix[0].size());
  assert(views > 0);
  assert(scans > 0);
  assert(min < max);
  assert(ang1 < ang2);

  //create vector of size views, initialize it with vector of size scans
  Matrix2DProj proj(views, std::vector<double>(scans));

  float phi = 0.;
  float stepsize = (ang2 - ang1) / views;
  assert(stepsize > 0); //maybe != 0 ?

  const int inputMatrixSize = emissionMatrix.size();
  const int center = inputMatrixSize / 2;

  //if no. scans is greater than the image width, then scale will be <1
  const double scale = inputMatrixSize / scans; //* 1.42
  const double sang = std::sqrt(2) / 2;

  int N = 0;
  int i = 0;

  double sinValue = 0., cosValue = 0.;
  double a = 0., aa = 0.;

  for (phi = ang1; phi < ang2; phi = phi + stepsize) {
    sinValue = std::sin((double)phi * M_PI / 180 - M_PI / 2);
    cosValue = std::cos((double)phi * M_PI / 180 - M_PI / 2);
    if (std::abs(sinValue) < 0.0000001) {
      sinValue = 0;
    }
    if (std::abs(cosValue) < 0.0000001) {
      cosValue = 0;
    }
    a = -cosValue / sinValue;
    if (std::isinf(a) || std::isinf(-a) || std::abs(a) < 0.0000001) {
      a = 0;
      aa = 0;
    } else aa = 1 / a;
    for (int scanNumber = 0; scanNumber < scans; scanNumber++) {
      N = scanNumber - scans / 2;
      proj[i][scans - 1 - scanNumber] = JPetRecoImageTools::calculateValue(emissionMatrix, std::abs(sinValue) > sang,
                                        N, cosValue, sinValue, scale, center,
                                        interpolationFunction, a, aa);
    }
    i++;
  }

  if (scaleResult) {
    JPetRecoImageTools::rescale(proj, min, max);
  }

  return proj;
}

double JPetRecoImageTools::calculateValue(Matrix2D& emissionMatrix, bool sang, int N, double cos,
    double sin, double scale, int center,
    InterpolationFunc& interpolationFunction,
    double a, double aa)
{
  double b = 0.;
  if (sang)
    b = (N - cos - sin) / sin;
  else
    b = (N - cos - sin) / cos;
  b *= scale;
  double value = 0.;
  int x = 0;
  int y = 0;
  if (sang) {
    for (x = -center; x < center; x++) {
      value += interpolationFunction(emissionMatrix, a, b, center, x, y, sang);
    }
    value /= std::abs(sin);
  } else {
    for (y = -center; y < center; y++) {
      value += interpolationFunction(emissionMatrix, aa, b, center, x, y, sang);
    }
    value /= std::abs(cos);
  }
  return value;
}

void JPetRecoImageTools::rescale(Matrix2DProj& matrix, double minCutoff, double rescaleFactor)
{

  double datamax = matrix[0][0];
  double datamin = matrix[0][0];
  for (unsigned int k = 0; k < matrix.size(); k++) {
    for (unsigned int j = 0; j < matrix[k].size(); j++) {
      ///Applying min Cutoff
      if (matrix[k][j] < minCutoff)
        matrix[k][j] = minCutoff;
      /// Finding the largest and the smallest element in the matrix.
      if (matrix[k][j] > datamax)
        datamax = matrix[k][j];
      if (matrix[k][j] < datamin)
        datamin = matrix[k][j];
    }
  }

  if (datamax == 0.) {
    WARNING("Maximum value in the matrix to rescale is 0. No rescaling performed.");
    return;
  }

  for (unsigned int k = 0; k < matrix.size(); k++) {
    for (unsigned int j = 0; j < matrix[0].size(); j++) {
      matrix[k][j] = (matrix[k][j] - datamin) * rescaleFactor / datamax;
    }
  }
}
