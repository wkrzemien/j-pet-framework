/**
 * @file HelperMathFunctions.h
 * @copyright Copyright (c) 2015, J-PET collaboration 
 * @brief Helper mathematical functions used in the Reconstruction
 * module. 
 */
#ifndef _HELPERMATHFUNCTIONS_H_
#define _HELPERMATHFUNCTIONS_H_
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <iostream>
#include <fstream>

using namespace boost::numeric::ublas;

vector<double> loadVector(const char* filename, int numElem)
{
  vector<double> A(numElem);
  double a;
  std::ifstream myFile(filename);  
  
  for (int i=0; i < numElem ; i++) {       
    if (myFile >> a) 
      A(i) = a;      
       
  }
  
  myFile.close();
    
  return A;
}

matrix<double> loadMatrix(const char* filename, int numRows, int numColumns)
{
  matrix<double> A(numRows, numColumns);
  double a;
  std::ifstream myFile(filename);  
  
  for (int i=0; i < numRows ; i++) {
    for (int j=0; j < numColumns ; j++) {
      if (myFile >> a) 
	A(i,j) = a;
      
    }
  }
  
  myFile.close();
    
  return A;
}

void printMatrix(matrix<double> m) 
{

  for (unsigned int i=0; i < m.size1(); i++) {
    for (unsigned int j=0; j < m.size2(); j++) {
      std::cout << m(i, j);
      if(j+1 != m.size2()) 
	std::cout << "\t";
      
    }
    std::cout << std::endl;
  }
}

matrix<double> gjinverse(matrix<double> &m, bool &singular)
{  
   unsigned const int size = m.size1();

   if (size != m.size2() || size == 0)
   {   
     singular = true;     
     matrix<double> A(0,0);     
     return A;          
   }
     
   if (size == 1)
   {         
     matrix<double> A(1, 1);
     if (m(0,0) == 0.0)
     {
       singular = true;       
       return A;            
    }
    singular = false;    
    A(0,0) = 1/m(0,0);
    
    return A;         
  }

  matrix<double> A(size, 2*size);
  matrix_range<matrix<double> > Aleft(A,range(0, size), range(0, size));
  Aleft = m;
  matrix_range<matrix<double> > Aright(A, range(0, size), range(size, 2*size));
  Aright = identity_matrix<double>(size);

  for (unsigned int k = 0; k < size; k++)
  {
    if ( A(k,k) == 0 ) // XXX: test for "small" instead
    {
      // Find a row(l) to swap with row(k)      
      int l = -1;
      
      for (unsigned int i = k+1; i < size; i++)       
      {      
	if ( A(i,k) != 0 )        
	{        
	  l = i;           
	  break;          	  
	}             	
      }

      if ( l < 0 )      
      {
	std::cerr << "Error:" <<  __FUNCTION__ << ":"        
	<< "Input matrix is singular, because cannot find"        
	<< " a row to swap while eliminating zero-diagonal.";        
	singular = true;        
	return Aleft;        	
      }       
      else        
      {      
	matrix_row<matrix<double> > rowk(A, k);        
	matrix_row<matrix<double> > rowl(A, l);       
	rowk.swap(rowl);  	
      }               
    }         
  }

  // Doing partial pivot 
  for (unsigned int k = 0; k < size; k++)
  {
    // normalize the current row   
    for (unsigned int j = k+1; j < 2*size; j++)
      A(k,j) /= A(k,k);
    A(k,k) = 1;

    // normalize other rows
    for (unsigned int i = 0; i < size; i++)
    {
      if ( i != k )  // other rows  // FIX: PROBLEM HERE   
      {
	if ( A(i,k) != 0 )       
	{
	  for (unsigned int j = k+1; j < 2*size; j++)         
	    A(i,j) -= A(k,j) * A(i,k);
	  A(i,k) = 0;
	  
	}          	
      }              
    }      
  }  
  singular = false;
  return Aright;
}


void InvertMatrix(matrix<double> InputMatrix, matrix<double> &InverseMatrix) 
{  
  // Create a duplicate of the input matrix
  int N = InputMatrix.size1();
  matrix<double> A = InputMatrix;
  
  // Create the permutation matrix
  typedef permutation_matrix<std::size_t> pmatrix;
  pmatrix P(N);
  
  // Assign the identity matrix to the inverse
  InverseMatrix.assign(identity_matrix<double>(N));

  //	LU factorization and substitution
  lu_factorize(A,P);
  lu_substitute(A,P,InverseMatrix);
}

vector<int> establish_index(vector<double> t)    
{ 
  int dlug = t.size();
  vector<int> ind(dlug);   
  int T = 50;              
  int stala = 20;
  double suma_x = 0;
  int suma_xp = 0;
  int index = 0;   
  
  
  for (int i = 0; i < dlug; i++)
  {
    suma_x = (t(i)-t(0))/T; 
    suma_xp = (int)round(suma_x); 
    index = suma_xp + stala;
    
    if(index < 300)
      ind(i) = index;
    
  }
  
  return ind;
}

matrix<double> create_matrix_A_omega(matrix<double> R, vector<int> omega_temp) 
{
  int N = R.size2();		
  int M = omega_temp.size();	
  
  matrix<double> R_omega(M,N); 
  
  vector<double> w(N);
  
  for(unsigned int i=0;i< omega_temp.size();i++)  
  {
    w = row(R,omega_temp(i));    
    row(R_omega,i) = w;     
  }   
  
  return R_omega;
}

vector<double> sig_recovery(matrix<double> A, matrix<double> P, vector<double> m, vector<int> omega, vector<double> yb)
{
  int dl  = omega.size();
  int roz = A.size1();         
  int r = A.size2();   
  
  vector<double> y_hat(roz);
  vector<double> m_omega(dl);   
  matrix<double> A_omega;
  vector<double> y;
  matrix<double> A_omega_Tr;
  matrix<double> B_inv(r,r);  
  vector<double> u(r);
  vector<double> x_hat(r);
  
  // 1. calculate A_omega, m_omega
  A_omega = create_matrix_A_omega(A, omega);   
  for(int i=0; i < dl; i++)               
    m_omega(i) = m(omega(i));
  
  // 2. calculate y = yb - m_omega
  y = yb - m_omega;
  
  // 3. calculate x_hat
  A_omega_Tr = trans(A_omega); 
  axpy_prod(A_omega_Tr, A_omega,P,false);
  InvertMatrix(P, B_inv);  
  axpy_prod(A_omega_Tr, y, u, true); 
  axpy_prod(B_inv,u,x_hat,true);
  
  // 4. calculate y_hat
  axpy_prod(A,x_hat, y_hat, true); 
  y_hat = y_hat + m;

  return y_hat;
}  

float polynomialFit(const vector<float>& t, const vector<float>& v_source, int alfa, float v0)
{
  vector<float> v(v_source);
  int K;
  float tSig = -1.0;
  float a, b;
  float meanT = 0.0;
  float meanV = 0.0;
  float sx  = 0.0;
  float sy  = 0.0;
  float sxy = 0.0;

  if (v.size() != t.size() )
    return tSig;

  K = v.size();

  if ((K < 2) || (alfa < 1)) {
    // no regression possible - return
    if (K == 1) {
      tSig = t(0);
    }

    return tSig;
  }

  // prepare the data to fit:
  // 1. change the sign to positive values
  // 2. calculate power 1/alfa
  for (int i = 0; i < K; i++)
    v(i) = pow(-v(i), 1.0 / alfa);

  // evaluate meanT, meanV
  for (int i = 0; i < K; i++) {
    meanT  = meanT + t(i);
    meanV  = meanV + v(i);
  }
  meanT = meanT / K;
  meanV = meanV / K;

  // evaluate sx, sy, sxy
  for (int i = 0; i < K; i++) {
    sx = sx + (t(i) - meanT) * (t(i) - meanT);
    sy = sy + (v(i) - meanV) * (v(i) - meanV);
    sxy = sxy + (t(i) - meanT) * (v(i) - meanV);
  }

  a = sxy / sx;
  b = meanV - a * meanT;
 

  if (fabs(a) < 1e-10)
    tSig = t(0);
  else {
    if (v0 > 0.0)
      v0 = 0.0;

    tSig = (pow(-v0, 1.0 / alfa) - b) / a;
  }
  return tSig;
}

#endif
