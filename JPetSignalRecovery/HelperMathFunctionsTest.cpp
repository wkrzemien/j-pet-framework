#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HelperMathFunctionsTest
#include <boost/test/unit_test.hpp>

#include "HelperMathFunctions.h"

BOOST_AUTO_TEST_SUITE(FirstSuite)
  
BOOST_AUTO_TEST_CASE( recoveryTest1  )
{
  double result;
  double epsilon = 1e-5;   

  // matrices, vectors declarations 
  matrix<double> A, P;      
  vector<double> m;
  vector<int> omega;
  vector<double> y_hat, yb;
  
  // load matrices, vectors
  A = loadMatrix("matrix_A.txt", 300, 45); 
  P = loadMatrix("matrix_P.txt", 45, 45);   
  m = loadVector("vector_m.txt", 300);
  
  // test data
  y_hat = loadVector("yhat_tst1.txt", 300);	
  yb = loadVector("yb_tst1.txt", 8);
  omega = loadVector("omega_tst1.txt", 8);
    
  // translate matlab indexing into C indexing
  for(unsigned int i = 0; i< omega.size();i++)
    omega(i) = omega(i) - 1;
          
  // recover signal
  vector<double> y_hat_JPET = sig_recovery(A, P, m, omega, yb);
 
  // Calculate infinity norm of difference
  result = 1.0 + norm_inf(y_hat-y_hat_JPET);
    
  BOOST_REQUIRE_CLOSE(result, 1.0, epsilon);
}

BOOST_AUTO_TEST_CASE( polynomialFitTest2  )
{
  vector<float> time(4);
  vector<float> volt(4);
    
  int alfa = 2;
  float v0 = -0.05;
    
  time(0) = 1035.0;   
  time(1) = 1542.0;
  time(2) = 2282.0;  
  time(3) = 2900.0;
    
  volt(0) = -0.06;
  volt(1) = -0.20;
  volt(2) = -0.35;  
  volt(3) = -0.50;
    
  float result = polynomialFit(time, volt,  alfa, v0);
  
  float epsilon = 0.1; 
  BOOST_REQUIRE_CLOSE(result, 793.1, epsilon);
  //BOOST_REQUIRE_EQUAL(result, 3);
}

BOOST_AUTO_TEST_SUITE_END()
