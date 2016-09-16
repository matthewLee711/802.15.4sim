/* Daniel R. Reynolds
   SMU Mathematics
   19 June 2015 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "matrix.hpp"
using namespace std;

// prototypes of other functions
int GramSchmidt(Matrix& X);


// Example routine to test the Mat class
int main(int argc, char* argv[]) {

  // create a row vector of length 5
  Matrix a(1,5);

  // create a row vec with an existing data array
  double dat1[5] = {0.1, 0.2, 0.3, 0.4, 0.5};
  Matrix b(1, 5, dat1);

  // create a column vec with an existing vector
  vector<double> dat2 = {0.1, 0.2, 0.3, 0.4, 0.5};
  Matrix b2(5, 1, dat2);

  // create a row vec with an existing valarray
  valarray<double> dat3(dat1, 5);
  Matrix b3(1, 5, dat3);

  // create a row vector using linspace
  Matrix c = Linspace(1.0, 5.0, 1, 5);

  // create a column vector using the single integer constructor
  Matrix h(7);

  // output vectors above to screen
  printf("writing array of zeros:\n");
  a.Write();
  printf("writing array of 0.1,0.2,0.3,0.4,0.5:\n");
  b.Write();
  printf("writing (column) array of 0.1,0.2,0.3,0.4,0.5:\n");
  b2.Write();
  printf("writing array of 0.1,0.2,0.3,0.4,0.5:\n");
  b3.Write();
  printf("writing array of 1,2,3,4,5:\n");
  c.Write();
  printf("writing a column vector of 7 zeros:\n");
  h.Write();

  // verify that b has size 5
  if (b.Rows()*b.Cols() != 5) 
    printf("error: incorrect matrix size\n");
  if (b.Size() != 5) 
    printf("error: incorrect matrix size\n");

  // edit entries of a in both matrix forms, and write each entry of a to screen
  a(0,0)  = 10.0;
  a(0,1)  = 15.0;
  a[2][0] = 20.0;
  a[3][0] = 25.0;
  a[4][0] = 30.0;
  printf("entries of a, one at a time: should give 10, 15, 20, 25, 30\n");
  for (size_t i=0; i<a.Cols(); i++) 
    printf("  %g\n",a(i));
  
  // write the values to file
  printf("writing this same vector to the file 'a_data':\n");
  a.Write("a_data");

  // Testing Read() constructor
  double tol = 1.0e-15;
  Matrix read_test1a = Random(3,4);
  read_test1a.Write("tmp.txt");
  Matrix read_test1b = Read("tmp.txt");
  Matrix read_test1_error = read_test1a-read_test1b;
  if (InfNorm(read_test1_error) < tol)
    cout << "Read() test 1 passed\n";
  else {
    cout << "Read() test 1 failed, ||error|| = " << InfNorm(read_test1_error) << endl;
    cout << "  read_test1a = \n" << read_test1a;
    cout << "  read_test1b = \n" << read_test1b;
  }

  Matrix read_test2a = Random(12,1);
  read_test2a.Write("tmp.txt");
  Matrix read_test2b = Read("tmp.txt");
  Matrix read_test2_error = read_test2a-read_test2b;
  if (InfNorm(read_test2_error) < tol)
    cout << "Read() test 2 passed\n";
  else {
    cout << "Read() test 2 failed, ||error|| = " << InfNorm(read_test2_error) << endl;
    cout << "  read_test2a = \n" << read_test2a;
    cout << "  read_test2b = \n" << read_test2b;
  }

  Matrix read_test3a = Random(1,7);
  read_test3a.Write("tmp.txt");
  Matrix read_test3b = Read("tmp.txt");
  Matrix read_test3_error = read_test3a-read_test3b;
  if (InfNorm(read_test3_error) < tol)
    cout << "Read() test 3 passed\n";
  else {
    cout << "Read() test 3 failed, ||error|| = " << InfNorm(read_test3_error) << endl;
    cout << "  read_test3a = \n" << read_test3a;
    cout << "  read_test3b = \n" << read_test3b;
  }
    

  
  // Testing copy constructor
  Matrix B = a;
  cout << "Matrix B = a uses copy constructor, should give 10, 15, 20, 25, 30\n";
  B.Write();
  // update one entry of a
  printf("updating the 4th entry of a to be 31:\n");
  a(4) = 31.0;
  cout << "   a = " << a;

  cout << "B should not have changed" << endl;
  B.Write();
  a(4) = 30.0;  // reset to original

  // Testing submatrix copy constructor
  Matrix B2 = a(0,0,1,3);
  cout << "Matrix B2 = a(0,0,1,3) uses submatrix copy constructor" << endl;
  B2.Write();
  // update entries of B2
  B2(0) = 4.0;
  B2(1) = 3.0;
  B2(2) = 2.0;
  // copy B2 back into a using submatrix copy
  a.Copy(B2,0,0,1,3);
  cout << "submatrix copy into a, should have entries 10 4 3 2 30" << endl;
  a.Write();
  a(1) = 15.0;  // reset to original
  a(2) = 20.0;
  a(3) = 25.0;

  // Testing string-based constructor (Matlab-like)
  Matrix C("1, 2, 3.14, -2.0; 5, -7, 1.2345678901234567890123, 0");
  cout << "Matrix C uses string constructor, should give\n";
  cout << "  1   2  3.14                      -2\n";
  cout << "  5  -7  1.2345678901234567890123   0\n";
  cout << " Actually gives:\n";
  C.Write();

  // Test arithmetic operators
  printf("Testing vector add, should give 1.1, 2.2, 3.3, 4.4, 5.5\n");
  b.Add(c);
  b.Write();

  printf("Testing scalar add, should give 2, 3, 4, 5, 6\n");
  c.Add(1.0);
  c.Write();

  printf("Testing vector subtract, should be 8, 12, 16, 20, 24\n");
  a.Sub(c);
  a.Write();

  printf("Testing scalar subtract, should be 0, 1, 2, 3, 4\n");
  c.Sub(2.0);
  c.Write();

  printf("Testing vector const, should all be -1\n");
  b.Const(-1.0);
  b.Write();

  printf("Testing vector copy, should be 0, 1, 2, 3, 4\n");
  a.Copy(c);
  a.Write();

  printf("Testing scalar mul, should be 0, 5, 10, 15, 20\n");
  c.Mul(5.0);
  c.Write();

  printf("Testing vector mul, should be 0, -1, -2, -3, -4\n");
  b.Mul(a);
  b.Write();

  printf("Testing vector div, should be 0, -2.5, -3.3333, -3.75, -4\n");
  Matrix j(c);
  b.Add(-1.0);
  j.Div(b);
  b.Add(1.0);
  j.Write();

  printf("Testing vector +=, should be 0, 4, 8, 12, 16\n");
  b += c;
  b.Write();

  printf("Testing scalar +=, should be 1, 6, 11, 16, 21\n");
  c += 1.0;
  c.Write();
  
  printf("Testing vector -=, should be 1, 2, 3, 4, 5\n");
  c -= b;
  c.Write();

  printf("Testing scalar -=, should be -2, -1, 0, 1, 2\n");
  a -= 2.0;
  a.Write();

  printf("Testing vector *=, should be 0, -4, 0, 12, 32\n");
  a *= b;
  a.Write();

  printf("Testing scalar *=, should be 2, 4, 6, 8, 10\n");
  c *= 2.0;
  c.Write();

  printf("Testing vector /=, should be 0, -1, 0, 1.5, 3.2\n");
  j = a;
  j /= c;
  j.Write();

  printf("Testing scalar /=, should be 1, 2, 3, 4, 5\n");
  j = c;
  j /= 2.0;
  j.Write();

  printf("Testing vector =, should be 2, 4, 6, 8, 10\n");
  b = c;
  b.Write();
  
  printf("Testing scalar =, should be 3, 3, 3, 3, 3\n");
  a = 3.0;
  a.Write();
  
  printf("Testing Norm, should be 14.8324\n");
  printf("  %g\n", Norm(b));
  
  printf("Testing InfNorm, should be 30\n");
  printf("  %g\n", InfNorm(b));
  
  printf("Testing OneNorm, should be 10\n");
  printf("  %g\n", OneNorm(b));
  
  printf("Testing Min, should be 2\n");
  printf("  %g\n", b.Min());
  
  printf("Testing Max, should be 10\n");
  printf("  %g\n", b.Max());
  
  printf("Testing Dot, should be 90\n");
  printf("  %g\n", Dot(a, c));
  
  printf("Testing Logspace, should be 0.01 0.1 1 10 100\n");
  Matrix e = Logspace(-2.0, 2.0, 1, 5);
  e.Write();
 
  ofstream out;
  out.open("e.txt");
  if(out.is_open())
  {
    out << e;
    cout << "Wrote to file e.txt:\n" << e;
  }
  out.close();
  
  printf("Testing Random\n");
  Matrix f = Random(3,3);
  f.Write();
  cout << "Testing Write with a temporary result" << endl;
  (f*f+f).Write();
  cout << "Testing overloading << with a temporary result (should be same as above)" << endl;
  cout << setprecision(16) << (f*f+f);
  cout << "Testing f==f, should be 1\n  " << (f==f) << endl;
  cout << "Testing e==f, should be 0\n  " << (e==f) << endl;

  // create and fill in a 10x5 matrix
  Matrix Y(10,5);
  for (size_t i=0; i<10; i++) {
    Y(i,0) = 1.0*i;
    Y(i,1) = -5.0 + 1.0*i;
    Y(i,2) = 2.0 + 2.0*i;
    Y(i,3) = 20.0 - 1.0*i;
    Y(i,4) = -20.0 + 1.0*i;
  }

  // extract columns from matrix (both ways)
  vector<double> Y0 = Y.Column(0);
  Matrix Y1(10, 1, Y.Column(1));
  Matrix Y2(10, 1, Y[2]);
  vector<double> Y3 = Y[3];
  vector<double> Y4 = Y[4];

  // check the LinearSum routine
  Y4 += Y3;
  printf("Testing column extraction, should be all zeros:\n");
  cout << Y4;

  Y2 = 0.0;
  Matrix Z(10,1);
  cout << "Testing the vector is all zeros, should be 1\n" << (Y2==Z) << endl;

  // check the LinearSum routine
  Matrix d = Linspace(0.0, 4.0, 1, 5);
  printf("Testing LinearSum, should be 0.02 1.2 4 23 204:\n");
  Matrix g(1,5);
  g.LinearSum(1.0, d, 2.0, e);
  g.Write();

  // check the Power routine
  d.Power(2.0);
  printf("Testing Power, should be 0 1 4 9 16:\n");
  d.Write();
  d ^= 0.5;
  printf("Testing Power, should be 0 1 2 3 4:\n");
  d.Write();

  // check the Abs routine
  Y1.Abs();
  printf("Testing Abs, should be the column 5 4 3 2 1 0 1 2 3 4:\n");
  Y1.Write();

  // check the Trans routine
  printf("Testing Trans, should be the row 5 4 3 2 1 0 1 2 3 4:\n");
  Y1.Trans();
  Y1.Write();

  printf("Testing GramSchmidt, should work\n");
  Matrix X = Random(20,3);
  int iret = GramSchmidt(X);
  vector<double> X0 = X.Column(0);
  vector<double> X1 = X.Column(1);
  vector<double> X2 = X.Column(2);
  printf("  GramSchmidt returned %i, dot-products are:\n",iret);
  printf("     <X0,X0> = %g\n", Dot(X0,X0));
  printf("     <X0,X1> = %g\n", Dot(X0,X1));
  printf("     <X0,X2> = %g\n", Dot(X0,X2));
  printf("     <X1,X1> = %g\n", Dot(X1,X1));
  printf("     <X1,X2> = %g\n", Dot(X1,X2));
  printf("     <X2,X2> = %g\n", Dot(X2,X2));
  
  printf("Testing GramSchmidt, should fail\n");
  Matrix V = Random(20,3);
  V[2] = V[1];
  iret = GramSchmidt(V);
  Matrix V0(20,1,V[0]);
  Matrix V1(20,1,V[1]);
  Matrix V2(20,1,V[2]);
  printf("  GramSchmidt returned %i, dot-products are:\n",iret);
  printf("     <V0,V0> = %g\n", Dot(V0,V0));
  printf("     <V0,V1> = %g\n", Dot(V0,V1));
  printf("     <V0,V2> = %g\n", Dot(V0,V2));
  printf("     <V1,V1> = %g\n", Dot(V1,V1));
  printf("     <V1,V2> = %g\n", Dot(V1,V2));
  printf("     <V2,V2> = %g\n", Dot(V2,V2));

  printf("Testing Matrix product, should be: 9 -1 9 -8 11 6\n");
  Matrix A_ = Eye(6);
  A_(0,3) = 2.0;
  A_(1,2) = -1.0;
  A_(2,5) = 1.0;
  A_(3,5) = -2.0;
  A_(4,5) = 1.0;
  Matrix xtrue_ = Linspace(1.0, 6.0, 1, 6).T();
  Matrix b_ = A_*xtrue_;
  b_.Write();

  printf("Testing BackSub with provided solution array:\n");
  Matrix x_(6,1);
  BackSub(A_, x_, b_);
  //x_.Write();
  cout << "  ||x - xtrue|| = " << InfNorm(x_ - xtrue_) << endl;

  printf("Testing FwdSub:\n");
  A_ = Eye(6);
  A_(3,0) = 2.0;
  A_(2,1) = -1.0;
  A_(5,2) = 1.0;
  A_(5,3) = -2.0;
  A_(5,4) = 1.0;
  b_ = A_*xtrue_;
  x_ = 0.0;
  FwdSub(A_, x_, b_);
  cout << "  ||x - xtrue|| = " << InfNorm(x_ - xtrue_) << endl;

  printf("Testing general solver with vector rhs & solution:\n");
  Matrix C_ = 100.0*Eye(9) + Random(9,9);
  Matrix z_ = Logspace(-4.0, 4.0, 9, 1);
  Matrix f_ = C_*z_;
  Matrix g_ = Solve(C_, f_);
  cout << "  ||x - xtrue|| = " << InfNorm(g_ - z_) << endl;

  printf("Testing subarray copy, should be: \n");
  printf("    0.01   0.02   0.04   0.08\n");
  printf("    0.1    0.2    0.4    0.8\n");
  printf("    1      2      4      8\n");
  printf("   10     20     40     80\n");
  printf(" Actually is:\n");
  Matrix z2_ = Logspace(-2.0, 1.0, 4, 1);
  Matrix B_(4,4);
  B_.Copy(z2_,0,3,0,0);
  z2_ *= 2.0;
  B_.Copy(z2_,0,3,1,1);
  z2_ *= 2.0;
  B_.Copy(z2_,0,3,2,2);
  z2_ *= 2.0;
  B_.Copy(z2_,0,3,3,3);
  B_.Write();

  printf("Testing general solver with matrix rhs & solution:\n");
  Matrix E_ = 100.0*Eye(4) + Random(4,4);
  Matrix F_ = E_*B_;
  Matrix X_ = Solve(E_, F_);
  cout << "  ||X - Xtrue|| = " << InfNorm(X_ - B_) << endl;

  printf("Testing matrix inverse:\n");
  Matrix D_ = 10.0*Eye(8) + Random(8,8);
  Matrix DDinv_(D_);
  Matrix Dinv_ = Inverse(D_);
  DDinv_ = D_*Dinv_;
  cout << "  ||I - D*Dinv|| = " << InfNorm(Eye(8) - DDinv_) << endl;
  DDinv_ = Dinv_*D_;
  cout << "  ||I - Dinv*D|| = " << InfNorm(Eye(8) - DDinv_) << endl;

  return 0;
} // end main

