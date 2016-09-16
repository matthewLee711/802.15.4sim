/* Daniel R. Reynolds
   SMU Mathematics
   2 July 2015 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "matrix.hpp"
using namespace std;


// Example routine to introduce our "Matrix" class
int main(int argc, char* argv[]) {

  cout << endl;

  cout << "\ncreating scalar variables (default to double in this class):\n";
  double a = 5.0;
  cout << "  double a = 5.0;\n";
  cout << "\n a: " << a << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\ncreating row vectors, editing entries, and printing to screen\n";
  Matrix y(1,3);
  y(0,0) = 1.0;
  y(1) = 4.0;
  y(2) = 6.0;
  cout << "  Matrix y(1,3);     // zero-valued matrix with 1 row, 3 columns\n";
  cout << "  y(0,0) = 1.0;      // access entries starting with 0, using () notation\n";
  cout << "  y(1) = 4.0;        // vectors also allow a single index accessor\n";
  cout << "  y(2) = 6.0;\n";
  cout << "  y.Write();         // the 'Write' routine will print to the screen\n";
  cout << "\n y:\n";
  y.Write();
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\ncreating column vectors and printing to file:\n";
  Matrix x(3);
  x(0,0) = 2.0;
  x(1) = 3.0;
  x.Write("x.txt");
  cout << "  Matrix x(3);       // single integer constructor, equivalent to 'Matrix x(3,1)'\n";
  cout << "  x(0,0) = 2.0;      // again, can access (row,column) entry, starting with 0\n";
  cout << "  x(1) = 3.0;        // again, vectors can also access with a single index\n";
  cout << "                     // Note: unset entries are zero-valued\n";
  cout << "  x.Write(\"x.txt\");  // if you supply a string to 'Write' it will write to that file\n";
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\ncreate Matrix from text file:\n";
  Matrix x2 = Read("x.txt");
  cout << "  Matrix x2 = Read(\"x.txt\");\n";
  cout << "\n x2:\n" << x2;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nother ways to create matrices:\n";
  double Adata[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  Matrix A(2,3,Adata);
  cout << "  double Adata[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};   // data array\n";
  cout << "  Matrix A(2,3,Adata);  // new 3x2 matrix, with a given data array \n";
  cout << "  cout << A;      // streaming output for matrices/vectors also works\n";
  cout << "                  // NOTE: column-major ordering of matrix data\n";
  cout << "\n A:\n" << A;
  cout << " Press [enter] to continue\n";  cin.get();

  double Bdata[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  valarray<double> Bva(Bdata, 6);
  Matrix B(3,2,Bva);
  cout << "  double Bdata[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};   // data array\n";
  cout << "  valarray<double> Bva(Bdata, 6);                    // data valarray\n";
  cout << "  Matrix B(3,2,Bva);  // new 3x2 matrix, with a given data array \n";
  cout << "\n B:\n" << B;
  cout << " Press [enter] to continue\n";  cin.get();
  
  vector<double> Cve = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  Matrix C(2,3,Cve);
  cout << "  vector<double> Cve = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};   // data vector\n";
  cout << "  Matrix C(2,3,Cve);  // new 2x3 matrix, with a given data array \n";
  cout << "\n C:\n" << C;
  cout << " Press [enter] to continue\n";  cin.get();
  
  cout << "\nchanging Matrix entries:\n";
  C(1,2) = 12.0;
  cout << "  C(1,2) = 12.0;      // can change existing matrix entries using ()\n";
  cout << "\n C:\n" << C;
  cout << " Press [enter] to continue\n";  cin.get();
  Cve[5] = 2.0;
  cout << "  Cve[5] = 2.0;       // changing the data array directly leaves Matrix unaffected\n";
  cout << "\n C:\n" << C;
  cout << " Press [enter] to continue\n";  cin.get();
  C[2][1] = -2.0;
  cout << "  C[2][1] = -2.0;     // can alternately change entries using [] (C-style, reverse order)\n";
  cout << "\n C:\n" << C;
  cout << " Press [enter] to continue\n";  cin.get();
  C(5) = 6.0;
  cout << "  C(5) = 6.0;         // or we can use single-index with () to specify the overall entry\n";
  cout << "\n C:\n" << C;
  cout << " Press [enter] to continue\n";  cin.get();
  
  cout << "\ncreating matrices from strings:\n";
  Matrix D("1, 2, 3; 4, 5, 6");
  cout << "  Matrix D(\"1, 2, 3; 4, 5, 6\");  // new 2x3 matrix, specified by a Matlab-like string\n";
  cout << "\n D:\n" << D;
  cout << " Press [enter] to continue\n";  cin.get();
  
  cout << "\nmatrix arithmetic operations:\n";
  Matrix tmp0 = 5*D;
  cout << "  Matrix tmp0 = 5*D;       // scalar multiplication\n";
  cout << "\n tmp0:\n" << tmp0;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nvector arithmetic operations, first let's make some vectors:\n";
  Matrix u(1,3);
  u(0) = 1.0;  u(1) = 2.0;  u(2) = 3.0;
  Matrix v("4, 5, 6");
  Matrix w("7.0; 8.0; 9.0");
  cout << "  Matrix u(1,3);\n";
  cout << "  u(0) = 1.0;  u(1) = 2.0;  u(2) = 3.0;\n";
  cout << "  Matrix v(\"4, 5, 6\");\n";
  cout << "  Matrix w(\"7.0; 8.0; 9.0\");\n";
  cout << "\n u:\n" << u;
  cout << "\n v:\n" << v;
  cout << "\n w:\n" << w;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nmatrix addition/subtraction:\n";
  cout << "  u-v:\n" << u-v;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nillegal matrix addition:\n";
  cout << "  u+w:\n";
  cout << u+w;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nmatrix-vector and/or matrix/matrix multiplication:\n";
  cout << "  A*w:\n" << A*w;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nother matrix creation operations:\n";
  cout << "\n Recall A:\n" << A << endl;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix tmp1 = A.T();
  cout << "  Matrix tmp1 = A.T();                   // matrix transpose\n";
  cout << "\n tmp1:\n" << tmp1;
  cout << " Press [enter] to continue\n";  cin.get();
  vector<double> tmp2 = A.Column(1);
  cout << "  vector<double> tmp2 = A.Column(1);     // column accessor (note: vectors print horizontally)\n";
  cout << "\n tmp2:\n" << tmp2;
  cout << " Press [enter] to continue\n";  cin.get();
  vector<double> tmp3 = A[2];
  cout << "  vector<double> tmp3 = A[2];            // simpler column accessor\n";
  cout << "\n tmp3:\n" << tmp3;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix tmp4 = Linspace(0.0, 1.0, 5, 1);
  cout << "  Matrix tmp4 = Linspace(0.0, 1.0, 5, 1);    // linear span constructor (5x1 result)\n";
  cout << "\n tmp4:\n" << tmp4;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix tmp5 = Logspace(-1.0, 1.0, 1, 5);
  cout << "  Matrix tmp5 = Logspace(-1.0, 1.0, 1, 5);   // log10 span constructor (1x5 result)\n";
  cout << "\n tmp5:\n" << tmp5;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix tmp6 = Random(3,2);
  cout << "  Matrix tmp6 = Random(3,2);              // random matrix constructor\n";
  cout << "\n tmp6:\n" << tmp6;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix tmp7 = Eye(4);
  cout << "  Matrix tmp7 = Eye(4);                   // identity matrix constructor\n";
  cout << "\n tmp7:\n" << tmp7;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nwe can even work with submatrices:\n";
  cout << "\n Recall A:\n" << A << endl;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix tmp8 = A(0,-1,0,1);
  cout << "Let's create a new matrix out of a submatrix of A:\n\n";
  cout << "  Matrix tmp8 = A(0,-1,0,1);              // negative index wraps around from end\n";
  cout << "\n tmp8:\n" << tmp8;
  cout << " Press [enter] to continue\n";  cin.get();
  cout << "Let's now modify tmp8:\n\n";
  tmp8(1,0) = 15.0;
  cout << "  tmp8(1,0) = 15.0;\n";
  cout << "\n tmp8:\n" << tmp8;
  cout << " Press [enter] to continue\n";  cin.get();
  A.Copy(tmp8,0,-1,-2,-1);
  cout << "Let's now stuff tmp8 back into a different part of A:\n\n";
  cout << "  A.Copy(tmp8,0,1,-2,-1);                 // negative values may be used here too\n";
  cout << "\n A:\n" << A;
  cout << " Press [enter] to continue\n";  cin.get();
  
  cout << "\nmatrix transformation operations:\n";
  cout << "\n Recall u:\n" << u;
  cout << "\n Recall v:\n" << v;
  cout << "\n Recall w:\n" << w << endl;
  cout << " Press [enter] to continue\n";  cin.get();
  w.Trans();
  cout << "  w.Trans();                    // transpose-in-place\n";
  cout << "\n w:\n" << w;
  cout << " Press [enter] to continue\n";  cin.get();
  u.LinearSum(2.0, 1.0, v);
  cout << "  u.LinearSum(2.0, 1.0, v);     // in-place linear combination (u = 2u + 1v)\n";
  cout << "\n u:\n" << u;
  cout << " Press [enter] to continue\n";  cin.get();
  y.LinearSum(1.0, w, -1.0, v);
  cout << "  y.LinearSum(1.0, w, -1.0, v); // in-place linear combination (y = 1w - 1v)\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Add(u);
  cout << "  y.Add(u);           // updates y with y+u\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Sub(u);
  cout << "  y.Sub(u);           // updates y with y-u\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Add(5.0);
  cout << "  y.Add(5.0);         // adds 5 to all entries of y\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Mul(u);
  cout << "  y.Mul(u);           // in-place product of corresponding entries\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Mul(2.0);
  cout << "  y.Mul(2.0);         // in-place scaling of y by 2\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Div(u);
  cout << "  y.Div(u);           // in-place quotient of corresponding entries\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Power(0.5);
  cout << "  y.Power(0.5);       // in-place exponentiation of an array\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Copy(u);
  cout << "  y.Copy(u);          // copies u into y\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Const(-2.0);
  cout << "  y.Const(-2.0);      // sets all entries to a constant\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Abs();
  cout << "  y.Abs();            // takes the absolute value of each entry\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y += u;
  cout << "  y += u;             // shortcut for Add()\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y -= 2.0;
  cout << "  y -= 2.0;           // shortcut for Sub()\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y *= v;
  cout << "  y *= v;             // shortcut for Mul()\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y /= v;
  cout << "  y /= v;             // shortcut for Div()\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y /= 3.0;
  cout << "  y /= 3.0;           // scalar /= operator\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y ^= 2.0;
  cout << "  y ^= 2.0;           // shortcut for Power()\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y = 3.0;
  cout << "  y = 3.0;            // shortcut for Const()\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nMatrix/vector inquiry functions:\n";
  cout << "\nRecall A:\n" << A;
  cout << "\n and u:\n" << u;
  cout << "\n and v:\n" << v;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nmatrix dimensions:\n";
  cout << "  A.Rows(): " << A.Rows() << endl;
  cout << "  A.Cols(): " << A.Cols() << endl;
  cout << "  A.Size(): " << A.Size() << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nvector 2-norm, matrix Frobenius norm:\n";
  cout << "  Norm(v): " << Norm(v) << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nmatrix/column vector one-norm, row vector inf-norm:\n";
  cout << "  OneNorm(v): " << OneNorm(v) << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nmatrix/column vector inf-norm, row vector one-norm:\n";
  cout << "  InfNorm(v): " << InfNorm(v) << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nequality comparison operator:\n";
  cout << "  u == v: " << (u == v) << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nsmallest entry of v:\n";
  cout << "  v.Min(): " << v.Min() << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nlargest entry of v:\n";
  cout << "  v.Max(): " << v.Max() << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\ndot-product of two vectors/matrices (must be the same shape):\n";
  cout << "  Dot(u,v): " << Dot(u,v) << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nlinear solver (Gaussian elimination with partial pivoting):\n";
  Matrix E = Random(4,4);
  Matrix E2 = E;
  cout << "  Matrix E = Random(4,4);\n";
  cout << "  Matrix E2 = E;\n";
  cout << "\n E:\n" << E;
  cout << "\n E2:\n" << E2;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix b = Random(4,1);
  Matrix b2 = b;
  cout << "  Matrix b = Random(4,1);\n";
  cout << "  Matrix b2 = b;\n";
  cout << "\n b:\n" << b;
  cout << "\n b2:\n" << b2;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix s = Solve(E,b);
  cout << "  Matrix s = Solve(E,b);                     // solver that allocates solution memory\n";
  cout << "\n s:\n" << s;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix s2(4,1);
  cout << "  Matrix s2(4,1);\n";
  cout << "\n s2:\n" << s2;
  cout << " Press [enter] to continue\n";  cin.get();
  Solve(E2,s2,b2);
  cout << "  Solve(E2,s2,b2);                           // solver that uses existing solution memory\n";
  cout << "\n s2:\n" << s2 << endl;

  cout << "\nSimilar vector operations and linear solvers are also implemented for vector<double> objects:\n";
  vector<double> p = Linspace(0.0, 2.0, 5);
  vector<double> q = Logspace(-2.0, 2.0, 5);
  vector<double> r = Random(5);
  cout << "  vector<double> p = Linspace(0.0, 2.0, 5)\n";
  cout << "  vector<double> q = Logspace(-2.0, 2.0, 5)\n";
  cout << "  vector<double> r = Random(5)\n";
  cout << "\n p:\n" << p;
  cout << "\n q:\n" << q;
  cout << "\n r:\n" << r;
  cout << " Press [enter] to continue\n";  cin.get();
  cout << "  Norm(p): " << Norm(p) << endl;
  cout << " Press [enter] to continue\n";  cin.get();
  cout << "  OneNorm(q): " << OneNorm(q) << endl;
  cout << " Press [enter] to continue\n";  cin.get();
  cout << "  InfNorm(r): " << InfNorm(r) << endl;
  cout << " Press [enter] to continue\n";  cin.get();
  cout << "  Dot(p,q): " << Dot(p,q) << endl;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix F = Random(4,4);
  Matrix F2 = F;
  cout << "  Matrix F = Random(4,4);\n";
  cout << "  Matrix F2 = F;\n";
  cout << "\n F:\n" << F;
  cout << "\n F2:\n" << F2;
  cout << " Press [enter] to continue\n";  cin.get();
  vector<double> rhs = Random(4);
  vector<double> rhs2 = rhs;
  cout << "  vector<double> rhs = Random(4);\n";
  cout << "  vector<double> rhs2 = rhs;\n";
  cout << "\n rhs:\n" << rhs;
  cout << "\n rhs2:\n" << rhs2;
  cout << " Press [enter] to continue\n";  cin.get();
  vector<double> sol = Solve(F,rhs);
  cout << "  vector<double> sol = Solve(F,rhs);\n";
  cout << "\n sol:\n" << sol;
  cout << " Press [enter] to continue\n";  cin.get();
  vector<double> sol2(4);
  cout << "  vector<double> sol2(4);\n";
  cout << "\n sol2:\n" << sol2;
  cout << " Press [enter] to continue\n";  cin.get();
  Solve(F2,sol2,rhs2);
  cout << "  Solve(F2,sol2,rhs2);\n";
  cout << "\n sol2:\n" << sol2 << endl;
  
  return 0;
} // end main

