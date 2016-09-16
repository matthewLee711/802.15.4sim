/* Daniel R. Reynolds
   SMU Mathematics
   19 June 2015 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <fstream>
#include "matrix.hpp"
using namespace std;

// singularity tolerance
#define STOL 1.e-15


// This file implements the operations defined in Matrix class, along with a
// variety of linear algebra functions defined based on matrices and vectors

//-------------------------------------------------------------

///// general Matrix class routines /////

// general constructor (initializes values to 0.0)
Matrix::Matrix(size_t m, size_t n) {
  nrows = m;
  ncols = n;
  data.resize(n);
  for (size_t j=0; j<n; j++)
    data[j].resize(m);
}

// column-vector matrix constructor (initializes values to 0.0)
Matrix::Matrix(size_t m) {
  nrows = m;
  ncols = 1;
  data.resize(1);
  data[0].resize(m);
}

// constructor that copies input data (double array)
Matrix::Matrix(size_t m, size_t n, double* vals) {
  nrows = m;
  ncols = n;
  data.resize(n);
  for (size_t j=0; j<n; j++)
    data[j].resize(m);
  size_t idx=0;
  for (size_t j=0; j<n; j++)
    for (size_t i=0; i<m; i++)
      data[j][i] = vals[idx++];
}

// constructor that copies input data (1D vector)
Matrix::Matrix(size_t m, size_t n, std::vector<double> vals) {
  if (m*n != vals.size())
    cerr << "Matrix constructor error: incompatible shape with vector length\n";
  nrows = m;
  ncols = n;
  data.resize(n);
  for (size_t j=0; j<n; j++)
    data[j].resize(m);
  size_t idx=0;
  for (size_t j=0; j<n; j++)
    for (size_t i=0; i<m; i++)
      data[j][i] = vals[idx++];
}

// constructor that copies input data (1D valarray)
Matrix::Matrix(size_t m, size_t n, std::valarray<double> vals) {
  if (m*n != vals.size())
    cerr << "Matrix constructor error: incompatible shape with valarray length\n";
  nrows = m;
  ncols = n;
  data.resize(n);
  for (size_t j=0; j<n; j++)
    data[j].resize(m);
  size_t idx=0;
  for (size_t j=0; j<n; j++)
    for (size_t i=0; i<m; i++)
      data[j][i] = vals[idx++];
}

// constructor that copies input data (1D vector) into a column vector
Matrix::Matrix(std::vector<double> vals) {
  nrows = vals.size();
  ncols = 1;
  data.resize(1);
  data[0].resize(nrows);
  for (size_t i=0; i<nrows; i++)
    data[0][i] = vals[i];
}

// constructor that copies input data (1D valarray) into a column vector
Matrix::Matrix(std::valarray<double> vals) {
  nrows = vals.size();
  ncols = 1;
  data.resize(1);
  for (size_t i=0; i<nrows; i++)
    data[0][i] = vals[i];
}

// constructor that copies input data (2D vector)
Matrix::Matrix(vector< vector<double> > vals) {
  ncols = vals.size();
  nrows = vals[0].size();
  for (size_t j=0; j<ncols; j++)
    if (vals[j].size() != nrows)
      cerr << "Matrix constructor error: rows in 2D vector must have the same length\n";
  for (size_t j=0; j<ncols; j++)
    for (size_t i=0; i<nrows; i++)
      data[j][i] = vals[j][i];
}

// constructor that copies input data (2D valarray)
Matrix::Matrix(valarray< valarray<double> > vals) {
  ncols = vals.size();
  nrows = vals[0].size();
  for (size_t j=0; j<ncols; j++)
    if (vals[j].size() != nrows)
      cerr << "Matrix constructor error: rows in 2D valarray must have the same length\n";
  for (size_t j=0; j<ncols; j++)
    for (size_t i=0; i<nrows; i++)
      data[j][i] = vals[j][i];
}

// string utility routines we'll need for string-based constructor
vector<string>& split(const string& s, char delim, vector<string>& elems) {
  stringstream ss(s);
  string item;
  while (getline(ss, item, delim)) 
    elems.push_back(item);
  return elems;
}
vector<string> split(const string& s, char delim) {
  vector<string> elems;
  split(s, delim, elems);
  return elems;
}

// string-based matrix constructor
Matrix::Matrix(string mat_spec) {

  // initialize empty matrix
  nrows = 0;
  ncols = 0;

  // parse string to determine matrix rows
  vector<string> rows = split(mat_spec, ';');
  size_t m = rows.size();
  if (m < 1) {
    cerr << "string Matrix constructor error: empty string!\n";
    return;
  }

  // verify that all rows have the same number of columns
  vector<string> strrow = split(rows[0], ',');
  size_t n = strrow.size();
  for (size_t i=1; i<m; i++) {
    strrow = split(rows[i], ',');
    if (strrow.size() != n) {
      cerr << "string Matrix constructor error: all rows must have the same number of columns!\n";
      return;
    }
  }
   
  // allocate matrix data
  data.resize(n);
  for (size_t j=0; j<n; j++)
    data[j].resize(m, 0.0);

  // fill Matrix structure
  nrows = m;
  ncols = n;
  for (size_t i=0; i<nrows; i++) {
    strrow = split(rows[i], ',');
    for (size_t j=0; j<ncols; j++) {
      stringstream ss(strrow[j]);
      ss >> data[j][i];
    }
  }
}

// copy constructor
Matrix::Matrix(const Matrix& A) {
  nrows = A.nrows;
  ncols = A.ncols;
  data.resize(ncols);
  for (size_t j=0; j<ncols; j++) {
    data[j].resize(nrows);
    data[j] = A.data[j];
  }
}

// C = A
Matrix& Matrix::operator=(const Matrix& A) {
  nrows = A.nrows;
  ncols = A.ncols;
  data.resize(ncols);
  for (size_t j=0; j<ncols; j++) {
    data[j].resize(nrows);
    data[j] = A.data[j];
  }
  return *this;
}

// dimension accessor routines
size_t Matrix::Size() const { 
  return ncols*nrows; 
}
size_t Matrix::Cols() const { 
  return ncols; 
}
size_t Matrix::Rows() const { 
  return nrows; 
}

// column accessor routines
vector<double>& Matrix::Column(size_t i) {
  return data[i];
}
vector<double>& Matrix::operator[](size_t i) {
  return data[i];
}

// Matlab/Fortran Matrix accessors (entry (i,j))
double& Matrix::operator()(size_t i, size_t j) {
  return data[j][i];
}
double Matrix::operator()(size_t i, size_t j) const {
  return data[j][i];
}
double& Matrix::operator()(size_t idx) {
  return data[idx/nrows][idx%nrows];
}
double Matrix::operator()(size_t idx) const {
  return data[idx/nrows][idx%nrows];
}

// write myself to stdout
int Matrix::Write() const {

  // print data to screen 
  for (size_t i=0; i<nrows; i++) {
    for (size_t j=0; j<ncols; j++)
      printf("  %.17g", data[j][i]);
    printf("\n");
  }
  return 0;
}

// write myself to a file
int Matrix::Write(const char *outfile) const {

  // return with failure if 'outfile' is empty
  if (strlen(outfile) < 1) {
    cerr << "Matrix::Write error, empty outfile\n";
    return 1;
  }

  // open output file
  FILE *fptr = NULL;
  fptr = fopen(outfile, "w");
  if (fptr == NULL) {
    cerr << "Matrix::Write error, unable to open " << outfile << " for writing\n";
    return 1;
  }

  // print data to file
  for (size_t i=0; i<nrows; i++) {
    for (size_t j=0; j<ncols; j++)
      fprintf(fptr, "  %.16g", data[j][i]);
    fprintf(fptr, "\n");
  }

  // close output file and return
  fclose(fptr);
  return 0;
}

// streaming output routine
std::ostream& operator<<(std::ostream& os, const Matrix& A) {
  for(size_t i=0; i<A.Rows(); i++) {
    for(size_t j=0; j<A.Cols(); j++)
      os << "  " << A.data[j][i];
    os << "\n";
  }
  return os;
}


///// Arithmetic operations defined on a given Mat /////

// C = A*X
int Matrix::Matvec(const Matrix& A, const Matrix& X) {

  // check that array sizes match
  if (A.nrows != nrows || A.ncols != X.nrows || X.ncols != ncols) {
    cerr << "Matrix::Axpy error, matrix size mismatch\n";
    cerr << "  Matrix 1 is " << A.nrows << " x " << A.ncols 
	 << "  Matrix 2 is " << X.nrows << " x " << X.ncols 
	 << ", output is " << nrows << " x " << ncols << endl;
    return 1;
  }
  
  // perform operation
  this->Const(0.0);
  for (size_t k=0; k<ncols; k++) 
    for (size_t j=0; j<A.Cols(); j++)
      for (size_t i=0; i<A.Rows(); i++) 
	(*this)(i,k) += A(i,j)*X(j,k);
    
  // return success
  return 0;

}


// C = a*C + b*B
int Matrix::LinearSum(double a, double b, const Matrix& B) {

  // check that array sizes match
  if (B.nrows != nrows || B.ncols != ncols) {
    cerr << "Matrix::LinearSum error, matrix size mismatch\n";
    cerr << "  Matrix 1 is " << nrows << " x " << ncols 
	 << ", Matrix 2 is " << B.nrows << " x " << B.ncols << endl;
    return 1;
  }
  
  // perform operation
  for (size_t j=0; j<B.ncols; j++)
    for (size_t i=0; i<B.nrows; i++)
      data[j][i] = a*data[j][i] + b*B.data[j][i];
    
  // return success
  return 0;
}

// C = a*A + b*B
int Matrix::LinearSum(double a, const Matrix& A, double b, const Matrix& B) {

  // check that array sizes match
  if (A.nrows != nrows || A.ncols != ncols || 
      B.nrows != nrows || B.ncols != ncols) {
    cerr << "Matrix::LinearSum error, matrix size mismatch\n";
    cerr << "  Matrix 1 is " << nrows << " x " << ncols 
	 << ",  Matrix 2 is " << A.nrows << " x " << A.ncols 
	 << ",  Matrix 3 is " << B.nrows << " x " << B.ncols << endl;
    return 1;
  }
  
  // perform operation
  for (size_t j=0; j<B.ncols; j++)
    for (size_t i=0; i<B.nrows; i++)
      data[j][i] = a*A.data[j][i] + b*B.data[j][i];
  
  // return success
  return 0;
}

// C = C+a  (adds scalar a to my data)
int Matrix::Add(double a) {
  for (size_t j=0; j<ncols; j++)
    for (size_t i=0; i<nrows; i++)
      data[j][i] += a;
  return 0;
}

// C = C.*A (component-wise multiply of my data by A)
int Matrix::Mul(const Matrix& A) {

  // check that array sizes match
  if (A.nrows != nrows || A.ncols != ncols) {
    cerr << "Matrix::Mul error, matrix size mismatch\n";
    cerr << "  Matrix 1 is " << nrows << " x " << ncols 
	 << ",  Matrix 2 is " << A.nrows << " x " << A.ncols << endl;
    return 1;
  }
  
  // perform operation
  for (size_t j=0; j<A.ncols; j++)
    for (size_t i=0; i<A.nrows; i++)
      data[j][i] *= A.data[j][i];
  
  // return success
  return 0;
}

// C = a*C  (scales my data by scalar a)
int Matrix::Mul(double a) {

  // perform operation
  for (size_t j=0; j<ncols; j++)
    for (size_t i=0; i<nrows; i++)
      data[j][i] *= a;
  
  // return success
  return 0;
}

// C = C./A (component-wise division of my data by A)
int Matrix::Div(const Matrix& A) {

  // check that array sizes match
  if (A.nrows != nrows || A.ncols != ncols) {
    cerr << "Matrix::Div error, matrix size mismatch\n";
    cerr << "  Matrix 1 is " << nrows << " x " << ncols 
	 << ",  Matrix 2 is " << A.nrows << " x " << A.ncols << endl;
    return 1;
  }
  
  // perform operation where A is nonzero, otherwise return with error
  for (size_t j=0; j<ncols; j++)
    for (size_t i=0; i<nrows; i++) {
      if (A(i,j) == 0.0)
	return 1;
      data[j][i] /= A.data[j][i];
    }

  // return success
  return 0;
}

//   C = A  (copies A into C)
int Matrix::Copy(const Matrix& A) {

  // check that array sizes match
  if (A.nrows != nrows || A.ncols != ncols) {
    cerr << "Matrix::Copy error, matrix size mismatch\n";
    cerr << "  Matrix 1 is " << nrows << " x " << ncols 
	 << ",  Matrix 2 is " << A.nrows << " x " << A.ncols << endl;
    return 1;
  }
  
  // perform operation
  for (size_t j=0; j<A.ncols; j++)
    data[j] = A.data[j];
  
  // return success
  return 0;
}

//   C(is:ie,js:je) = A  (copies A into a submatrix of C)
//     is,ie,js,je negative  =>  offset from end of dimension (-1 == end)
int Matrix::Copy(const Matrix& A, long int is, long int ie, long int js, long int je) {

  // update is,ie,js,je if any are negative
  is = (is < 0) ? is+nrows : is;
  ie = (ie < 0) ? ie+nrows : ie;
  js = (js < 0) ? js+ncols : js;
  je = (je < 0) ? je+ncols : je;

  // check that array sizes match
  if (A.nrows != (ie-is+1) || A.ncols != (je-js+1)) {
    cerr << "Matrix::Copy error, matrix size mismatch\n";
    cerr << "  supplied Matrix is " << A.nrows << " x " << A.ncols 
	 << ", but requested submatrix is " << ie-is+1 << " x " 
	 << je-js+1 << endl;
    return 1;
  }
  // check for valid submatrix
  if (is < 0 || is >= nrows) {
    cerr << "Matrix::Copy error, requested submatrix does not exist\n";
    cerr << "  illegal is = " << is << " (matrix has " << nrows << " rows)\n";
    return 1;
  }
  if (ie < 0 || ie >= nrows) {
    cerr << "Matrix::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  illegal ie = " << ie << " (matrix has " << nrows << " rows)\n";
    return 1;
  }
  if (js < 0 || js >= ncols) {
    cerr << "Matrix::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  illegal js = " << js << " (matrix has " << ncols << " columns)\n";
    return 1;
  }
  if (je < 0 || je >= ncols) {
    cerr << "Matrix::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  illegal je = " << je << " (matrix has " << ncols << " columns)\n";
    return 1;
  }
  if (ie < is || je < js) {
    cerr << "Matrix::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  upper index is below lower index: is = " << is << ", ie = " 
	 << ie << ", js = " << js << ", je = " << je << endl;
    return 1;
  }
  
  // perform operation
  for (size_t j=0; j<A.ncols; j++)  
    for (size_t i=0; i<A.nrows; i++)  
      data[j+js][i+is] = A.data[j][i];
  
  // return success
  return 0;
}

//   C = a  (sets all entries of C to the scalar a)
int Matrix::Const(double a) {
  for (size_t j=0; j<ncols; j++)  
    for (size_t i=0; i<nrows; i++)  
      data[j][i] = a;
  return 0;
}

// C = C.^p (component-wise exponentiation of my data to the power p)
int Matrix::Power(double p) {
  for (size_t j=0; j<ncols; j++)
    for (size_t i=0; i<nrows; i++)
      data[j][i] = pow(data[j][i], p);
  return 0;
}

// Cij = |Cij| (component-wise absolute value of my data)
int Matrix::Abs() {
  for (size_t j=0; j<ncols; j++)
    for (size_t i=0; i<nrows; i++)
      data[j][i] = abs(data[j][i]);
  return 0;
}

// Cij = Cji
int Matrix::Trans() {

  // perform operation in place if matrix is square
  if (nrows == ncols) {
    double tmp;
    for (size_t i=0; i<nrows; i++)
      for (size_t j=0; j<i; j++)
	std::swap( data[j][i], data[i][j] );

  // otherwise we need a new data array to ensure a clean transpose
  } else {

    // create temporary data array
    vector< vector<double> > newdata;
    newdata.resize(nrows);
    for (size_t j=0; j<nrows; j++)
      newdata[j].resize(ncols);

    // copy transposed data over 
    for (size_t j=0; j<ncols; j++)
      for (size_t i=0; i<nrows; i++)
	newdata[i][j] = data[j][i];

    // copy newdata values into existing data array
    data = newdata;

    // swap matrix dimensions
    std::swap(nrows, ncols);
  }
  
  // return success
  return 0;
}

// computes the inverse of a nonsingular matrix 
int Matrix::Inverse() {

  // check that matrix sizes match
  if (nrows != ncols) {
    cerr << "Inverse error, non-square matrix\n";
    cerr << "  Matrix is " << nrows << " x " << ncols << endl;
    return 1;
  }

  // create two temporary matrices for operation
  Matrix B = Eye(nrows);
  Matrix A = *this;

  // call existing Solve routine for computations
  if (Solve(A, *this, B) != 0)
    return 1;

  // return success
  return 0;
}


///// Derived matrix creation operations /////

// C = A^T
Matrix Matrix::T() {
  Matrix C(*this);
  C.Trans();
  return C;
}

// submatrix copy routine (creates a matrix from a portion of an existing Mat)
//     is,ie,js,je negative  =>  offset from end of dimension (-1 == end)
Matrix Matrix::operator()(long int is, long int ie, long int js, long int je) {

  // update is,ie,js,je if any are negative
  is = (is < 0) ? is+nrows : is;
  ie = (ie < 0) ? ie+nrows : ie;
  js = (js < 0) ? js+ncols : js;
  je = (je < 0) ? je+ncols : je;

  // check that requested submatrix exists
  if (is < 0 || is >= nrows) {
    cerr << "Matrix::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  illegal is = " << is << " (matrix has " << nrows << " rows)\n";
  }
  if (ie < 0 || ie >= nrows) {
    cerr << "Matrix::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  illegal ie = " << ie << " (matrix has " << nrows << " rows)\n";
  }
  if (js < 0 || js >= ncols) {
    cerr << "Matrix::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  illegal js = " << js << " (matrix has " << ncols << " columns)\n";
  }
  if (je < 0 || je >= ncols) {
    cerr << "Matrix::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  illegal je = " << je << " (matrix has " << ncols << " columns)\n";
  }
  if (ie < is || je < js) {
    cerr << "Matrix::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  upper index is below lower index: is = " << is << ", ie = " 
	 << ie << ", js = " << js << ", je = " << je << endl;
  }

  // create new matrix of desired size
  Matrix C(ie-is+1, je-js+1);

  // copy requested data
  for (size_t j=js; j<=je; j++) 
    for (size_t i=is; i<=ie; i++) 
      C.data[j-js][i-is] = data[j][i];

  // return object
  return C;
}


///// Scalar output operators on matrices /////

// minimum entry in the matrix
double Matrix::Min() const {
  double mn=data[0][0];
  for (size_t j=0; j<ncols; j++)
    for (size_t i=0; i<nrows; i++)
      mn = std::min( mn, data[j][i] );
  return mn;
}

// maximum entry in the matrix
double Matrix::Max() const {
  double mx=data[0][0];
  for (size_t j=0; j<ncols; j++)
    for (size_t i=0; i<nrows; i++)
      mx = std::max( mx, data[j][i] );
  return mx;
}

// equivalence-checking operator
bool Matrix::operator==(const Matrix& A) const {

  // quick check for compatible sizes
  if (A.nrows != nrows || A.ncols != ncols)
    return false;

  // detailed check on values
  bool equal = true;
  for (size_t j=0; j<ncols; j++)  
    for (size_t i=0; i<nrows; i++)
      equal &= (A.data[j][i] == data[j][i]);
  return equal;
}


//--- supplementary matrix arithmetic routines ---

// column/row vector dot product, something else for matrices?
double Dot(const Matrix& A, const Matrix& B) {
  double sum=0.0;
  if ((A.Cols() != B.Cols()) || (A.Rows() != B.Rows())) {
    cerr << "Dot error, matrix objects must be the same size\n";
  } else {
    for (size_t j=0; j<A.Cols(); j++) 
      for (size_t i=0; i<A.Rows(); i++)
	sum += A(i,j) * B(i,j);
  }
  return sum;
}

// matrix Frobenius norm (column/row vector 2-norm)
double Norm(const Matrix& A) {
  double sum=0.0;
  for (size_t j=0; j<A.Cols(); j++) 
    for (size_t i=0; i<A.Rows(); i++)
      sum += A(i,j)*A(i,j);
  return sqrt(sum);
}

// matrix infinity norm (column vector infinity norm, row vector one norm)
double InfNorm(const Matrix& A) {
  double mx=0.0;
  for (size_t i=0; i<A.Rows(); i++) {
    double sum=0.0;
    for (size_t j=0; j<A.Cols(); j++) 
      sum += std::abs(A(i,j));
    mx = std::max(mx,sum);
  }
  return mx;
}

// matrix one norm (column vector one norm, row vector infinity norm)
double OneNorm(const Matrix& A) {
  double mx=0.0;
  for (size_t j=0; j<A.Cols(); j++) {
    double sum=0.0;
    for (size_t i=0; i<A.Rows(); i++)
      sum += std::abs(A(i,j));
    mx = std::max(mx,sum);
  }
  return mx;
}


//--- new matrix creation routines (C is the output, A and B are the operands) ---

// C = A+B
Matrix operator+(const Matrix& A, const Matrix& B) {

  // check that array sizes match
  if (B.Rows() != A.Rows() || B.Cols() != A.Cols()) {
    cerr << "Matrix operator+ error, matrix size mismatch\n";
    cerr << "  Matrix 1 is " << A.Rows() << " x " << A.Cols()
	 << ",  Matrix 2 is " << B.Rows() << " x " << B.Cols() << endl;
    return Matrix(0,0);
  }

  // create new Mat for output, and do operation
  Matrix C(A.Rows(),A.Cols());
  for (size_t j=0; j<A.Cols(); j++)
    for (size_t i=0; i<A.Rows(); i++)
      C.data[j][i] = A.data[j][i] + B.data[j][i];

  // return result
  return C;
}

// C = A-B
Matrix operator-(const Matrix& A, const Matrix& B) {

  // check that array sizes match
  if (B.Rows() != A.Rows() || B.Cols() != A.Cols()) {
    cerr << "Matrix operator- error, matrix size mismatch\n";
    cerr << "  Matrix 1 is " << A.Rows() << " x " << A.Cols()
	 << ",  Matrix 2 is " << B.Rows() << " x " << B.Cols() << endl;
    return Matrix(0,0);
  }

  // create new Mat for output, and do operation
  Matrix C(A.Rows(),A.Cols());
  for (size_t j=0; j<A.Cols(); j++)
    for (size_t i=0; i<A.Rows(); i++)
      C.data[j][i] = A.data[j][i] - B.data[j][i];

  // return result
  return C;
}

// C = A*B
Matrix operator*(const Matrix& A, const Matrix& B) {

  // determine if either matrix is a scalar
  bool A_scalar = ((A.Rows()==1) && (A.Cols()==1));
  bool B_scalar = ((B.Rows()==1) && (B.Cols()==1));

  // scalar-times-matrix
  if (A_scalar) {

    // create new Matrix for output, and do operation
    Matrix C(B.Rows(),B.Cols());
    for (size_t j=0; j<B.Cols(); j++)
      for (size_t i=0; i<B.Rows(); i++)
	C.data[j][i] = A.data[0][0] * B.data[j][i];
    return C;
    
  // scalar-times-matrix
  } else if (B_scalar) {
    
    // create new Mat for output, and do operation
    Matrix C(A.Rows(),B.Cols());
    for (size_t j=0; j<A.Cols(); j++)
      for (size_t i=0; i<A.Rows(); i++)
	C.data[j][i] = A.data[j][i] * B.data[0][0];
    return C;
    
  // normal matrix product
  } else {

    // check that array sizes are acceptable
    if (B.Rows() != A.Cols()) {
      cerr << "Matrix operator* error, inner dimension mismatch\n";
      cerr << "  Matrix 1 is " << A.Rows() << " x " << A.Cols()
	   << ",  Matrix 2 is " << B.Rows() << " x " << B.Cols() << endl;
      return Matrix(0,0);
    } else {
    
      // create new Mat for output, and do operation
      Matrix C(A.Rows(),B.Cols());
      for (size_t i=0; i<A.Rows(); i++) 
	for (size_t j=0; j<B.Cols(); j++)
	  for (size_t k=0; k<A.Cols(); k++)
	    C.data[j][i] += A.data[k][i] * B.data[j][k];
      return C;
    }
  }
}


// C = A*b
Matrix operator*(const Matrix& A, const double b) {
  Matrix C(A.Rows(),A.Cols());
  for (size_t j=0; j<A.Cols(); j++)
    for (size_t i=0; i<A.Rows(); i++)
      C.data[j][i] = A.data[j][i] * b;
  return C;
}

// C = a*B
Matrix operator*(const double a, const Matrix& B) {
  Matrix C(B.Rows(),B.Cols());
  for (size_t j=0; j<B.Cols(); j++)
    for (size_t i=0; i<B.Rows(); i++)
      C.data[j][i] = a * B.data[j][i];
  return C;
}

// create a new matrix of linearly spaced data
Matrix Linspace(double a, double b, size_t m, size_t n) {
  Matrix C(m,n);
  double h = (b-a)/(m*n-1);
  size_t idx=0;
  for (size_t j=0; j<n; j++)
    for (size_t i=0; i<m; i++)
      C.data[j][i] = a + (idx++)*h;
  return C;
}

// create a new column-vector matrix of logarithmically spaced data
Matrix Logspace(double a, double b, size_t m, size_t n) {
  Matrix C(m,n);
  double h = (b-a)/(m*n-1);
  size_t idx=0;
  for (size_t j=0; j<n; j++)
    for (size_t i=0; i<m; i++)
      C.data[j][i] = pow(10.0, a + (idx++)*h);
  return C;
}

// create a matrix with uniformly-distributed random numbers in [0,1]
Matrix Random(size_t m, size_t n) {
  Matrix C(m,n);
  for (size_t j=0; j<n; j++)
    for (size_t i=0; i<m; i++)
      C.data[j][i] = random() / (pow(2.0,31.0) - 1.0);
  return C;
}

// create a new n by n identity matrix
Matrix Eye(size_t n) {
  Matrix I(n,n);
  for (size_t i=0; i<n; i++)  I.data[i][i] = 1.0;
  return I;
}

// creates a matrix from a specified input file
Matrix Read(const char *infile) {
  
  // determine matrix size
  size_t nrows=0, ncols=0;
  ifstream ifs;
  string line;
  ifs.open(infile);
  while (getline(ifs, line)) {
    istringstream iss(line);   // convert line to stringstream
    float value;               // determine the number of columns on this row
    size_t n=0;
    while (iss >> value)  n++;
    if ((n > 0) && (nrows == 0))  // first row, set ncols
      ncols = n;
    if ((n > 0) && (n != ncols)) {  // later row, with bad number of columns
      cerr << "Read() error, not all rows in file " << infile 
	   << " have the same number of cols, "
	   << n << " != " << ncols << endl;
      return Matrix(0,0);
    }
    if (n > 0) nrows++;          // legal row, increment counter
  }
  ifs.close();

  // create matrix of desired size
  Matrix A(nrows,ncols);

  // load matrix based on data from file
  ifs.open(infile);   // reopen input file
  for (size_t i=0; i<nrows; i++) {
    getline( ifs, line ); 
    istringstream iss(line);   // convert line to stringstream
    for (size_t j=0; j<ncols; j++) 
      iss >> A.data[j][i];
  }
  ifs.close();

  // return result
  return A;
}

// computes the inverse of a nonsingular matrix A
Matrix Inverse(const Matrix& A) {
  
  // copy A into a new output matrix
  Matrix X(A);

  // call existing Inverse routine for computations
  if (X.Inverse() != 0)
    cerr << "Inverse: error inverting matrix\n";

  // return result
  return X;
}


//--- supplementary matrix-vector arithmetic routines ---

// standard matrix-vector product
vector<double> MatVec(const Matrix& A, const vector<double>& v) {
  vector<double> res(0.0, A.Rows());
  if (A.Cols() != v.size()) {
    cerr << "MatVec: incompatible matrix/vector sizes in A*v\n";
  } else {
    for (size_t i=0; i<A.Rows(); i++) 
      for (size_t j=0; j<A.Cols(); j++)
	res[i] += A(i,j)*v[j];
  }
  return res;
}

vector<double> operator*(const Matrix& A, const std::vector<double>& v) {
  return MatVec(A,v);
}


//--- supplementary vector<double> vector-arithmetic routines ---

// inner product between two vectors
double Dot(const vector<double>& v1, const vector<double>& v2) {
  if (v1.size() != v2.size()) {
    cerr << "Dot: incompatible vector sizes\n";
    return 0.0;
  }
  double res = 0.0;
  for (size_t i=0; i<v2.size(); i++)  res += v1[i]*v2[i];
  return res;
}

// vector 2-norm
double Norm(const vector<double>& v) {
  double sum=0.0;
  for (size_t i=0; i<v.size(); i++)  sum += v[i]*v[i];
  return sqrt(sum);
}

// vector infinity norm
double InfNorm(const vector<double>& v) {
  double mx=0.0;
  for (size_t i=0; i<v.size(); i++)  
    mx = std::max(mx,std::abs(v[i]));
  return mx;
}

// vector one-norm
double OneNorm(const vector<double>& v) {
  double sum=0.0;
  for (size_t i=0; i<v.size(); i++)  sum += std::abs(v[i]);
  return sum;
}

// create a new vector of linearly spaced data
vector<double> Linspace(double a, double b, size_t n) {
  if (n<2) cerr << "Linspace::length must be > 1\n";
  vector<double> v(n);
  double h = (b-a)/(n-1);
  for (size_t i=0; i<n; i++)
    v[i] = a + i*h;
  return v;
}

// create a new vector of logarithmically spaced data
vector<double> Logspace(double a, double b, size_t n) {
  if (n<2) cerr << "Logspace::length must be > 1\n";
  vector<double> v(n);
  double h = (b-a)/(n-1);
  for (size_t i=0; i<n; i++)
    v[i] = pow(10.0, a + i*h);
  return v;
}

// create a new vector with uniformly-distributed random numbers in [0,1]
vector<double> Random(size_t n) {
  if (n<1) cerr << "Random::length must be > 0\n";
  vector<double> v(n);
  for (size_t i=0; i<n; i++)
    v[i] = random() / (pow(2.0,31.0) - 1.0);
  return v;
}

// streaming output routine
std::ostream& operator<<(std::ostream& os, const std::vector<double>& v)
{
  for (size_t i=0; i<v.size(); i++)
    os << "  " << v[i];
  os << "\n";
  return os;
}

// arithmetic operators
vector<double>& operator+=(vector<double>& v, const double c) {
  for (size_t i=0; i<v.size(); i++)
    v[i] += c;
  return v;
}
vector<double>& operator+=(vector<double>& v, const vector<double>& w) {
  if (v.size() != w.size())
    cerr << "vector<double> += error: incompatible vector sizes!";
  else
    for (size_t i=0; i<v.size(); i++)
      v[i] += w[i];
  return v;
}
vector<double>& operator-=(vector<double>& v, const double c) {
  for (size_t i=0; i<v.size(); i++)
    v[i] -= c;
  return v;
}
vector<double>& operator-=(vector<double>& v, const vector<double>& w) {
  if (v.size() != w.size())
    cerr << "vector<double> -= error: incompatible vector sizes!";
  else
    for (size_t i=0; i<v.size(); i++)
      v[i] -= w[i];
  return v;
}
vector<double>& operator*=(vector<double>& v, const double c) {
  for (size_t i=0; i<v.size(); i++)
    v[i] *= c;
  return v;
}
vector<double>& operator*=(vector<double>& v, const vector<double>& w) {
  if (v.size() != w.size())
    cerr << "vector<double> *= error: incompatible vector sizes!";
  else
    for (size_t i=0; i<v.size(); i++)
      v[i] *= w[i];
  return v;
}
vector<double>& operator/=(vector<double>& v, const double c) {
  for (size_t i=0; i<v.size(); i++)
    v[i] /= c;
  return v;
}
vector<double>& operator/=(vector<double>& v, const vector<double>& w) {
  if (v.size() != w.size())
    cerr << "vector<double> /= error: incompatible vector sizes!";
  else
    for (size_t i=0; i<v.size(); i++)
      v[i] /= w[i];
  return v;
}
vector<double>& operator^=(vector<double>& v, const double c) {
  for (size_t i=0; i<v.size(); i++)
    v[i] = pow(v[i], c);
  return v;
}
vector<double>& operator^=(vector<double>& v, const vector<double>& w) {
  if (v.size() != w.size())
    cerr << "vector<double> /= error: incompatible vector sizes!";
  else
    for (size_t i=0; i<v.size(); i++)
      v[i] = pow(v[i], w[i]);
  return v;
}
vector<double> operator+(const vector<double>& v, const double c) {
  vector<double> x(v.size());
  for (size_t i=0; i<v.size(); i++)
    x[i] = v[i] + c;
  return x;
}
vector<double> operator+(const double c, const vector<double>& v) {
  vector<double> x(v.size());
  for (size_t i=0; i<v.size(); i++)
    x[i] = v[i] + c;
  return x;
}
vector<double> operator+(const vector<double>& v, const vector<double>& w) {
  if (v.size() != w.size()) {
    cerr << "vector<double> + error: incompatible vector sizes!";
    return vector<double>(0);
  }
  vector<double> x(v.size());
  for (size_t i=0; i<v.size(); i++)
    x[i] = v[i] + w[i];
  return x;
}
vector<double> operator-(const vector<double>& v, const double c) {
  vector<double> x(v.size());
  for (size_t i=0; i<v.size(); i++)
    x[i] = v[i] - c;
  return x;
}
vector<double> operator-(const double c, const vector<double>& v) {
  vector<double> x(v.size());
  for (size_t i=0; i<v.size(); i++)
    x[i] = c - v[i];
  return x;
}
vector<double> operator-(const vector<double>& v, const vector<double>& w) {
  if (v.size() != w.size()) {
    cerr << "vector<double> - error: incompatible vector sizes!";
    return vector<double>(0);
  }
  vector<double> x(v.size());
  for (size_t i=0; i<v.size(); i++)
    x[i] = v[i] - w[i];
  return x;
}
vector<double> operator*(const vector<double>& v, const double c) {
  vector<double> x(v.size());
  for (size_t i=0; i<v.size(); i++)
    x[i] = v[i] * c;
  return x;
}
vector<double> operator*(const double c, const vector<double>& v) {
  vector<double> x(v.size());
  for (size_t i=0; i<v.size(); i++)
    x[i] = v[i] * c;
  return x;
}
vector<double> operator*(const vector<double>& v, const vector<double>& w) {
  if (v.size() != w.size()) {
    cerr << "vector<double> * error: incompatible vector sizes!";
    return vector<double>(0);
  }
  vector<double> x(v.size());
  for (size_t i=0; i<v.size(); i++)
    x[i] = v[i] * w[i];
  return x;
}
vector<double> operator/(const vector<double>& v, const double c) {
  vector<double> x(v.size());
  for (size_t i=0; i<v.size(); i++)
    x[i] = v[i] / c;
  return x;
}
vector<double> operator/(const double c, const vector<double>& v) {
  vector<double> x(v.size());
  for (size_t i=0; i<v.size(); i++)
    x[i] = c / v[i];
  return x;
}
vector<double> operator/(const vector<double>& v, const vector<double>& w) {
  if (v.size() != w.size()) {
    cerr << "vector<double> / error: incompatible vector sizes!";
    return vector<double>(0);
  }
  vector<double> x(v.size());
  for (size_t i=0; i<v.size(); i++)
    x[i] = v[i] / w[i];
  return x;
}
vector<double> operator^(const vector<double>& v, const double c) {
  vector<double> x(v.size());
  for (size_t i=0; i<v.size(); i++)
    x[i] = pow(v[i], c);
  return x;
}
vector<double> operator^(const double c, const vector<double>& v) {
  vector<double> x(v.size());
  for (size_t i=0; i<v.size(); i++)
    x[i] = pow(c, v[i]);
  return x;
}
vector<double> operator^(const vector<double>& v, const vector<double>& w) {
  if (v.size() != w.size()) {
    cerr << "vector<double> ^ error: incompatible vector sizes!";
    return vector<double>(0);
  }
  vector<double> x(v.size());
  for (size_t i=0; i<v.size(); i++)
    x[i] = pow(v[i], w[i]);
  return x;
}


//--- linear algebra routines ---

// backward substitution on the linear system U*X = B, filling in an existing Matrix X
int BackSub(const Matrix& U, Matrix& X, const Matrix& B) {

  // check that matrix sizes match
  if (U.Rows() != B.Rows() || U.Rows() != U.Cols() || 
      B.Cols() != X.Cols() || X.Rows() != U.Rows()) {
    cerr << "BackSub error, incompatible matrix/vector dimensions\n";
    cerr << "  Matrix is " << U.Rows() << " x " << U.Cols() 
	 << ",  rhs is " << B.Rows() << " x " << B.Cols()
	 << ",  solution is " << X.Rows() << " x " << X.Cols() << endl;
    return 1;
  }
  
  // copy B into X
  X.Copy(B);

  // analyze matrix for magnitude
  double Umax = InfNorm(U);

  // perform column-oriented Backward Substitution algorithm
  for (long int j=U.Rows()-1; j>=0; j--) {

    // check for nonzero matrix diagonal
    if (fabs(U.data[j][j]) < STOL*Umax) {
      cerr << "BackSub error: numerically singular matrix!\n";
      return 1;
    }

    // solve for this row of solution
    for (long int k=0; k<X.Cols(); k++) 
      X.data[k][j] /= U.data[j][j];

    // update all remaining rhs
    for (long int k=0; k<X.Cols(); k++)
      for (long int i=0; i<j; i++)
	X.data[k][i] -= U.data[j][i]*X.data[k][j];

  }

  // return success
  return 0;
}

// backward substitution on the linear system U*X = B, returning X as a new Matrix; 
//    U and B remain unchanged in this operation
Matrix BackSub(const Matrix& U, const Matrix& B) {

  // check that matrix sizes match
  if (U.Rows() != B.Rows() || U.Rows() != U.Cols()) {
    cerr << "BackSub error, illegal matrix/vector dimensions\n";
    cerr << "  Matrix is " << U.Rows() << " x " << U.Cols() 
	 << ",  rhs is " << B.Rows() << " x " << B.Cols() << endl;
    return Matrix(0,0);
  } else {
    // create new Matrix for output and call existing BackSub routine
    Matrix X(U.Rows(),B.Cols());
    if (BackSub(U, X, B) != 0)
      cerr << "BackSub Warning: error in BackSub call\n";
    return X;
  }
}

// backward substitution on U*x = b, filling in an existing vector<double> x
int BackSub(const Matrix& U, vector<double>& x, const vector<double>& b) {

  // check that matrix sizes match
  if (U.Rows() != b.size() || U.Rows() != U.Cols() || x.size() != U.Rows()) {
    cerr << "BackSub error, incompatible matrix/vector dimensions\n";
    cerr << "  Matrix is " << U.Rows() << " x " << U.Cols() 
	 << ",  rhs is " << b.size() << " x 1"
	 << ",  solution is " << x.size() << " x 1\n";
    return 1;
  }
  
  // copy b into x
  x = b;

  // analyze matrix for magnitude
  double Umax = InfNorm(U);

  // perform column-oriented Backward Substitution algorithm
  for (long int j=U.Rows()-1; j>=0; j--) {

    // check for nonzero matrix diagonal
    if (fabs(U.data[j][j]) < STOL*Umax) {
      cerr << "BackSub error: numerically singular matrix!\n";
      return 1;
    }

    // solve for this row of solution
    x[j] /= U.data[j][j];

    // update all remaining rhs
    for (long int i=0; i<j; i++)
      x[i] -= U.data[j][i]*x[j];

  }

  // return success
  return 0;
}

// backward substitution on U*x = b, returning x as a new vector<double>
//    U and b remain unchanged in this operation
vector<double> BackSub(const Matrix& U, const vector<double>& b) {
  if (U.Rows() != b.size() || U.Rows() != U.Cols()) {
    cerr << "BackSub error, illegal matrix/vector dimensions\n";
    cerr << "  Matrix is " << U.Rows() << " x " << U.Cols() 
	 << ",  rhs is " << b.size() << " x 1\n";
    return vector<double>(0);
  }

  // create output vector, call existing BackSub routine, and return
  vector<double> x(0.0, U.Cols());
  if (BackSub(U, x, b) != 0)
    cerr << "BackSub Warning: error in BackSub call\n";
  return x;
}

// forward substitution on the linear system L*X = B, filling in the input Matrix X
//    L and B remain unchanged in this operation; X holds the result
//    B and X may have multiple columns
int FwdSub(const Matrix& L, Matrix& X, const Matrix& B) {

  // check that matrix sizes match
  if (L.Rows() != B.Rows() || L.Rows() != L.Cols() || 
      B.Cols() != X.Cols() || X.Rows() != L.Rows()) {
    cerr << "FwdSub error, illegal matrix/vector dimensions\n";
    cerr << "  Matrix is " << L.Rows() << " x " << L.Cols() 
	 << ",  rhs is " << B.Rows() << " x " << B.Cols()
	 << ",  solution is " << X.Rows() << " x " << X.Cols() << endl;
    return 1;
  }
  
  // copy B into X
  X.Copy(B);

  // analyze matrix magnitude
  double Lmax = InfNorm(L);

  // perform column-oriented Forwards Substitution algorithm
  for (long int j=0; j<L.Rows(); j++) {

    // check for nonzero matrix diagonal
    if (fabs(L.data[j][j]) < STOL*Lmax) {
      cerr << "FwdSub error: singular matrix!\n";
      return 1;
    }

    // solve for this row of solution
    for (long int k=0; k<X.Cols(); k++)
      X.data[k][j] /= L.data[j][j];

    // update all remaining rhs
    for (long int k=0; k<X.Cols(); k++)
      for (long int i=j+1; i<L.Rows(); i++)
	X.data[k][i] -= L.data[j][i]*X.data[k][j];

  }

  // return success
  return 0;
}

// forward substitution on the linear system L*X = B, returning X as a new Matrix
//    L and B remain unchanged in this operation
Matrix FwdSub(const Matrix& L, const Matrix& B) {

  // check that matrix sizes match
  if (L.Rows() != B.Rows() || L.Rows() != L.Cols()) {
    cerr << "FwdSub error, illegal matrix/vector dimensions\n";
    cerr << "  Matrix is " << L.Rows() << " x " << L.Cols() 
	 << ",  rhs is " << B.Rows() << " x " << B.Cols() << endl;
    return Matrix(0,0);
  } else {
    // create new Matrix for output and call existing BackSub routine
    Matrix X(L.Rows(),B.Cols());
    if (FwdSub(L, X, B) != 0)
      cerr << "FwdSub Warning: error in FwdSub call\n";
    return X;
  }
}

// forward substitution on L*x = b, filling in an existing vector<double) x
//    L and b remain unchanged in this operation; x holds the result
int FwdSub(const Matrix& L, vector<double>& x, const vector<double>& b) {

  // check that matrix sizes match
  if (L.Rows() != b.size() || L.Rows() != L.Cols() || 
      x.size() != L.Rows()) {
    cerr << "FwdSub error, illegal matrix/vector dimensions\n";
    cerr << "  Matrix is " << L.Rows() << " x " << L.Cols() 
	 << ",  rhs is " << b.size() << " x 1"
	 << ",  solution is " << x.size() << " x 1\n";
    return 1;
  }
  
  // copy B into X
  x = b;

  // analyze matrix for magnitude
  double Lmax = InfNorm(L);

  // perform column-oriented Forwards Substitution algorithm
  for (long int j=0; j<L.Rows(); j++) {

    // check for nonzero matrix diagonal
    if (fabs(L.data[j][j]) < STOL*Lmax) {
      cerr << "FwdSub error: singular matrix!\n";
      return 1;
    }

    // solve for this row of solution
    x[j] /= L.data[j][j];

    // update all remaining rhs
    for (long int i=j+1; i<L.Rows(); i++)
      x[i] -= L.data[j][i]*x[j];

  }

  // return success
  return 0;
}

// forward substitution on L*x = b, returning x as a new vector<double>
//    L and b remain unchanged in this operation
vector<double> FwdSub(const Matrix& L, const vector<double>& b) {
  // check that matrix sizes match
  if (L.Rows() != b.size() || L.Rows() != L.Cols()) {
    cerr << "FwdSub error, illegal matrix/vector dimensions\n";
    cerr << "  Matrix is " << L.Rows() << " x " << L.Cols() 
	 << ",  rhs is " << b.size() << " x 1\n";
    return vector<double>(0);
  }

  // create output vector and return
  vector<double> x(0.0, L.Cols());
  if (FwdSub(L, x, b) != 0)
    cerr << "FwdSub Warning: error in FwdSub call\n";
  return x;
}

// solves a linear system A*X = B, filling in the input Mat X
int Solve(Matrix& A, Matrix& X, Matrix& B) {

  // check that matrix sizes match
  if (A.Rows() != B.Rows() || A.Rows() != A.Cols() ||
      B.Cols() != X.Cols() || X.Rows() != A.Rows()) {
    cerr << "Solve error, illegal matrix/vector dimensions\n";
    cerr << "  Matrix is " << A.Rows() << " x " << A.Cols() 
	 << ",  rhs is " << B.Rows() << " x " << B.Cols()
	 << ",  solution is " << X.Rows() << " x " << X.Cols() << endl;
    return 1;
  }
  
  // create temporary variables
  long int i, j, k, p;
  double Amax;

  // determine magnitude of entries in A (for singularity check later)
  Amax = InfNorm(A);

  // perform Gaussian elimination to convert A,B to an upper-triangular system
  for (k=0; k<A.Rows()-1; k++) {   // loop over diagonals

    // find the pivot row p
    p=k;
    for (i=k; i<A.Rows(); i++)  
      if (fabs(A.data[k][i]) > std::abs(A.data[k][p]))
	p=i;

    // swap rows in A
    for (j=k; j<A.Rows(); j++) 
      std::swap(A.data[j][p], A.data[j][k]);

    // swap rows in B
    for (j=0; j<B.Cols(); j++)
      std::swap(B.data[j][p], B.data[j][k]);
		
    // check for nonzero matrix diagonal
    if (fabs(A.data[k][k]) < STOL*Amax) {
      cerr << "Solve error: numerically singular matrix!\n";
      return 1;
    }

    // perform elimination using row k
    for (i=k+1; i<A.Rows(); i++)      // store multipliers in column below pivot
      A.data[k][i] /= A.data[k][k];
    for (j=k+1; j<A.Rows(); j++)      // loop over columns of A, to right of pivot 
      for (i=k+1; i<A.Rows(); i++)    // update rows in column
	A.data[j][i] -= A.data[k][i]*A.data[j][k];
    for (j=0; j<B.Cols(); j++)
      for (i=k+1; i<A.Rows(); i++)      // update entries in B
	B.data[j][i] -= A.data[k][i]*B.data[j][k];
  }

  // check for singularity at end (only need to check final diagonal entry)
  if (fabs(A.data[A.Rows()-1][A.Rows()-1]) < STOL*Amax) {
    cerr << "Solve error: numerically singular matrix!\n";
    return 1;
  }

  // perform Backward Substitution on result
  if (BackSub(A, X, B) != 0) {
    cerr << "Solve: error in BackSub call\n";
    return 1;
  }

  // return success
  return 0;
}

// solves a linear system A*X = B, returning X as a new Matrix
Matrix Solve(Matrix& A, Matrix& B) {

  // check that matrix sizes match
  if (A.Rows() != B.Rows() || A.Rows() != A.Cols()) {
    cerr << "Solve error, illegal matrix/vector dimensions\n";
    cerr << "  Matrix is " << A.Rows() << " x " << A.Cols() 
	 << ",  rhs is " << B.Rows() << " x " << B.Cols() << endl;
    return Matrix(0,0);
  } else {
    // create new Mat for output, and call existing Solve routine
    Matrix X(A.Rows(),B.Cols());
    if (Solve(A, X, B) != 0)
      cerr << "Solve: error in in-place Solve call\n";
    return X;
  }
}

// solves a linear system A*x = b, filling in the input vector<double> x
//    A and b are modified in this operation; x holds the result
int Solve(Matrix& A, vector<double>& x, vector<double>& b) {

  // check that matrix sizes match
  if (A.Rows() != b.size() || A.Rows() != A.Cols() ||
      x.size() != A.Cols()) {
    cerr << "Solve error, illegal matrix/vector dimensions\n";
    cerr << "  Matrix is " << A.Rows() << " x " << A.Cols() 
	 << ",  rhs is " << b.size() << " x 1"
	 << ",  solution is " << x.size() << " x 1\n";
    return 1;
  }
  
  // create temporary variables
  long int i, j, k, p;
  double Amax;

  // determine magnitude of entries in A (for singularity check later)
  Amax = InfNorm(A);

  // perform Gaussian elimination to convert A,B to an upper-triangular system
  for (k=0; k<A.Rows()-1; k++) {   // loop over diagonals

    // find the pivot row p
    p=k;
    for (i=k; i<A.Rows(); i++)  
      if (fabs(A.data[k][i]) > fabs(A.data[k][p]))
	p=i;

    // swap rows in A
    for (j=k; j<A.Rows(); j++) 
      std::swap(A.data[j][p], A.data[j][k]);

    // swap rows in b
    std::swap(b[p], b[k]);
		
    // check for nonzero matrix diagonal
    if (fabs(A.data[k][k]) < STOL*Amax) {
      cerr << "Solve error: numerically singular matrix!\n";
      return 1;
    }

    // perform elimination using row k
    for (i=k+1; i<A.Rows(); i++)      // store multipliers in column below pivot
      A.data[k][i] /= A.data[k][k];
    for (j=k+1; j<A.Rows(); j++)      // loop over columns of A, to right of pivot 
      for (i=k+1; i<A.Rows(); i++)    // update rows in column
	A.data[j][i] -= A.data[k][i]*A.data[j][k];
    for (i=k+1; i<A.Rows(); i++)      // update entries in b
      b[i] -= A.data[k][i]*b[k];
  }

  // check for singularity at end (only need to check final diagonal entry)
  if (fabs(A.data[A.Rows()-1][A.Rows()-1]) < STOL*Amax) {
    cerr << "Solve error: numerically singular matrix!\n";
    return 1;
  }

  // perform Backward Substitution on result
  if (BackSub(A, x, b) != 0) {
    cerr << "Solve: error in BackSub call\n";
    return 1;
  }

  // return success
  return 0;
}

// solves a linear system A*x = b, returning x as a new vector<double>
//    A and b are modified in this operation; x holds the result
vector<double> Solve(Matrix& A, vector<double>& b) {

  // check that matrix sizes match
  if (A.Rows() != b.size() || A.Rows() != A.Cols()) {
    cerr << "Solve error, illegal matrix/vector dimensions\n";
    cerr << "  Matrix is " << A.Rows() << " x " << A.Cols() 
	 << ",  rhs is " << b.size() << " x 1\n";
    return vector<double>(0);
  }

  // create output, call existing Solve routine and return
  vector<double> x = b;
  if (Solve(A, x, b) != 0)
    cerr << "Solve: error in in-place Solve call\n";
  return x;
}


// end of file
