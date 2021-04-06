#include "matrix.h"
class Matrixsol:public Matrix  {
    public:
    static Matrix Gauss_Jacobi(Matrix A,Matrix b);
    static Matrix Gauss_Seidel_under_relaxation(Matrix A,Matrix b);
    static Matrix Gauss_Seidel_over_relaxation(Matrix A,Matrix b);
   
   /*
    Matrix Gauss_Elimination(Matrix A,Matrix b);
    Matrix Conjugate_Gradient(Matrix A,Matrix b);
    */
};

