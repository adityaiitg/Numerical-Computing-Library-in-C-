//  matrix.cpp

#include <stdexcept>
#include "matrix.h"

#define EPS 1e-10

using std::ostream;  using std::istream;  using std::endl;
using std::domain_error;

/* PUBLIC MEMBER FUNCTIONS
 ********************************/
// help for member initialising
/*
Member initialization in constructors
When a constructor is used to initialize other members, these other members can be initialized directly, without resorting to statements in its body. This is done by inserting, before the constructor's body, a colon (:) and a list of initializations for class members. For example, consider a class with the following declaration:

1
2
3
4
5
6
class Rectangle {
    int width,height;
  public:
    Rectangle(int,int);
    int area() {return width*height;}
};


The constructor for this class could be defined, as usual, as:

 
Rectangle::Rectangle (int x, int y) { width=x; height=y; }


But it could also be defined using member initialization as:

 
Rectangle::Rectangle (int x, int y) : width(x) { height=y; }


Or even:

Rectangle::Rectangle (int x, int y) : width(x), height(y) { }


Note how in this last case, the constructor does nothing else than initialize its members, hence it has an empty function body.

For members of fundamental types, it makes no difference which of the ways above the constructor is defined, because they are not initialized by default, but for member objects (those whose type is a class), if they are not initialized after the colon, they are default-constructed.

Default-constructing all members of a class may or may always not be convenient: in some cases, this is a waste (when the member is then reinitialized otherwise in the constructor), but in some other cases, default-construction is not even possible (when the class does not have a default constructor). In these cases, members shall be initialized in the member initialization list. For example:

// member initialization
#include <iostream>
using namespace std;

class Circle {
    double radius;
  public:
    Circle(double r) : radius(r) { }
    double area() {return radius*radius*3.14159265;}
};

class Cylinder {
    Circle base;
    double height;
  public:
    Cylinder(double r, double h) : base (r), height(h) {}
    double volume() {return base.area() * height;}
};

int main () {
  Cylinder foo (10,20);

  cout << "foo's volume: " << foo.volume() << '\n';
  return 0;
}
foo's volume: 6283.19
 Edit & Run


In this example, class Cylinder has a member object whose type is another class (base's type is Circle). Because objects of class Circle can only be constructed with a parameter, Cylinder's constructor needs to call base's constructor, and the only way to do this is in the member initializer list.

These initializations can also use uniform initializer syntax, using braces {} instead of parentheses ():

 
Cylinder::Cylinder (double r, double h) : base{r}, height{h} { }
*/
Matrix::Matrix(int rows=1, int cols=1) : Number_of_rows(rows), Number_of_cols(cols)
{
    SpaceAllocator();
    for (int i = 0; i < Number_of_rows; ++i) {
        for (int j = 0; j < Number_of_cols; ++j) {
            p[i][j] = 0;
        }
    }
}


Matrix::Matrix(double** a, int rows, int cols) : Number_of_rows(rows), Number_of_cols(cols)
{
    SpaceAllocator();
    for (int i = 0; i < Number_of_rows; ++i) {
        for (int j = 0; j < Number_of_cols; ++j) {
            p[i][j] = a[i][j];
        }
    }
}
/*
Matrix::Matrix() : Number_of_rows(1), Number_of_cols(1)
{
    SpaceAllocator();
    p[0][0] = 0;
}

*/
Matrix::~Matrix()
{
    for (int i = 0; i < Number_of_rows; ++i) {
        delete[] p[i];
    }
    delete[] p;
}

Matrix::Matrix(const Matrix& m) : Number_of_rows(m.Number_of_rows), Number_of_cols(m.Number_of_cols)
{
    SpaceAllocator();
    for (int i = 0; i < Number_of_rows; ++i) {
        for (int j = 0; j < Number_of_cols; ++j) {
            p[i][j] = m.p[i][j];
        }
    }
}

Matrix& Matrix::operator=(const Matrix& m)
{
    if (this == &m) {
        return *this;
    }

    if (Number_of_rows != m.Number_of_rows || Number_of_cols != m.Number_of_cols) {
        for (int i = 0; i < Number_of_rows; ++i) {
            delete[] p[i];
        }
        delete[] p;

        Number_of_rows = m.Number_of_rows;
        Number_of_cols = m.Number_of_cols;
        SpaceAllocator();
    }

    for (int i = 0; i < Number_of_rows; ++i) {
        for (int j = 0; j < Number_of_cols; ++j) {
            p[i][j] = m.p[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator+=(const Matrix& m)
{
    if(this->Number_of_cols!=m.Number_of_cols||this->Number_of_rows!=m.Number_of_rows){
        throw "Left and Right size diffent sizes!";
    }
    for (int i = 0; i < Number_of_rows; ++i) {
        for (int j = 0; j < Number_of_cols; ++j) {
            p[i][j] += m.p[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& m)
{
    if(this->Number_of_cols!=m.Number_of_cols||this->Number_of_rows!=m.Number_of_rows){
        throw "Left and Right size diffent sizes!";
    }
    for (int i = 0; i < Number_of_rows; ++i) {
        for (int j = 0; j < Number_of_cols; ++j) {
            p[i][j] -= m.p[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator*=(const Matrix& m)
{
    Matrix temp(Number_of_rows, m.Number_of_cols);
    for (int i = 0; i < temp.Number_of_rows;i++) {
        for (int j = 0; j < temp.Number_of_cols;j++) {
            for (int k = 0; k < Number_of_cols; k++) {
                temp.p[i][j] += (p[i][k] * m.p[k][j]);
            }
        }
    }
    return (*this = temp);
}
// element wise multiplication
Matrix& Matrix::operator*=(double num)
{
    for (int i = 0; i < Number_of_rows; ++i) {
        for (int j = 0; j < Number_of_cols; ++j) {
            p[i][j] *= num;
        }
    }
    return *this;
}
// element wise division
Matrix& Matrix::operator/=(double num)
{
    for (int i = 0; i < Number_of_rows; ++i) {
        for (int j = 0; j < Number_of_cols; ++j) {
            p[i][j] /= num;
        }
    }
    return *this;
}

// element wise powers
Matrix Matrix::operator^(int num)
{
    Matrix temp(*this);
    return expHelper(temp, num);
}
// Swap the rows
void Matrix::swapRows(int r1, int r2)
{
    double *temp = p[r1];
    p[r1] = p[r2];
    p[r2] = temp;
}
// Transpose
Matrix Matrix::transpose()
{
    Matrix ret(Number_of_cols, Number_of_rows);
    for (int i = 0; i < Number_of_rows; ++i) {
        for (int j = 0; j < Number_of_cols; ++j) {
            ret.p[j][i] = p[i][j];
        }
    }
    return ret;
}


/* STATIC CLASS FUNCTIONS
 ********************************/
int Matrix::get_Col(){
    return Number_of_cols; 
}
int Matrix::get_Row(){
    return Number_of_rows;
}
Matrix Matrix::createIdentity(int size)
{
    Matrix temp(size, size);
    for (int i = 0; i < temp.Number_of_rows; ++i) {
        for (int j = 0; j < temp.Number_of_cols; ++j) {
            if (i == j) {
                temp.p[i][j] = 1;
            } else {
                temp.p[i][j] = 0;
            }
        }
    }
    return temp;
}

Matrix Matrix::solve(Matrix A, Matrix b)
{
    // Gaussian elimination
    for (int i = 0; i < A.Number_of_rows; ++i) {
        if (A.p[i][i] == 0) {
            // pivot 0 - throw error
            throw domain_error("Error: the coefficient matrix has 0 as a pivot. Please fix the input and try again.");
        }
        for (int j = i + 1; j < A.Number_of_rows; ++j) {
            for (int k = i + 1; k < A.Number_of_cols; ++k) {
                A.p[j][k] -= A.p[i][k] * (A.p[j][i] / A.p[i][i]);
                if (A.p[j][k] < EPS && A.p[j][k] > -1*EPS)
                    A.p[j][k] = 0;
            }
            b.p[j][0] -= b.p[i][0] * (A.p[j][i] / A.p[i][i]);
            if (A.p[j][0] < EPS && A.p[j][0] > -1*EPS)
                A.p[j][0] = 0;
            A.p[j][i] = 0;
        }
    }

    // Back substitution
    Matrix x(b.Number_of_rows, 1);
    x.p[x.Number_of_rows - 1][0] = b.p[x.Number_of_rows - 1][0] / A.p[x.Number_of_rows - 1][x.Number_of_rows - 1];
    if (x.p[x.Number_of_rows - 1][0] < EPS && x.p[x.Number_of_rows - 1][0] > -1*EPS)
        x.p[x.Number_of_rows - 1][0] = 0;
    for (int i = x.Number_of_rows - 2; i >= 0; --i) {
        int sum = 0;
        for (int j = i + 1; j < x.Number_of_rows; ++j) {
            sum += A.p[i][j] * x.p[j][0];
        }
        x.p[i][0] = (b.p[i][0] - sum) / A.p[i][i];
        if (x.p[i][0] < EPS && x.p[i][0] > -1*EPS)
            x.p[i][0] = 0;
    }

    return x;
}

Matrix Matrix::bandSolve(Matrix A, Matrix b, int k)
{
    // optimized Gaussian elimination
    int bandsBelow = (k - 1) / 2;
    for (int i = 0; i < A.Number_of_rows; ++i) {
        if (A.p[i][i] == 0) {
            // pivot 0 - throw exception
            throw domain_error("Error: the coefficient matrix has 0 as a pivot. Please fix the input and try again.");
        }
        for (int j = i + 1; j < A.Number_of_rows && j <= i + bandsBelow; ++j) {
            int k = i + 1;
            while (k < A.Number_of_cols && A.p[j][k]) {
                A.p[j][k] -= A.p[i][k] * (A.p[j][i] / A.p[i][i]);
                k++;
            }
            b.p[j][0] -= b.p[i][0] * (A.p[j][i] / A.p[i][i]);
            A.p[j][i] = 0;
        }
    }

    // Back substitution
    Matrix x(b.Number_of_rows, 1);
    x.p[x.Number_of_rows - 1][0] = b.p[x.Number_of_rows - 1][0] / A.p[x.Number_of_rows - 1][x.Number_of_rows - 1];
    for (int i = x.Number_of_rows - 2; i >= 0; --i) {
        int sum = 0;
        for (int j = i + 1; j < x.Number_of_rows; ++j) {
            sum += A.p[i][j] * x.p[j][0];
        }
        x.p[i][0] = (b.p[i][0] - sum) / A.p[i][i];
    }

    return x;
}

// functions on VECTORS
double Matrix::dotProduct(Matrix a, Matrix b)
{
    double sum = 0;
    for (int i = 0; i < a.Number_of_rows; ++i) {
        sum += (a(i, 0) * b(i, 0));
    }
    return sum;
}

// functions on AUGMENTED matrices
Matrix Matrix::augment(Matrix A, Matrix B)
{
    Matrix AB(A.Number_of_rows, A.Number_of_cols + B.Number_of_cols);
    for (int i = 0; i < AB.Number_of_rows; ++i) {
        for (int j = 0; j < AB.Number_of_cols; ++j) {
            if (j < A.Number_of_cols)
                AB(i, j) = A(i, j);
            else
                AB(i, j) = B(i, j - B.Number_of_cols);
        }
    }
    return AB;
}

Matrix Matrix::gaussianEliminate()
{
    Matrix Ab(*this);
    int rows = Ab.Number_of_rows;
    int cols = Ab.Number_of_cols;
    int Acols = cols - 1;

    int i = 0; // row tracker
    int j = 0; // column tracker

    // iterate through the rows
    while (i < rows)
    {
        // find a pivot for the row
        bool pivot_found = false;
        while (j < Acols && !pivot_found)
        {
            if (Ab(i, j) != 0) { // pivot not equal to 0
                pivot_found = true;
            } else { // check for a possible swap
                int max_row = i;
                double max_val = 0;
                for (int k = i + 1; k < rows; ++k)
                {
                    if(Ab(k,j)!=0){
                        max_row=k;
                        max_val=Ab(k,j)>0?Ab(k,j):-1*Ab(k,j);
                    }
                }
             /*   for (int k = i + 1; k < rows; ++k)
                {
                    double cur_abs = Ab(k, j) >= 0 ? Ab(k, j) : -1 * Ab(k, j);
                    if (cur_abs > max_val)
                    {
                        max_row = k;
                        max_val = cur_abs;
                    }
                } */
                if (max_row != i) {
                    Ab.swapRows(max_row, i);
                    pivot_found = true;
                } else {
                    j++;
                }
            }
        }

        // perform elimination as normal if pivot was found
        if (pivot_found)
        {
            for (int t = i + 1; t < rows; ++t) {
                for (int s = j + 1; s < cols; ++s) {
                    Ab(t, s) = Ab(t, s) - Ab(i, s) * (Ab(t, j) / Ab(i, j));
                    // in case of error in calculation of floating numbers forcing the values to 0.
                    if (Ab(t, s) < EPS && Ab(t, s) > -1*EPS)
                        Ab(t, s) = 0;
                }
                Ab(t, j) = 0;
            }
        }
        i++;
        j++;
    }
    return Ab;
}


int Matrix::deteminant()
{
    Matrix Ab(*this);
    int rows = Ab.Number_of_rows;
    int cols = Ab.Number_of_cols;
    if(rows!=cols){
        std::cout<<"Not Defined. Deteminant are defined only for Squared Matrix"<<endl;
        return 0;
    }
    int Acols = cols - 1;
    int i = 0; // row tracker
    int j = 0; // column tracker
    int sign=1;
    // iterate through the rows
    while (i < rows)
    {
        // find a pivot for the row
        bool pivot_found = false;
        while (j < Acols && !pivot_found)
        {
            if (Ab(i, j) != 0) { // pivot not equal to 0
                pivot_found = true;
            } else { // check for a possible swap
                int max_row = i;
                double max_val = 0;
                for (int k = i + 1; k < rows; ++k)
                {
                    if(Ab(k,j)!=0){
                        max_row=k;
                        max_val=Ab(k,j)>0?Ab(k,j):-1*Ab(k,j);
                        break;
                    }
                }
                if (max_row != i) {
                    Ab.swapRows(max_row, i);
                    sign=-1*sign;
                    pivot_found = true;
                } else {
                    j++;
                }
            }
        }

        // perform elimination as normal if pivot was found
        if (pivot_found)
        {
            for (int t = i + 1; t < rows; ++t) {
                for (int s = j + 1; s < cols; ++s) {
                    Ab(t, s) = Ab(t, s) - Ab(i, s) * (Ab(t, j) / Ab(i, j));
                    // in case of error in calculation of floating numbers forcing the values to 0.
                    if (Ab(t, s) < EPS && Ab(t, s) > -1*EPS)
                        Ab(t, s) = 0;
                }
                Ab(t, j) = 0;
            }
        }
        i++;
        j++;
    }
    double ans=1;
    for(int i=0;i<rows;i++){
        ans=ans*Ab(i,i);
    }
    ans=ans*sign;
    return ans;
}

Matrix Matrix::rowReduceFromGaussian()
{
    Matrix R(*this);
    int rows = R.Number_of_rows;
    int cols = R.Number_of_cols;

    int i = rows - 1; // row tracker
    int j = cols - 2; // column tracker

    // iterate through every row
    while (i >= 0)
    {
        // find the pivot column
        int k = j - 1;
        while (k >= 0) {
            if (R(i, k) != 0)
                j = k;
            k--;
        }

        // zero out elements above pivots if pivot not 0
        if (R(i, j) != 0) {
       
            for (int t = i - 1; t >= 0; --t) {
                for (int s = 0; s < cols; ++s) {
                    if (s != j) {
                        R(t, s) = R(t, s) - R(i, s) * (R(t, j) / R(i, j));
                        if (R(t, s) < EPS && R(t, s) > -1*EPS)
                            R(t, s) = 0;
                    }
                }
                R(t, j) = 0;
            }

            // divide row by pivot
            for (int k = j + 1; k < cols; ++k) {
                R(i, k) = R(i, k) / R(i, j);
                if (R(i, k) < EPS && R(i, k) > -1*EPS)
                    R(i, k) = 0;
            }
            R(i, j) = 1;
        }
        i--;
        j--;
    }

    return R;
}

void Matrix::readSolutionsFromRREF(ostream& os)
{
    Matrix R(*this);

    // print number of solutions
    bool hasSolutions = true;
    bool doneSearching = false;
    int i = 0;
    while (!doneSearching && i < Number_of_rows)
    {
        bool allZeros = true;
        for (int j = 0; j < Number_of_cols - 1; ++j) {
            if (R(i, j) != 0)
                allZeros = false;
        }
        if (allZeros && R(i, Number_of_cols - 1) != 0) {
            hasSolutions = false;
            os << "NO SOLUTIONS" << endl << endl;
            doneSearching = true;
        } else if (allZeros && R(i, Number_of_cols - 1) == 0) {
            os << "INFINITE SOLUTIONS" << endl << endl;
            doneSearching = true;
        } else if (Number_of_rows < Number_of_cols - 1) {
            os << "INFINITE SOLUTIONS" << endl << endl;
            doneSearching = true;
        }
        i++;
    }
    if (!doneSearching)
        os << "UNIQUE SOLUTION" << endl << endl;

    // get solutions if they exist
    if (hasSolutions)
    {
        Matrix particular(Number_of_cols - 1, 1);
        Matrix special(Number_of_cols - 1, 1);

        for (int i = 0; i < Number_of_rows; ++i) {
            bool pivotFound = false;
            bool specialCreated = false;
            for (int j = 0; j < Number_of_cols - 1; ++j) {
                if (R(i, j) != 0) {
                    // if pivot variable, add b to particular
                    if (!pivotFound) {
                        pivotFound = true;
                        particular(j, 0) = R(i, Number_of_cols - 1);
                    } else { // otherwise, add to special solution
                        if (!specialCreated) {
                            special = Matrix(Number_of_cols - 1, 1);
                            specialCreated = true;
                        }
                        special(j, 0) = -1 * R(i, j);
                    }
                }
            }
            os << "Special solution:" << endl << special << endl;
        }
        os << "Particular solution:" << endl << particular << endl;
    }
}

Matrix Matrix::inverse()
{
    Matrix I = Matrix::createIdentity(Number_of_rows);
    Matrix AI = Matrix::augment(*this, I);
    Matrix U = AI.gaussianEliminate();
    Matrix IAInverse = U.rowReduceFromGaussian();
    Matrix AInverse(Number_of_rows, Number_of_cols);
    for (int i = 0; i < AInverse.Number_of_rows; ++i) {
        for (int j = 0; j < AInverse.Number_of_cols; ++j) {
            AInverse(i, j) = IAInverse(i, j + Number_of_cols);
        }
    }
    return AInverse;
}


/* PRIVATE HELPER FUNCTIONS
 ********************************/

void Matrix::SpaceAllocator()
{
    p = new double*[Number_of_rows];
    for (int i = 0; i < Number_of_rows; ++i) {
        p[i] = new double[Number_of_cols];
    }
}

Matrix Matrix::expHelper(const Matrix& m, int num)
{
    if (num == 0) { 
        return createIdentity(m.Number_of_rows);
    } else if (num == 1) {
        return m;
    } else if (num % 2 == 0) {  // num is even
        return expHelper(m * m, num/2);
    } else {                    // num is odd
        return m * expHelper(m * m, (num-1)/2);
    }
}

/* NON-MEMBER FUNCTIONS
 ********************************/

Matrix operator+(const Matrix& m1, const Matrix& m2)
{
    Matrix temp(m1);
    return (temp += m2);
}

Matrix operator-(const Matrix& m1, const Matrix& m2)
{
    Matrix temp(m1);
    return (temp -= m2);
}

Matrix operator*(const Matrix& m1, const Matrix& m2)
{
    Matrix temp(m1);
    return (temp *= m2);
}

Matrix operator*(const Matrix& m, double num)
{
    Matrix temp(m);
    return (temp *= num);
}

Matrix operator*(double num, const Matrix& m)
{
    return (m * num);
}

Matrix operator/(const Matrix& m, double num)
{
    Matrix temp(m);
    return (temp /= num);
}

ostream& operator<<(ostream& os, const Matrix& m)
{
    for (int i = 0; i < m.Number_of_rows; ++i) {
        os << m.p[i][0];
        for (int j = 1; j < m.Number_of_cols; ++j) {
            os << " " << m.p[i][j];
        }
        os << endl;
    }
    return os;
}

istream& operator>>(istream& is, Matrix& m)
{
    for (int i = 0; i < m.Number_of_rows; ++i) {
        for (int j = 0; j < m.Number_of_cols; ++j) {
            is >> m.p[i][j];
        }
    }
    return is;
}