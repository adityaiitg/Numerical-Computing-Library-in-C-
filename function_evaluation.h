#include "matrixsol.h"
#include<iostream>
#define EPS 1e-10
//  x(n+1) = x(n)-f(xn)/f'(xn)
class function_evaluation
{
public:
static double* derivative(double* polynomial,int size);
static double evaluate(double x,double* polynomial,int size);
static double Newton_Raphson(double* polynomial,int size);  
static void print(double* polynomial,int size);
static double trapezoidal(double a,double b,double* polynomial,int size);
static double dy_by_dx(double x,double y);
static double forward_euler(double x_0,double y,double step_size,double x);
};

