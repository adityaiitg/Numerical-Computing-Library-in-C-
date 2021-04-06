#include<iostream>
#include "function_evaluation.h"
using namespace std;
int main(){
    /*
    Matrix A(3,3);
    cin>>A;
    Matrix b(3,1);
    cin>>b;
    // cout<<A*b;
    cout<<"Gauss_Jacobi(A,b)"<<Matrixsol::Gauss_Jacobi(A,b);
    cout<<"Gauss_Seidel(A,b)"<<Matrixsol::Gauss_Seidel_under_relaxation(A,b);
    cout<<"Gauss_Seidel(A,b)"<<Matrixsol::Gauss_Seidel_over_relaxation(A,b);
    */
  /* int size=0;
   cin>>size;
    double* A=new double[size];
    for(int i=0;i<size;i++){
        cin>>A[i];
    }
     cout<<function_evaluation::evaluate(0,A,size)<<endl;
    cout<<function_evaluation::Newton_Raphson(A,size);
    */
   int size=0;
   cin>>size;
    double* A=new double[size];
    for(int i=0;i<size;i++){
        cin>>A[i];
    }
     cout<<function_evaluation::trapezoidal(0,1,A,size)<<endl;
}

/*
// Driver program 
int main() 
{ 
	// Initial Values 
	float x0 = 0; 
	float y0 = 1; 
	float h = 0.025; 

	// Value of x at which we need approximation 
	float x = 0.1; 

	euler(x0, y0, h, x); 
	return 0; 
}*/