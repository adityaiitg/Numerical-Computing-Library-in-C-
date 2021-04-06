#include"function_evaluation.h"
#define EPS 1e-10

double* function_evaluation::derivative(double* polynomial,int size){
    int new_size=size-1;
    double *p=new double[new_size];
    for(int i=0;i<new_size;i++){
        p[i]=polynomial[i+1]*(i+1);
    }
    return p;
}
double function_evaluation::evaluate(double x,double* polynomial,int size){
    double ans=0;
    double x_term=1;
    for(int i=0;i<size;i++){
        ans=ans+x_term*polynomial[i];
        x_term=x_term*x;
       // std::cout<<ans<<std::endl;
    }
    return ans;
}
void function_evaluation::print(double* polynomial,int size){
    for(int i=0;i<size-1;i++){
        std::cout<<polynomial[i]<<"x^"<<i<<"+";
    }
    std::cout<<polynomial[size-1]<<"x^"<<size-1<<std::endl;
}
//  x(n+1) = x(n)-f(xn)/f'(xn)
double function_evaluation::Newton_Raphson(double* polynomial,int size){
    print(polynomial,size);
    print(derivative(polynomial,size),size-1);
    double root_old=0;
    double root=root_old-evaluate(root_old,polynomial,size)/evaluate(root_old,derivative(polynomial,size),size-1);
    for(int i=0;i<100;i++){
   //     std::cout<<i<<" "<<root<<std::endl;
       root_old=root; 
       root=root_old-evaluate(root_old,polynomial,size)/evaluate(root_old,derivative(polynomial,size),size-1);
    }
    return root;
}



 
// evalaluation of integral by trapezoidal rule
double function_evaluation::trapezoidal(double a,double b,double* polynomial,int size)
{ 
    double n=100000;
    double step_size =(b-a)/n; 
    double s=0;
    for (double i = a; i < b; i=i+step_size) 
        s+=step_size*(evaluate(i,polynomial,size)+evaluate(i+step_size,polynomial,size))/2; 
    return s; 
} 




// Consider a differential equation 
// dy/dx=(x + y + xy) 
// double func(double x,double y){
double function_evaluation::dy_by_dx(double x,double y){ 
	return (x+y+x*y);
}
// Function for Euler formula 
double function_evaluation::forward_euler(double x_0,double y,double step_size,double x){
	double temp=0;
	while (x_0<x){
		y =y+step_size*dy_by_dx(x_0, y);
		x_0=x_0+step_size; 
	}
	return y;
}
/*
modifications
double backward_euler(double x_0,double y,double step_size,double x){
	double temp=0;
	while (x_0<x){
		y=y+step_size*dy_by_dx(x_0, y);
		x_0 = x_0 + step_size; 
	}
	return y;
}
*/