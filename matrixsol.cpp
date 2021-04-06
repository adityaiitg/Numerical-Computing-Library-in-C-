#include "matrixsol.h"

#define EPS 1e-10
Matrix Matrixsol::Gauss_Jacobi(Matrix A,Matrix b){
    //A=D-L-U
    // CALCULATING INVERSE OF D
    double **p=new double*[A.get_Row()];
    for(int i=0;i<A.get_Row();i++){
        p[i]=new double[A.get_Col()];
    }
    for(int i=0;i<A.get_Row();i++){
      for(int j=0;j<A.get_Col();j++){
        if(i==j){
            p[i][j]=1/A(i,j);
        }else{
            p[i][j]=0;
        }
    }  
    }
    Matrix D_inv(p,A.get_Row(),A.get_Col());

    for(int i=0;i<A.get_Row();i++){
      for(int j=0;j<A.get_Col();j++){
        if(i>j){
            p[i][j]=-A(i,j);
        }else{
            p[i][j]=0;
        }
    }  
    }
    Matrix L(p,A.get_Row(),A.get_Col());
    for(int i=0;i<A.get_Row();i++){
      for(int j=0;j<A.get_Col();j++){
        if(i<j){
            p[i][j]=-A(i,j);
        }else{
            p[i][j]=0;
        }
    }  
    }
    Matrix U(p,A.get_Row(),A.get_Col());
    for(int i=0;i<A.get_Row();i++){
        delete p[i];
    }
    delete p;
    Matrix x_old(b.get_Row(),1); 
    std::cin>>x_old;
    Matrix T=D_inv*(L+U);
    Matrix C=D_inv*b;
    Matrix x=T*x_old+C;
    while(Matrix::dotProduct(x_old-x,x_old-x)>EPS){
        x_old=x;
        x=T*x_old+C;
    }
    return x;
}
Matrix Matrixsol::Gauss_Seidel_under_relaxation(Matrix A,Matrix b){
    double **p=new double*[A.get_Row()];
    for(int i=0;i<A.get_Row();i++){
        p[i]=new double[A.get_Col()];
    }
    for(int i=0;i<A.get_Row();i++){
      for(int j=0;j<A.get_Col();j++){
        if(i==j){
            p[i][j]=A(i,j);
        }else{
            p[i][j]=0;
        }
    }  
    }
    Matrix D(p,A.get_Row(),A.get_Col());

    for(int i=0;i<A.get_Row();i++){
      for(int j=0;j<A.get_Col();j++){
        if(i>j){
            p[i][j]=-A(i,j);
        }else{
            p[i][j]=0;
        }
    }  
    }
    Matrix L(p,A.get_Row(),A.get_Col());
    for(int i=0;i<A.get_Row();i++){
      for(int j=0;j<A.get_Col();j++){
        if(i<j){
            p[i][j]=-A(i,j);
        }else{
            p[i][j]=0;
        }
    }  
    }
    Matrix U(p,A.get_Row(),A.get_Col());
    for(int i=0;i<A.get_Row();i++){
        delete p[i];
    }
    delete p;
    Matrix x_old(A.get_Col(),1);    
    Matrix temp=(D-L).inverse();
    Matrix T=temp*U;
    Matrix C=temp*b;
    Matrix x=T*x_old+C;
    while(Matrix::dotProduct(x_old-x,x_old-x)>EPS){
        x_old=x;
        x=T*x_old+C;
    }
    return x;
}
// help from
// https://www.maa.org/press/periodicals/loci/joma/iterative-methods-for-solving-iaxi-ibi-the-sor-method#:~:text=A%20third%20iterative%20method%2C%20called,x(k%2B1).
Matrix Matrixsol::Gauss_Seidel_over_relaxation(Matrix A,Matrix b){
    double **p=new double*[A.get_Row()];
    for(int i=0;i<A.get_Row();i++){
        p[i]=new double[A.get_Col()];
    }
    for(int i=0;i<A.get_Row();i++){
      for(int j=0;j<A.get_Col();j++){
        if(i==j){
            p[i][j]=A(i,j);
        }else{
            p[i][j]=0;
        }
    }  
    }
    Matrix D(p,A.get_Row(),A.get_Col());

    for(int i=0;i<A.get_Row();i++){
      for(int j=0;j<A.get_Col();j++){
        if(i>j){
            p[i][j]=-A(i,j);
        }else{
            p[i][j]=0;
        }
    }  
    }
    Matrix L(p,A.get_Row(),A.get_Col());
    for(int i=0;i<A.get_Row();i++){
      for(int j=0;j<A.get_Col();j++){
        if(i<j){
            p[i][j]=-A(i,j);
        }else{
            p[i][j]=0;
        }
    }  
    }
    Matrix U(p,A.get_Row(),A.get_Col());
    for(int i=0;i<A.get_Row();i++){
        delete p[i];
    }
    delete p;
    Matrix x_old(A.get_Col(),1);    
    Matrix temp=(D-L).inverse();
    Matrix T=temp*U;
    Matrix C=temp*b;
    Matrix x=T*x_old+C;
    while(Matrix::dotProduct(x_old-x,x_old-x)>EPS){
        x_old=x;
        x=1.5*T*x_old+x_old;
    }
    return x;
}

Matrixsol::Matrix Gauss_Elimination(Matrix A,Matrix b){
    // Already implemented in matrix
return Matrix::solve(A,b);
}
/*
#include <iostream>
#include <vector>
using namespace std;

//Naive matrix multiplications!! Anyways, to compile g++ -std=c++14 -o main conjugate_gradient.cpp
//gionuno

vector<double> conjugate_gradient(const vector<vector<double> > & A, const vector<double> & b,int T)
{
	int N = b.size();
	vector<double> r(N,0.0);
	vector<double> p(N,0.0);
	vector<double> x(N,0.0);
	for(int i=0;i<N;i++)
		p[i] = r[i] = b[i];
	int t = 0;
	while(t < T)
	{
		double rtr = 0.0;
		double ptAp = 0.0;
		for(int i=0;i<N;i++)
			rtr += r[i]*r[i];
		for(int i=0;i<N;i++)
			for(int j=0;j<N;j++)
				ptAp += A[i][j]*p[i]*p[j];
		double alpha = rtr / (ptAp + 1e-10);
		vector<double> rn(N,0.0);
		for(int i=0;i<N;i++)
		{
			x[i] += alpha * p[i];
			rn[i] = r[i];
			for(int j=0;j<N;j++)
				rn[i] -= alpha*A[i][j]*p[j];
		}
		double rntrn = 0.0;
		for(int i=0;i<N;i++)
			rntrn += rn[i]*rn[i];
		if(rntrn < 1e-10) break;
		double beta = rntrn / rtr;
		for(int i=0;i<N;i++)
		{
			p[i] = beta*p[i] + rn[i];
			r[i] = rn[i];
		} 
		t++;
	}
	return x;
}

int main()
{
	//Only well defined for symmetric positive def matrices.
	vector<vector<double> > A(3,vector<double>(3,0.0));
	A[0][0] =  7.0; A[0][1] =  3.0; A[0][2] =  1.0;
	A[1][0] =  3.0;	A[1][1] =  7.0;
	A[2][0] =  1.0;                 A[2][2] = 10.0;
	vector<double> b(3,0.0);
	b[0] = 1.0;
	b[1] = -5.0;
	b[2] = 2.0;
	vector<double> x = conjugate_gradient(A,b,1000);
	for(int i=0;i<b.size();i++)
		cout << x[i] << endl;
	return 0;
}


*/