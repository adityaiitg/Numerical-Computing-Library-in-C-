#include<global.h>
template <class T>
class ComplexNumbers {
    T r;
    T i;
 public:
    T getr(){
        return this->r;
}
    T geti(){
        return this->i;
    }
    ComplexNumbers(T x,T y){
    this->r=x;
    this->i=y;
}
    void plus(ComplexNumbers c){
        r+=c.getr();
        i+=c.geti();
    }
    void multiply(ComplexNumbers c){
        T real=r;
        r=c.getr()*real-c.geti()*i;
        i=c.geti()*real+i*c.getr();
    }
    void print(){
        if(i>0){
        cout<<r<<" + i"<<i<<endl;
        }else{
          cout<<r<<" - i"<<abs(i)<<endl;  
        }
    } 
};