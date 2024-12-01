#include <iostream>
#include <iomanip>
#include "matclpro/cmatrix"
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
using namespace techsoft;
using namespace std;

double const k = 4500;
double const m = 2000;
double const pi = 4*atan(1);

//función normmax
double normmax(matrix<double> p){
    int n = p.rowno();
    double sum = 0;
    for(int i=0;i<n;i++){
        sum += fabs(p(i,0));
    }
return sum;
}
// función que devuelve los indices del mayor elemento fuera de la diagonal

matrix<double> ind(matrix<double> A){
double max = 0;
int n = 4;
matrix<double> ind(1,2);ind.null();
int im = 0,jm = 0;
 for(int i=0; i<n; ++i){
        for(int j=0; j<n;++j){
            if (i<j){


            if(fabs(A(i,j))>max){
                max = fabs(A(i,j));
                im = i;
                jm = j;
                
            }

        }
        
    }

}
ind(0,0) = im;
ind(0,1) = jm;
return ind;
}
// función que diagonaliza
void Adiag( matrix<double> (*ind)(matrix<double>), matrix<double> A, matrix<double>& Ad,int& nit, double tol, matrix<double>& U){

int im = 0, jm = 0;
double theta = 0;
int it = 0;
matrix<double> c(1,2);c.null();
U.unit();


do{
  
   c = ind(A);
   im = c(0,0);
   jm = c(0,1);
 
if(A(im,im)!= A(jm,jm)){
    theta = 0.5*atan(2*A(im,jm)/(A(im,im)-A(jm,jm)));
}
else{
    theta = pi/4;
}

matrix<double> R(4,4);R.unit();
R(im,im) = cos(theta);
R(jm,jm) = cos(theta);
R(im,jm) = -sin(theta);
R(jm,im) = sin(theta);

A = ~(R) * A * R;
U = U*R;
it += 1;
A(im,jm) = A(jm,im);
    c = ind(A);
   im = c(0,0);
   jm = c(0,1);

}while(fabs(A(im,jm))>tol);

 Ad = A;
nit = it;



}

int main(){
    int n = 4;
//creamos el fichero de salida
string infile = "Kr.txt";
    ofstream ff(infile);

    if (ff.is_open()){


    
    matrix<double> Kr(4,4);Kr.null();
    matrix<double> K(4,4);K.null();
    matrix<double> M(4,4);M.null();


// defino K y M
for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){

    if(i==j){
        K(i,j) = 2*k;
        M(i,j) = m;
    }
    else if((i == j-1)||(i == j+1)){
        K(i,j) = -k;
    }  
    }
     
}
K(n-1,n-1) = k;


 Kr = !(M) * K;
 cout<<"A matriz Kr é: "<<endl<<Kr<<endl;

ff<< 4<<" "<<4<<endl;
        ff<< Kr; 





ff.close();

matrix<double> Ad(4,4);Ad.unit();
matrix<double> U(4,4);U.unit();
int nit = 0;
double tol = pow(10,-10);

 Adiag(&ind,Kr,Ad,nit,tol,U);
cout<<"A matriz diagonalizada é: "<<endl<< Ad<<endl;
cout<< "A matriz de autovetor é: "<<endl<<U<<endl;
cout<< "O número de iterações é: "<<nit<<endl;

//prueba
double err =0;
err = normmax(~U*Kr*U- Ad);
cout<<"O error é: "<<err<<endl;

    }
return 0; 
}