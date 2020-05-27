#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
/*Bisection Method¶þ·Ö·¨ 
Aouthor:Xipnig.Yu
Date:2020.05.27 
------------ 
min F(X)
s.t. X¡Ê[a,b] 
*/ 
double Function(double x);
double Func_derivative(double X,double step);

int main(int argc, char *argv[]) {
	/*Data + End gap*/
	double Left = 0;
	double Right = 6;
	double deta = 0.001;
	double step = 0.001; 
	/*Operation*/
	double (*Func)(double) = Function;
	double X_left = Left, X_right = Right, X_mid;
	while((X_right-X_left)>deta){
		X_mid = (X_right+X_left)/2;
		if (Func_derivative(X_mid,step)>0){
			X_right = X_mid;
		}else{
			if (Func_derivative(X_mid,step)<0){
				X_left = X_mid;
			}else{
				X_left = X_mid - step;
				X_right = X_mid + step;				
			}
		}
	}
		double X = (X_right+X_left)/2;
		double Value = (*Func)(X);
		printf("X* = %lf, F(X*) = %lf",X,Value);
	return 0;
}

double Function(double x){
    return pow(x,3)-5*pow(x,2)-x;
}

double Func_derivative(double X,double step){
	double (*Func)(double) = Function;
	return ((*Func)(X+step)-(*Func)(X))/step;
}
