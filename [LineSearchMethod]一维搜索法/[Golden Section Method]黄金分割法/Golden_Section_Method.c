#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
/*Golden Section Method]黄金分割法
Aouthor:Xipnig.Yu
Date:2020.05.27 
------------ 
min F(X)
s.t. X¡Ê[a,b] 
*/ 
double Function(double x);

int main(int argc, char *argv[]) {
	/*Data + End gap*/
	double Left = 4;
	double Right = 6;
	double deta = 0.001;
	/*Golden rate*/
	double Golden = 0.618;
	/*Operation*/
	double (*Func)(double) = Function;
	double LEFT = Left, RIGHT = Right;
	double X_left = Golden*LEFT+(1-Golden)*RIGHT, X_right = (1-Golden)*LEFT+Golden*RIGHT;
	double val_left = (*Func)(X_left),val_right = (*Func)(X_right);
	while((X_right-X_left)>deta){
		if (val_left<val_right){
			RIGHT = X_right;
			X_right = X_left;
			X_left = Golden*LEFT+(1-Golden)*RIGHT;
			val_left = (*Func)(X_left);
			val_right = val_left;
		}else{
			if(val_left>val_right){
				LEFT = X_left;
				X_left = X_right;
				X_right = (1-Golden)*LEFT+Golden*RIGHT;
				val_left = val_right;
				val_right = (*Func)(X_right);
			}else{
				LEFT = X_left;
				RIGHT = X_right;
				X_left = Golden*LEFT+(1-Golden)*RIGHT;
				X_right = (1-Golden)*LEFT+Golden*RIGHT;
				val_left = (*Func)(X_left);
				val_right = (*Func)(X_right);
			}
		}
	}
		double X = (X_right+X_left)/2;
		double Value = (*Func)(X);
		printf("X* = %lf, F(X*) = %lf",X,Value);
	return 0;
}

double Function(double x){
    return pow(x,3)-5*pow(x,2);//
}
