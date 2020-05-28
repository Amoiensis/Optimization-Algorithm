#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include "matrix.h"
/*Univariate_Search_Method
多元函数[高/二次]-坐标轮换法 
Aouthor:Xipnig.Yu
Date:2020.05.27 
------------ 
min F(X)
s.t. X∈R^n
*/ 
double Function(Matirx* mat_X);
Matirx* Func_derivative(Matirx* mat_X, double step);
double Function_sub(double lamda, Matirx* _mat_X, Matirx* _mat_gradient);
double Func_derivative_sub(double X,double step,Matirx* _mat_X, Matirx* _mat_gradient);
double Bisection_Method(Matirx* _mat_X, Matirx* _mat_gradient) ;

int main(int argc, char *argv[]) {
	/*Data + End gap*/
	double (*Func)(Matirx*) = Function;
//	MATRIX_TYPE _mat_X [1][2] = {2,1}; //Two-dimensional function
	MATRIX_TYPE _mat_X [1][3] = {2,1,1}; //Three-dimensional function 
	double deta = 0.0018;
	double step_gradient = 0.01;
	double step_lamda = 0.01;
	/*Generate mat_X*/
	int row = sizeof(_mat_X) / sizeof(_mat_X[0]);
	int column = sizeof(_mat_X[0]) / sizeof(_mat_X[0][0]);
	Matirx*  mat_X = Matrix_gen(row,column,_mat_X);
	MATRIX_TYPE Mat_norm = 10*deta;
	Matirx* mat_grandient_last = NULL,* mat_X_gap = NULL,* mat_temp = NULL;
	Matirx* mat_direction = Matrix_copy(mat_X);
	double beta = 0;
	int times = 0,i = 0;
	/*Operation*/
	while(Mat_norm > deta){/*使用共轭梯度法，理论上只需要3次即收敛*/ 
		i = times%(column);
		if (times==0){
			mat_direction->data[i] = 1;
		}else{
			mat_direction->data[i-1] = 0;
			mat_direction->data[i] = 1;
		}
		step_lamda = Bisection_Method(mat_X,mat_direction);
		mat_temp = mat_X;
		M_print(mat_X);
		mat_X = M_add_sub(1,mat_temp,step_lamda,mat_direction);
		M_print(mat_X);
		mat_X_gap = M_add_sub(1,mat_X,1,mat_temp);
		Mat_norm = M_norm(mat_X_gap,2);
		M_free(mat_temp);
		times ++;
	}
	printf("X* =\n");
	M_print(mat_X);
	printf("F(X*) = %lf",(*Func)(mat_X));
	system("pause");
	return 0;
}

double Function(Matirx* mat_X){/*F(X)
X为参数，多元函数； 
*/ 
	/*Y=XAX'+BX+C*/
	/*mat_A*/
//	MATRIX_TYPE _mat_A [2][2] = {1,0,0,1}; //Two-dimensional function
	MATRIX_TYPE _mat_A [3][3] = {2,-2,0,0,3,-2,0,0,2}; //Three-dimensional function 
	/*-具体生成A-*/
	int row = sizeof(_mat_A) / sizeof(_mat_A[0]);
	int column = sizeof(_mat_A[0]) / sizeof(_mat_A[0][0]);
	Matirx*  mat_A = Matrix_gen(row,column,_mat_A);
	/*mat_B*/
//	MATRIX_TYPE _mat_B [2][1] = {-2,-4}; //Two-dimensional function
	MATRIX_TYPE _mat_B [3][1] = {-2,-2,-2}; //Three-dimensional function 
	/*-具体生成B-*/
	row = sizeof(_mat_B) / sizeof(_mat_B[0]);
	column = sizeof(_mat_B[0]) / sizeof(_mat_B[0][0]);
	Matirx*  mat_B = Matrix_gen(row,column,_mat_B);
	/*XAX'+BX*/
	Matirx* mat_Result = M_mul(M_mul(mat_X,mat_A),M_T(mat_X));
	Matirx*  mat_temp = mat_Result;
	mat_Result = M_add_sub(1,mat_temp,-1,M_mul(mat_X,mat_B));
	M_free(mat_temp);
	/*XAX'+BX + C*/ 
	double Result = (mat_Result->data)[0] -2;
	M_free(mat_Result);
	M_free(mat_A);
	M_free(mat_B);
    return Result;
}

Matirx* Func_derivative(Matirx* mat_X, double step){/*F'(X)
X为参数，多元函数的一阶导数；
*/ 
	double (*Func)(Matirx*) = Function;
	int Dimension = (mat_X->row)*(mat_X->column);
	MATRIX_TYPE* _mat_grandient = (MATRIX_TYPE*)malloc(Dimension*sizeof(MATRIX_TYPE));
	Matirx* mat_tempX = Matrix_gen(mat_X->row,mat_X->column,mat_X->data);
	double Val_ini = (*Func)(mat_tempX);
	int i;
	/*Find Grandient*/
	for(i=0;i<Dimension;i++){
		if (i==0){
			mat_tempX->data[i] +=  step;
		}else{
			mat_tempX->data[i-1] -=  step;
			mat_tempX->data[i] +=  step;
		}
		_mat_grandient[i] = ((*Func)(mat_tempX)-Val_ini)/step;
	}
	/*Generate mat_grandient*/
	int row = mat_X->row;
	int column = mat_X->column;
	Matirx*  mat_grandient = Matrix_gen(row,column,_mat_grandient);
	return mat_grandient;
}

double Function_sub(double lamda, Matirx* _mat_X, Matirx* _mat_gradient){/*F'(x_k+lamda*d_k) 
以lamda为变量的函数的导数 
*/ 
	Matirx* mat_X = M_add_sub(1,_mat_X,lamda,_mat_gradient);
	/*Y=XAX'+BX+C*/
	/*mat_A*/
//	MATRIX_TYPE _mat_A [2][2] = {1,0,0,1}; //Two-dimensional function
	MATRIX_TYPE _mat_A [3][3] = {2,-2,0,0,3,-2,0,0,2}; //Three-dimensional function 
	/*-具体生成A-*/
	int row = sizeof(_mat_A) / sizeof(_mat_A[0]);
	int column = sizeof(_mat_A[0]) / sizeof(_mat_A[0][0]);
	Matirx*  mat_A = Matrix_gen(row,column,_mat_A);
	/*mat_B*/
//	MATRIX_TYPE _mat_B [2][1] = {-2,-4}; //Two-dimensional function
	MATRIX_TYPE _mat_B [3][1] = {-2,-2,-2}; //Three-dimensional function 
	/*-具体生成B-*/
	row = sizeof(_mat_B) / sizeof(_mat_B[0]);
	column = sizeof(_mat_B[0]) / sizeof(_mat_B[0][0]);
	Matirx*  mat_B = Matrix_gen(row,column,_mat_B);
	/*XAX'+BX*/
	Matirx* mat_Result = M_mul(M_mul(mat_X,mat_A),M_T(mat_X));
	Matirx*  mat_temp = mat_Result;
	mat_Result = M_add_sub(1,mat_temp,-1,M_mul(mat_X,mat_B));
	M_free(mat_temp);
	/*XAX'+BX + C*/ 
	double Result = (mat_Result->data)[0] -2; 
	M_free(mat_Result);
	M_free(mat_A);
	M_free(mat_B);
	M_free(mat_X);
    return Result;
}

double Func_derivative_sub(double X,double step,Matirx* _mat_X, Matirx* _mat_gradient){/*F(x_k+lamda*d_k)
以lamda为变量的函数
*/ 
	double (*Func)(double,Matirx*,Matirx*) = Function_sub;
	return ((*Func)(X+step,_mat_X,_mat_gradient)-(*Func)(X,_mat_X,_mat_gradient))/step;
}

double Bisection_Method(Matirx* _mat_X, Matirx* _mat_gradient){/*Bisection_Method for Lamda
min[lamda] F(x_k+lamda*d_k)
二分法求解最优lamda值(X更新步长)*/
	/*Data + End gap*/
	double Left = -20;
	double Right = 20;
	double deta = 0.1;
	double step = 0.1; 
	/*Operation*/
	double (*Func)(double,Matirx*,Matirx*) = Function_sub;
	double X_left = Left, X_right = Right, X_mid;
	while((X_right-X_left)>deta){
		X_mid = (X_right+X_left)/2;
		double derivative = Func_derivative_sub(X_mid,step,_mat_X,_mat_gradient);
		if (derivative>0){
			X_right = X_mid;
		}else{
			if (derivative<0){
				X_left = X_mid;
			}else{
				X_left = X_mid - step;
				X_right = X_mid + step;				
			}
		}
	}
		double X = (X_right+X_left)/2;
		double Value = (*Func)(X,_mat_X, _mat_gradient);
		printf("X* = %lf, F(X*) = %lf",X,Value);
	return X;
}

