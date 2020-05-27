#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include "matrix.h"
/*Grandient Method
多元函数[高/二次]-梯度算法 
Aouthor:Xipnig.Yu
Date:2020.05.27 
------------ 
min F(X)
s.t. X∈R^n
*/ 
double Function(Matirx* mat_X);
Matirx* Func_derivative(Matirx* mat_X, double step);

int main(int argc, char *argv[]) {
	/*Data + End gap*/
	double (*Func)(Matirx*) = Function;
//	MATRIX_TYPE _mat_X [1][2] = {2,1}; //Two-dimensional function
	MATRIX_TYPE _mat_X [1][3] = {2,1,1}; //Three-dimensional function 
	double deta = 0.01;
	double step_gradient = 0.01;
	double step_lamda = 0.01;
	/*Generate mat_X*/
	int row = sizeof(_mat_X) / sizeof(_mat_X[0]);
	int column = sizeof(_mat_X[0]) / sizeof(_mat_X[0][0]);
	Matirx*  mat_X = Matrix_gen(row,column,_mat_X);
	MATRIX_TYPE Mat_norm = 10*deta;
	Matirx* mat_grandient = NULL,* mat_temp = NULL;
	/*Operation*/
	while(Mat_norm>deta){
		mat_grandient = Func_derivative(mat_X,step_gradient);
		mat_temp = mat_X;
		mat_X = M_add_sub(1,mat_temp,step_lamda,mat_grandient);
		M_free(mat_temp);
		M_free(mat_grandient);
		Mat_norm = M_norm(mat_grandient,2);
	}
	printf("X* =\n");
	M_print(mat_X);
	printf("F(X*) = %lf",(*Func)(mat_X));
	return 0;
}

double Function(Matirx* mat_X){
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
    return Result;
}

Matirx* Func_derivative(Matirx* mat_X, double step){
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
