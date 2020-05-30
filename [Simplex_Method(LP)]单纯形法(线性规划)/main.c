#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include "Simplex_state.h"

/*Simplex method
线性规划-单纯形法 
Aouthor:Xipnig.Yu
Date:2020.05.29 
------------ 
min CX 
s.t. AX=b,X>=0
*/ 

int main(int argc, char *argv[]) {
	/*Data + End gap*/
	MATRIX_TYPE _mat_A [4][7]  = {
	1,1,1,1,0,0,0,
    1,0,0,0,1,0,0,
    0,0,1,0,0,1,0,
    0,3,1,0,0,0,1};
    MATRIX_TYPE _mat_B [4][1]  =
    {4,2,3,6}; 
    MATRIX_TYPE _mat_C [1][7]  = 
	{1,14,6,0,0,0,0};
	/*Generate Matrix*/
	/*M A*/
	int row = sizeof(_mat_A) / sizeof(_mat_A[0]);
	int column = sizeof(_mat_A[0]) / sizeof(_mat_A[0][0]);
	Matirx*  mat_A = Matrix_gen(row,column,_mat_A);
	/*M B*/
	row = sizeof(_mat_B) / sizeof(_mat_B[0]);
	column = sizeof(_mat_B[0]) / sizeof(_mat_B[0][0]);
	Matirx*  mat_B = Matrix_gen(row,column,_mat_B);
	/*M C*/
	row = sizeof(_mat_C) / sizeof(_mat_C[0]);
	column = sizeof(_mat_C[0]) / sizeof(_mat_C[0][0]);
	Matirx*  mat_C = Matrix_gen(row,column,_mat_C);
	/*调用单纯形法求解*/ 
	Matirx*  mat_test = Simplex_method(mat_A, mat_B,mat_C); 
	printf("$$最终所得解X*：\n");
	M_print(mat_test);
	return 0;
}

