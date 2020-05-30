#include "matrix.h"

#define PUTOUT_R_001	 	"<<ROW>>C_B   --  ��   --  b  --   x_1~x_n (PUTOUT_R)\n"
#define PUTOUT_C_001	 	"<<Column>> x_1~x_m   --  deta (PUTOUT_C)\n"

Matirx* Simplex_method(Matirx* mat_A,Matirx* mat_B,Matirx* mat_C);


Matirx* Simplex_method(Matirx* mat_A,Matirx* mat_B,Matirx* mat_C){
	double _sita_ = 0.001;
	Matirx* Result = M_Zeros(1,mat_A->column);//������ 
	
	Matirx* deta = M_Zeros(1,mat_A->column);//�������б� 
	Matirx* sita_arry = M_Zeros(1,mat_A->row);
	MATRIX_TYPE sita =0;
	Matirx* Coef = M_Zeros(1,mat_A->row);//ϵ��ֵ
	Matirx* Base = NULL;//M_Zeros(1,mat_A->row);//ѡ��Ļ�
	Matirx* Decsion  = Matrix_copy(mat_B);//����ȡֵ
	Matirx* END_deta = M_Ones(1,mat_A->column);//ʹ��detaֵ�궨�������
	Matirx* ZEROS_ceofA = M_Zeros(mat_A->row,1);//�����ж�����⣬�޽�
	Matirx* Position = NULL, * index_I, *mat_temp =NULL, *mat_temp2 =NULL, *mat_temp3 =NULL;
	Matirx* PRINT = NULL; //���ڴ�ӡ�����α� 
	int FLAG = 0;//��������ѭ����0-����//1-Ψһ��(11)or������(12)//2-�޽��//3-�޿��н�
	int Step = 1;
	int i,temp_1,temp_2,temp,Position_in,Position_out;
	MATRIX_TYPE temp_num, temp_deta, temp_val;
	
	while(1){ 
		if (mat_A->row == 1){
			/*ʵ�� Position(find(A == 1)) = 1;*/
			Matirx* temp_position = M_find(mat_A,1);
			Position = M_Zeros(1,(mat_A->column)*(mat_A->row));
			for (i=0;i<temp_position->row;i++){
				temp_1 = temp_position->data[i];
				Position->data[temp_1] = 1;
			}
			M_free(temp_position);
		}else{
			temp_1 = M_logic(M_add_sub(1,M_sum(mat_A),1,M_Ones(1,mat_A->column)),NULL,_NOT_);
			temp_2 = M_logic(M_minax_val(mat_A,M_min(mat_A)),NULL,_NOT_);
			Position = M_logic(temp_1,temp_2,_AND_);
		}
		index_I = M_Zeros(1,mat_A->row);
		for(i=0;i<mat_A->column;i++){
			if ((Position->data[i])!=0){
				M_print(M_Cut(mat_A,_HEAD_,_END_,i+1,i+1));
				temp = (M_max(M_Cut(mat_A,_HEAD_,_END_,i+1,i+1)))->data[0];
				index_I->data[temp] = i;
			}
		}
		M_print(index_I);
		
		Base = index_I; 
		M_print(mat_C);
		Coef =  M_setval(Coef,mat_C,index_I,_ORD4VAL_);
		M_print(Coef);
		
		//Ѱ�ң���������/detaֵ �� sitaֵ
    	//��һ������deta�����sita��������Ƿ����Ż��ռ䣿
    	/*deta*/
    	M_print(mat_A);
    	for (i=0;i<mat_A->column;i++){
    		mat_temp3 = M_Cut(mat_A,_HEAD_,_END_,i+1,i+1);
    		M_print(mat_temp3);
    		M_print(Coef);
    		mat_temp = M_mul(Coef,mat_temp3);
    		M_print(mat_temp);
			deta->data[i] = mat_C->data[i] - mat_temp->data[0];
			M_print(deta);
			M_free(mat_temp3);
			M_free(mat_temp);
		}
		Position_in = (M_max(deta))->data[0];
		temp_deta = deta->data[Position_in];
		/**************/
		/*��ӡ�����α�*/
		PRINT = M_Zeros(1+(mat_A->row)+1,3+(mat_A->column));
 		PRINT = M_matFull(PRINT,0,3,mat_C);
 		PRINT = M_matFull(PRINT,1,0,M_T(Coef));
 		PRINT = M_matFull(PRINT,1,1,M_T(Base));
		PRINT = M_matFull(PRINT,1,2,Decsion);
		PRINT = M_matFull(PRINT,1,3,mat_A);
		PRINT = M_matFull(PRINT,(PRINT->row)-1,3,deta);
		M_print(PRINT);
		/*��ʾ*/
		printf("$$���� = %d\n",Step);
		printf(PUTOUT_C_001);
	 	printf(PUTOUT_R_001);
	 	/*******/
	 	Step = Step +1;
	 	/*******/
	 	/*sita*/
	 	M_print(Decsion); 
	 	M_print(M_Cut(mat_A,_HEAD_,_END_,Position_in+1,Position_in+1));
	 	sita_arry = M_pmuldiv(Decsion,M_Cut(mat_A,_HEAD_,_END_,Position_in+1,Position_in+1),_DIV_);
	 	M_print(sita_arry);
	 	for(i=1;i<(sita_arry->row);i++){
	 		if((sita_arry->data[i])<=0){
	 			sita_arry->data[i] = _FULLY_BIG_;
			}
 		}
	 	/*������Ҫע�⣬����ѡ���˹�����*/
	 	M_print(sita_arry);
	 	Position_out = (M_min(sita_arry))->data[0];
		sita = sita_arry->data[Position_out];
		for(i=1;i<(sita_arry->row);i++){
			if((sita_arry->data[i])==sita){
				if(Base->data[Position_out] < Base->data[i]);
				Position_out = i;
			} 
		} 
	 	/*�޽��*/
	 	if (sita>_INFINITE_){
	 		FLAG = 2;
            break;
		}
		/*�������-Ψһ��or������*/
		temp_1 = (M_max(deta))->data[0];
		temp_val = deta->data[temp_1];
		M_print(deta);
		if (temp_val <= _sita_){
	 		FLAG = 1;
            break;
		}
		/**********/
		/*�г��ȱ任*/
		temp_num = mat_A->data[(Position_out)*(mat_A->column)+Position_in];
		Decsion->data[Position_out] = (Decsion->data[Position_out])/temp_num;
		M_print(Decsion);
		mat_A = M_matFull(mat_A,(Position_out),0,M_numul(M_Cut(mat_A,Position_out+1,Position_out+1,_HEAD_,_END_),1/temp_num));
		M_print(mat_A);
		for (i=0;i<mat_A->row;i++){
			if (i!=Position_out){
				Decsion->data[i] = Decsion->data[i] - (mat_A->data[(i*(mat_A->column)+Position_in)])*(Decsion->data[Position_out]);
				M_print(Decsion);
				mat_temp2 = M_add_sub(1,M_Cut(mat_A,i+1,i+1,_HEAD_,_END_),(mat_A->data[i*(mat_A->column)+Position_in]),M_Cut(mat_A,Position_out+1,Position_out+1,_HEAD_,_END_));
				M_print(mat_temp2);
				mat_A = M_matFull(mat_A,i,0,mat_temp2);
				M_print(mat_A);
			}
		}
	}
	/*STEP 2 �������*/
	/*��Ϊ��Ψһ���Ž⡢��������Ž⡢�޽�⡢�޿��н�*/
	Result = M_setval(Result,Coef,Base,_ORD4INI_);
	if (FLAG == 1){
		for(i=0;i<(mat_A->column);i++){
			if(((deta->data[i])==0)&&(sita>0)){
				FLAG = 12;   //�����
             	break;
			}
		}
		if(FLAG ==12){
			printf("�����Թ滮���⣬���ڡ�����ࡱ���Ž�.\n");
			printf("����һ�������Ž����Ϊ��(x1~xn)\n");
			M_print(Result);
			printf("����Ŀ�꺯��ֵΪ��\n");
			M_print(M_mul(mat_C,M_T(Result)));
		}else{
			printf("�����Թ滮���⣬����Ψһ���Ž�.\n");
			("���Ž����Ϊ��(x1~xn)\n");
			M_print(Result);
	 		printf("����Ŀ�꺯��ֵΪ��\n");
  	 		M_print(M_mul(Coef,M_T(M_Cut(Result,1,1,1,Coef->column))));
		}
	}else{
		if(FLAG == 2){
			printf("�����Թ滮���⣬����޽磨���޴�.\n");
		}else{
			printf("�����Թ滮���⣬�޿��н�.");
		}
	}
	return Result;
}
