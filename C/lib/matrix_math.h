#include <stdio.h>
#include <malloc.h>
#include <math.h>
#define Pi 3.141592653


#define EYES {{1,0,0},{0,1,0},{0,0,1}}
#define ZEROS {{0,0,0},{0,0,0},{0,0,0}} 

#define SIZE 3 

#define memery (double(*)[SIZE])malloc(sizeof(double)*SIZE*SIZE)
#define vemery (double * )malloc(sizeof(double)*SIZE)

typedef double (*matrix)[SIZE];  
typedef double  *vector;  

matrix eyes()
{
	matrix Mtemp;Mtemp=memery; 
	int i,j;
	for (i=0;i<3;i++)
	{
		for(j=0;j<3;j++){
			if (i==j)
			  Mtemp[i][j]=1;
			else
			  Mtemp[i][j]=0;
		}
	}
	return Mtemp;
}

matrix zeros()
{
	matrix Mtemp;Mtemp=memery; 
	int i,j;
	for (i=0;i<3;i++)
	{
		for(j=0;j<3;j++){
			  Mtemp[i][j]=0;
		}
	}
	return Mtemp;
}


void mprint(matrix M)
{
  int i,j,k;
  
  for(i=0;i<SIZE;i++)
	{
		for(j=0;j<SIZE;j++)
		{
			printf("%lf ",M[i][j]);
		}
		printf("\n");
	}	
	printf("\n");
}
void vprint(vector v,int len)
{
	int i,j,k;
	for(i=0;i<len;i++)
	{
		printf("%lf ",v[i]);
	}
	printf("\n");
}
matrix mpm(matrix M1,matrix M2)
{
	int i,j,k;
	matrix Mtemp;Mtemp=zeros(); 
	
	for(i=0;i<SIZE;i++)
	{
		for(j=0;j<SIZE;j++)
		{
			Mtemp[i][j]=M1[i][j]+M2[i][j];
		}
	}
	return Mtemp;
 } 

matrix mm2m(matrix M1,matrix M2)
{
	int i,j,k;
	matrix Mtemp;Mtemp=zeros(); 

	for(i=0;i<SIZE;i++)
	{
		for(j=0;j<SIZE;j++)
		{
			for(k=0;k<SIZE;k++)
			{
				Mtemp[i][j]+=M1[i][k] * M2[k][j];
			}	
		}
	}	
    return Mtemp;
    free(Mtemp);
}

matrix vv2m(vector v1,vector v2)
{
	int i,j,k;
	matrix Mtemp;Mtemp=zeros(); 
	for(i=0;i<SIZE;i++)
	{
		for(j=0;j<SIZE;j++)
		{
			Mtemp[i][j]=v1[i]*v2[j];	
		}
	}	
    return Mtemp;
}

double vv2n(vector v1,vector v2)
{
	double out;
	out=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
	return out;
}

vector mv2v(matrix M,vector v)
{
	vector vec;vec=vemery;
	int i,j;
	for(i=0;i<SIZE;i++)
	{
		vec[i]=0;
		for(j=0;j<SIZE;j++)
		{
			vec[i]+=M[i][j]*v[j];	
		}
	}	
	
	return vec;
}

vector vm2v(vector v,matrix M)
{
	vector vec;vec=vemery;
	int i,j;
	for(i=0;i<SIZE;i++)
	{
		vec[i]=0;
		for(j=0;j<SIZE;j++)
		{
			vec[i]+=v[j]*M[j][i];	
		}
	}	
	
	return vec;
}

matrix nm2m(double num, matrix M)
{
	int i,j,k;
	matrix Mtemp;Mtemp=zeros(); 
	
	for(i=0;i<SIZE;i++)
	{
		for(j=0;j<SIZE;j++)
		{
			Mtemp[i][j]=num*M[i][j];
		}
	}
	return Mtemp;
}
matrix mn2m(matrix M,double num)
{
	int i,j,k;
	matrix Mtemp;Mtemp=zeros(); 
	
	for(i=0;i<SIZE;i++)
	{
		for(j=0;j<SIZE;j++)
		{
			Mtemp[i][j]=num*M[i][j];
		}
	}
	return Mtemp;
}

matrix trans(matrix M)
{
	int i,j,k;
	matrix Mtemp;Mtemp=zeros(); 
	
	for(i=0;i<SIZE;i++)
	{
		for(j=0;j<SIZE;j++)
		{
			Mtemp[i][j]=M[j][i];
		}
	}
	return Mtemp;
 } 
vector v_v(vector v1,vector v2)
{
	vector vec;vec=vemery;
	int i;
	for(i=0;i<3;i++)
	{
		vec[i]=v1[i]-v2[i];
	}
	return vec;
}
vector vpv(vector v1,vector v2)
{
	vector vec;vec=vemery;
	int i;
	for(i=0;i<3;i++)
	{
		vec[i]=v1[i]+v2[i];
	}
	return vec;
}

matrix inverse(matrix M)
{
	matrix mtemp;
	double out[3][3]={
	{(-M[1][2])*M[2][1] + M[1][1]*M[2][2],   M[0][2] *M[2][1] - M[0][1]*M[2][2],(-M[0][2])*M[1][1] + M[0][1]*M[1][2]},
	{  M[1][2] *M[2][0] - M[1][0]*M[2][2], (-M[0][2])*M[2][0] + M[0][0]*M[2][2],  M[0][2] *M[1][0] - M[0][0]*M[1][2]},
	{(-M[1][1])*M[2][0] + M[1][0]*M[2][1],   M[0][1] *M[2][0] - M[0][0]*M[2][1],(-M[0][1])*M[1][0] + M[0][0]*M[1][1]}
	};
	double temp;
	temp= -M[0][2]*M[1][1]*M[2][0] + M[0][1]*M[1][2]*M[2][0] + 
           M[0][2]*M[1][0]*M[2][1] - M[0][0]*M[1][2]*M[2][1] - 
           M[0][1]*M[1][0]*M[2][2] + M[0][0]*M[1][1]*M[2][2];
    if (temp==0)
      {
      	printf("matrix is singular\n");
      	temp=1;
	  }
    return mn2m(out,1/temp);
} 

matrix condition_inverse(matrix M,int num)
{
    double out[3][3],temp;
	if (num==1){	
	double out[3][3]={
	{(-M[1][2])*M[2][1] + M[1][1]*M[2][2],   M[0][2] *M[2][1] - M[0][1]*M[2][2],(-M[0][2])*M[1][1] + M[0][1]*M[1][2]},
	{  M[1][2] *M[2][0] - M[1][0]*M[2][2], (-M[0][2])*M[2][0] + M[0][0]*M[2][2],  M[0][2] *M[1][0] - M[0][0]*M[1][2]},
	{(-M[1][1])*M[2][0] + M[1][0]*M[2][1],   M[0][1] *M[2][0] - M[0][0]*M[2][1],(-M[0][1])*M[1][0] + M[0][0]*M[1][1]}
	};
	double temp;
	temp= -M[0][2]*M[1][1]*M[2][0] + M[0][1]*M[1][2]*M[2][0] + 
           M[0][2]*M[1][0]*M[2][1] - M[0][0]*M[1][2]*M[2][1] - 
           M[0][1]*M[1][0]*M[2][2] + M[0][0]*M[1][1]*M[2][2];
    if (temp==0)
      {
      	printf("matrix is singular\n");
      	mprint(M);
      	system("pause");
	  }
	}
	else{
	double out[3][3]={
	{M[1][1]*M[2][2],0,-M[0][2]*M[1][1]},
	{0,0,0},
	{-M[1][1]*M[2][0],0,M[0][0]*M[1][1]}
	};
	double temp;
	temp= -M[0][2]*M[2][0] + M[0][0]*M[2][2];
    if (temp==0)
      {
      	printf("matrix is singular\n");
      	mprint(M);
      	system("pause");
	  }
	}
    return mn2m(out,1/temp);
} 

//int self_condition_inverse(matrix M,int num)
//{
//    matrix mtemp;
//	double out[3][3],temp;
//	if (num==1){	
//	double out[3][3]={
//	{(-M[1][2])*M[2][1] + M[1][1]*M[2][2],   M[0][2] *M[2][1] - M[0][1]*M[2][2],(-M[0][2])*M[1][1] + M[0][1]*M[1][2]},
//	{  M[1][2] *M[2][0] - M[1][0]*M[2][2], (-M[0][2])*M[2][0] + M[0][0]*M[2][2],  M[0][2] *M[1][0] - M[0][0]*M[1][2]},
//	{(-M[1][1])*M[2][0] + M[1][0]*M[2][1],   M[0][1] *M[2][0] - M[0][0]*M[2][1],(-M[0][1])*M[1][0] + M[0][0]*M[1][1]}
//	};
//	temp= -M[0][2]*M[1][1]*M[2][0] + M[0][1]*M[1][2]*M[2][0] + 
//           M[0][2]*M[1][0]*M[2][1] - M[0][0]*M[1][2]*M[2][1] - 
//           M[0][1]*M[1][0]*M[2][2] + M[0][0]*M[1][1]*M[2][2];
//    if (temp==0){mark=1;}
//    else memcpy(M,mtemp=mn2m(out,1/temp),mem);free(mtemp);
//	}
//	else{
//	double out[3][3]={
//	{M[1][1]*M[2][2],0,-M[0][2]*M[1][1]},
//	{0,0,0},
//	{-M[1][1]*M[2][0],0,M[0][0]*M[1][1]}
//	};
//	temp= -M[0][2]*M[2][0] + M[0][0]*M[2][2];
//    if (temp==0){mark=1;}
//    else memcpy(M,mtemp=mn2m(out,1/temp),mem);free(mtemp);
//	}
//    return mark;
//} 
