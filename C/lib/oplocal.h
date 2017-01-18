#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "universe_math.h"


#define Pi 3.141592653
#define EYES {{1,0,0},{0,1,0},{0,0,1}}
#define ZEROS {{0,0,0},{0,0,0},{0,0,0}}

#define memery (double(*)[SIZE])malloc(sizeof(double)*SIZE*SIZE)
#define vemery (double * )malloc(sizeof(double)*SIZE)
#define gemery (double(*)[3])malloc(sizeof(double)*3*NUM)
#define gem sizeof(double)*3*NUM
#define mem sizeof(double)*3*3
#define vem sizeof(double)*3

typedef double (*matrix)[3];  
typedef double  *vector;  
typedef double (*GEO)[3];  


double Project(double x,int k,double sup,double inf)
{
	
	if(k==1)
	{
	if(x<inf) x=inf;
	if(x>sup) x=sup;
	}
	return x;	  
} 

double junp(double x,int k,double sup,double inf)
{
	double box;
	box=sup-inf;
	if(k==1)
	{
	while(x<inf) x+=box;
	while(x>sup) x-=box;
	}
	return x;	  
}
void gProject(GEO geo,int NUM,double sup,double inf)
{
	int i,j;
	for(i=0;i<NUM;i++)
	{
		for(j=0;j<3;j++)
		{
			geo[i][j]=junp(geo[i][j],j,sup,inf);
		}
	}
}
//=============================================================
vector local_findmin(GEO geo,int NUM,int ele,vector dir,double sup,double inf)
{
   double x2[21][3]={0},val[21]={0};
   double scale;
   int position,i,j,k,count=0;
   vector x1;x1=vemery;
   
   memcpy(x1,geo[ele],vem);
   scale=norm(x1);
   
   j=1;
   while(j<4)
   {
   	if (norm(dir)==0) break;
   	for(i=0;i<21;i++)
   	{
   	  for(k=0;k<3;k++){geo[ele][k]=x2[i][k]=Project(x1[k]+pow(0.1,j)*(i-10)*(dir[k]),k,sup,inf);}
	  val[i]=vvalue(geo,NUM);	  
    }
    
    memcpy(x1,x2[position=getminnum(val,21)],vem);

    if (position==0||position==20) j--;
	else j++; 
   }
    
	for(k=0;k<3;k++) geo[ele][k]=x1[k];
      
   return x1;
}
//=============================================================

vector local_findmin_Armijo(GEO geo,int NUM,int ele,vector dir,double sup,double inf)
{
   double x1[3],x2[3]={0},delta[3];
   int i=2,k;
   vector out;
   double sigma=0.0001,lambda;
   double val1,val2;
   for(k=0;k<3;k++){x1[k]=geo[ele][k];}
   val1=vvalue(geo,NUM);
    
   do{
   i--;
   lambda=pow(2,i); 
   for(k=0;k<3;k++)
	   {
	   geo[ele][k]=x2[k]=Project(x1[k]-lambda*dir[k],k,sup,inf);
	   delta[k]=x2[k]-x1[k];
	   }
   val2=vvalue(geo,NUM);
   
   if(i<-20)
       break; 
   }while(val2-val1<-sigma*pow(norm(delta),2)/lambda);
   
   lambda=pow(2,i+1); 
   for(k=0;k<3;k++){x2[k]=Project(x1[k]-lambda*dir[k],k,sup,inf);}  
   
	   	     
      
   out=vemery;memcpy(out,x2,sizeof(double)*3);
   return out;
}
//=============================================================

matrix matrixchange(matrix B,double rho,vector y,vector s)
{
	double left[3][3];
	double right[3][3];
	double extra[3][3];
	int i,j,k;
	
	matrix H;H=memery;
	matrix mtemp;
	
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			left[i][j]=-rho*y[i]*s[j];
			right[i][j]=-rho*s[i]*y[j];
			extra[i][j]=rho*s[i]*s[j];
			if(i==j)
			{
				left[i][j]++;
				right[i][j]++;
			}
		}
	}
    
    memcpy(H,mtemp=mm2m(left,B),mem);free(mtemp);
    memcpy(H,mtemp=mm2m(H,right),mem);free(mtemp);
    memcpy(H,mtemp=mpm(H,extra),mem);free(mtemp);
	
	return H;  
}


matrix matrixchange2(matrix B,double rho,vector y,vector s,matrix I)
{
	double left[3][3];
	double right[3][3];
	double extra[3][3];
	int i,j,k;
	
	matrix H;H=memery;
	matrix mtemp;
	
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			left[i][j]=-rho*y[i]*s[j];
			right[i][j]=-rho*s[i]*y[j];
			extra[i][j]=rho*s[i]*s[j];
			if(i==j)
			{
				left[i][j]++;
				right[i][j]++;
			}
		}
	}
    
    memcpy(H,mtemp=mm2m(left,I),mem);free(mtemp);
    memcpy(H,mtemp=mm2m(H,B),mem);free(mtemp);
    memcpy(H,mtemp=mm2m(H,I),mem);free(mtemp);
    memcpy(H,mtemp=mm2m(H,right),mem);free(mtemp);
    memcpy(H,mtemp=mpm(H,extra),mem);free(mtemp);
	
	return H; 
}
//=============================================================
matrix matrixchange_b_p(matrix Ac,vector y,vector s,matrix P)
{		
	matrix H;H=memery;
	matrix left;left=memery;
	matrix right;right=memery;
	matrix mid;mid=memery;
	matrix mtemp;
	vector vmid=vemery;
	vector vtemp;
	
	memcpy(left,mtemp=mm2m(P,Ac),mem);free(mtemp);
    memcpy(left,mtemp=mm2m(left,P),mem);free(mtemp);
    
    memcpy(mid,mtemp=vv2m(y,y),mem);free(mtemp);
    memcpy(mid,mtemp=mn2m(mid,vv2n(y,s)),mem);free(mtemp);
    
    memcpy(vmid,vtemp=mv2v(Ac,s),vem);free(vtemp);
    
    memcpy(right,mtemp=vv2m(vmid,vmid),mem);free(mtemp);
    memcpy(right,mtemp=mn2m(right,-vv2n(s,vmid)),mem);free(mtemp);
    memcpy(right,mtemp=mm2m(P,right),mem);free(mtemp);
    memcpy(right,mtemp=mm2m(right,P),mem);free(mtemp);
    

    memcpy(H,mtemp=mpm(left,mid),mem);free(mtemp);
	memcpy(H,mtemp=mpm(H,right),mem);free(mtemp);

    free(left);free(right);free(mid);
    
	return H; 
}
//=============================================================
void op_local(GEO geo,int NUM)
{
	int nobest=0;
	int i,j,k;
	double vlast=0,vmin=0;
	int ele;
	double val1,val2,sup,inf;

	double rho,temp1,temp2;
	
	vector g1=vemery,g2=vemery,x1=vemery,x2=vemery;	
	vector s=vemery,y=vemery;
	vector dir=vemery;
	matrix B;B=memery;
	int *xu;xu=(int * )malloc(sizeof(int)*NUM);
	int *xtemp;
	vector vtemp;
	matrix mtemp;
	
	while(nobest<(NUM/2))
	{  
	    
       memcpy(xu,xtemp=randperm(NUM),sizeof(int)*NUM);free(xtemp);
       
	   for(i=0;i<NUM;i++)
	   {
	   	  ele=xu[i];
	   	  memcpy(g1,vtemp=df(geo,NUM,ele),vem);free(vtemp);
	   	  val1=vvalue(geo,NUM);
	   	  val2=val1;
	   	  memcpy(x1,geo[ele],vem);
	   	  memcpy(B,mtemp=eyes(),mem);free(mtemp);
	   	  do
	   	 {	
	   	     val1=val2;
	   	  	 memcpy(dir,vtemp=mv2v(B,g1),vem);free(vtemp);
	   	  	 memcpy(x2,vtemp=local_findmin(geo,NUM,ele,dir,sup,inf),vem);free(vtemp);
			 val2=vvalue(geo,NUM);
			 memcpy(g2,vtemp=df(geo,NUM,ele),vem);free(vtemp);
			 memcpy(s,vtemp=v_v(x2,x1),vem);free(vtemp);
			 memcpy(y,vtemp=v_v(g2,g1),vem);free(vtemp);
			 if (norm(s)<0.001)
			     break;
			 rho=vv2n(s,y);rho=1/rho;
			 memcpy(B,mtemp=matrixchange(B,rho,y,s),mem);free(mtemp);
			 memcpy(g1,g2,vem);
			 memcpy(x1,x2,vem);
	     }while (val1-val2>0.001);
		  vmin=val2;
		  
	   }
	   if (vlast-vmin<0.001)
	      nobest=nobest+1;
	   vlast=vmin;
	}
    free(g1);free(g2);free(x1);free(x2);free(dir);free(s);free(y);free(B);free(xu);
} 
//=====================================TEMP=====================================

double op_local_b(GEO geo,int NUM,double sup,double inf)
{
	int nobest=0;
	int i,j,k,breaknum,breaknum2=0;
	double vlast=0,vmin=0;
	int ele;
	double val1,val2;

	double rho,temp1,temp2;
	
	vector g1=vemery,g2=vemery,x1=vemery,x2=vemery;	
	vector s=vemery,y=vemery;
	vector dir=vemery;
	matrix B=memery;
	double pg[3],epsilon;
	int *xu;xu=(int * )malloc(sizeof(int)*NUM);
	int *xtemp;
	vector vtemp;
	matrix mtemp;
	matrix H=memery,A=memery,I=memery;
  
	val1=vvalue(geo,NUM);
	breaknum2=0;
	while(nobest<3)
	{  
	    
       memcpy(xu,xtemp=randperm(NUM),sizeof(int)*NUM);free(xtemp);
       
	   for(i=0;i<NUM;i++)
	   {
	   	
	   	  ele=xu[i];
	   	   
	   	   
	   	  val1=vvalue(geo,NUM);
	   	  memcpy(g1,vtemp=df(geo,NUM,ele),vem);free(vtemp); 
			 
	   	  val2=val1;
	   	  
	   	   memcpy(x1,geo[ele],vem);
	   	   memcpy(I,mtemp=eyes(),mem);free(mtemp);
		   memcpy(A,mtemp=zeros(),mem);free(mtemp);
	   	   
	       if (x1[1]>sup||x1[1]<inf){I[1][1]=0;A[1][1]=1;}
	   	   memcpy(B,I,mem);
	   	   
	   	   breaknum=0;
	   	  do
	   	 {				
			 val1=val2;
	   	     
	   	      for(k=0;k<3;k++){pg[k]=x1[k]-Project(x1[k]-g1[k],k,sup,inf);}
	   	      epsilon=norm(pg);if (epsilon>(sup-inf)/2) epsilon=(sup-inf)/2;
	   	      
	   	         
	   	  	 memcpy(H,mtemp=mpm(A,B),mem);free(mtemp);

			 memcpy(dir,vtemp=mv2v(H,g1),vem);free(vtemp);
			  
	   	  	 memcpy(x2,vtemp=local_findmin(geo,NUM,ele,dir,sup,inf),vem);free(vtemp);
	   	  	 
			 val2=vvalue(geo,NUM);
			 
			  memcpy(g2,vtemp=df(geo,NUM,ele),vem);free(vtemp);
			  
			  memcpy(I,mtemp=eyes(),mem);free(mtemp);
	   	      memcpy(A,mtemp=zeros(),mem);free(mtemp);
	   	      
			  if (x2[1]>sup-epsilon||x2[1]<inf+epsilon){I[1][1]=0;A[1][1]=1;}
			  
			 memcpy(s,vtemp=v_v(x2,x1),vem);free(vtemp);
			 memcpy(y,vtemp=v_v(g2,g1),vem);free(vtemp);
			 
			 memcpy(s,vtemp=mv2v(I,s),vem);free(vtemp);	 
			 memcpy(y,vtemp=mv2v(I,y),vem);free(vtemp);
			 
//			system("pause");   
			rho=vv2n(s,y);
			
//			rho=1/rho;
//			memcpy(B,mtemp=matrixchange2(B,rho,y,s,I),mem);free(mtemp);
			 
			 if (rho>0)
			    {
			     rho=1/rho;
			     memcpy(B,mtemp=matrixchange2(B,rho,y,s,I),mem);free(mtemp);
			     
	         	}
			 else 
			    {
				  memcpy(B,I,mem);
				}
			
			 memcpy(g1,g2,vem);
			 memcpy(x1,x2,vem);
//			 printf("epsilon=%lf\n",epsilon);
//			 printf("val2=%lf\n",val2);
             breaknum++;
             if (breaknum>1000)
                 {
//                 	 printf("bingo\n"); 
					 break;
				 }
	     }while (epsilon>0.01);
		  vmin=val2;
	   }
	   if (vlast-vmin<0.001)
	      nobest=nobest+1;
	   breaknum2++;
	      
	   vlast=vmin;
	   if (breaknum2>10)
	     {
			break;
		}
	}
	
    free(g1);free(g2);free(x1);free(x2);free(dir);free(s);free(y);free(B);free(xu);free(A);free(I);
    
    return vlast;
    
} 
//=====================================TEMP=====================================
double op_local_b_f_para(GEO geo,int NUM,double sup,double inf,double mineps,int maxloop1,int maxloop2)
{
	int nobest=0;
	int i,j,k,breaknum,breaknum2=0;
	int ele;
	int *xu;xu=(int * )malloc(sizeof(int)*NUM);
	int *xtemp;
	
	double val1,val2;
	double rho,temp1,temp2;
	double vlast=0,vmin=0;
	double pg[3],epsilon;
	
	vector g1=vemery,g2=vemery,x1=vemery,x2=vemery;	
	vector s=vemery,y=vemery;
	vector dir=vemery;
	vector vtemp;
	
	matrix B=memery;
	matrix mtemp;
	matrix H=memery,A=memery,I=memery;

	
	val1=vvalue(geo,NUM);val2=val1;
 	memcpy(I,mtemp=eyes(),mem);free(mtemp);
	memcpy(A,mtemp=zeros(),mem);free(mtemp);     
    memcpy(xu,xtemp=randperm(NUM),sizeof(int)*NUM);free(xtemp);
    
	while(nobest<2)
	{     
	   for(i=0;i<NUM;i++)
	   {	
	   	   ele=xu[i];	  
			   	   
	   	   memcpy(g1,vtemp=df(geo,NUM,ele),vem);free(vtemp);   
	   	   memcpy(x1,geo[ele],vem);
	   	    
	       if (x1[1]>sup||x1[1]<inf){I[1][1]=0;A[1][1]=1;}
	       else{I[1][1]=1;A[1][1]=0;}
	   	   memcpy(B,I,mem);
	   	   
	   	   breaknum=0;
	   	  do
	   	 {		
	   	      for(k=0;k<3;k++){pg[k]=x1[k]-Project(x1[k]-g1[k],k,sup,inf);}
	   	      epsilon=norm(pg);if (epsilon>(sup-inf)/2) epsilon=(sup-inf)/2;
	   	      
	   	  	 memcpy(H,mtemp=mpm(A,B),mem);free(mtemp);
             memcpy(dir,vtemp=mv2v(H,g1),vem);free(vtemp);		  
	   	  	 memcpy(x2,vtemp=local_findmin(geo,NUM,ele,dir,sup,inf),vem);free(vtemp);
	   	  	 
			 val2=vvalue(geo,NUM);
			 
			 memcpy(g2,vtemp=df(geo,NUM,ele),vem);free(vtemp);
			  
	   	     if (x2[1]>sup-epsilon||x2[1]<inf+epsilon){I[1][1]=0;A[1][1]=1;}
			 else{I[1][1]=1;A[1][1]=0;}
			  
			 memcpy(s,vtemp=v_v(x2,x1),vem);free(vtemp);memcpy(s,vtemp=mv2v(I,s),vem);free(vtemp);
			 memcpy(y,vtemp=v_v(g2,g1),vem);free(vtemp);memcpy(y,vtemp=mv2v(I,y),vem);free(vtemp);
 
			 rho=vv2n(s,y);
			 if (rho>0) {rho=1/rho;	memcpy(B,mtemp=matrixchange2(B,rho,y,s,I),mem);free(mtemp);}
			 else {memcpy(B,I,mem);}
			
			 memcpy(g1,g2,vem);
			 memcpy(x1,x2,vem);
             val1=val2;
             
             breaknum++;if (breaknum>maxloop1) break;
                   
         }while (epsilon>mineps);
		     if(val2<=vmin)
			    vmin=val2;
	   }  
	   
	   if (vlast-vmin<0.001) nobest=nobest+1;
	   breaknum2++;
	   vlast=vmin;
	   if (breaknum2>maxloop2) break;
	}  

    free(g1);free(g2);free(x1);free(x2);free(dir);free(s);free(y);free(B);free(xu);free(A);free(I);
    
    return vlast;
    
} 
double op_local_b_f(GEO geo,int NUM,double sup,double inf)
{
	int nobest=0;
	int i,j,k,breaknum,breaknum2=0;
	int ele;
	int *xu;xu=(int * )malloc(sizeof(int)*NUM);
	int *xtemp;
	
	double val1,val2;
	double rho,temp1,temp2;
	double vlast=0,vmin=0;
	double pg[3],epsilon;
	
	vector g1=vemery,g2=vemery,x1=vemery,x2=vemery;	
	vector s=vemery,y=vemery;
	vector dir=vemery;
	vector vtemp;
	
	matrix B=memery;
	matrix mtemp;
	matrix H=memery,A=memery,I=memery;

	
	val1=vvalue(geo,NUM);val2=val1;
 	memcpy(I,mtemp=eyes(),mem);free(mtemp);
	memcpy(A,mtemp=zeros(),mem);free(mtemp);     
    memcpy(xu,xtemp=randperm(NUM),sizeof(int)*NUM);free(xtemp);
    
	while(nobest<2)
	{     
	   for(i=0;i<NUM;i++)
	   {	
	   	   ele=xu[i];	  
			   	   
	   	   memcpy(g1,vtemp=df(geo,NUM,ele),vem);free(vtemp);   
	   	   memcpy(x1,geo[ele],vem);
	   	    
	       if (x1[1]>sup||x1[1]<inf){I[1][1]=0;A[1][1]=1;}
	       else{I[1][1]=1;A[1][1]=0;}
	   	   memcpy(B,I,mem);
	   	   
	   	   breaknum=0;
	   	  do
	   	 {		
	   	      for(k=0;k<3;k++){pg[k]=x1[k]-Project(x1[k]-g1[k],k,sup,inf);}
	   	      epsilon=norm(pg);if (epsilon>(sup-inf)/2) epsilon=(sup-inf)/2;
	   	      
	   	  	 memcpy(H,mtemp=mpm(A,B),mem);free(mtemp);
             memcpy(dir,vtemp=mv2v(H,g1),vem);free(vtemp);		  
	   	  	 memcpy(x2,vtemp=local_findmin(geo,NUM,ele,dir,sup,inf),vem);free(vtemp);
	   	  	 
			 val2=vvalue(geo,NUM);
			 
			 memcpy(g2,vtemp=df(geo,NUM,ele),vem);free(vtemp);
			  
	   	     if (x2[1]>sup-epsilon||x2[1]<inf+epsilon){I[1][1]=0;A[1][1]=1;}
			 else{I[1][1]=1;A[1][1]=0;}
			  
			 memcpy(s,vtemp=v_v(x2,x1),vem);free(vtemp);memcpy(s,vtemp=mv2v(I,s),vem);free(vtemp);
			 memcpy(y,vtemp=v_v(g2,g1),vem);free(vtemp);memcpy(y,vtemp=mv2v(I,y),vem);free(vtemp);
 
			 rho=vv2n(s,y);
			 if (rho>0) {rho=1/rho;	memcpy(B,mtemp=matrixchange2(B,rho,y,s,I),mem);free(mtemp);}
			 else {memcpy(B,I,mem);}
			
			 memcpy(g1,g2,vem);
			 memcpy(x1,x2,vem);
             val1=val2;
  
             breaknum++;if (breaknum>100) break;
                   
         }while (epsilon>0.01);
		     if(val2<=vmin)
			    vmin=val2;
	   }  
	   
	   if (vlast-vmin<0.01) nobest=nobest+1;
	   breaknum2++;
	   vlast=vmin;
	   if (breaknum2>10) break;
	}  

    free(g1);free(g2);free(x1);free(x2);free(dir);free(s);free(y);free(B);free(xu);free(A);free(I);
    
    return vlast;
    
} 

double op_local_b_veryfast(GEO geoo,int NUM,double sup,double inf)
{
	GEO geo;geo=gemery;

	memcpy(geo,geoo,gem);
	
	int nobest=0;
	int i,j,k,breaknum,breaknum2;
	int ele;
	int *xu;xu=(int * )malloc(sizeof(int)*NUM);
	int *xtemp;
	
	double val1,val2;
	double rho,temp1,temp2;
	double vlast=0,vmin=0;
	double pg[3],epsilon;
	
	vector g1=vemery,g2=vemery,x1=vemery,x2=vemery;	
	vector s=vemery,y=vemery;
	vector dir=vemery;
	vector vtemp;
	
	matrix B=memery;
	matrix mtemp;
	matrix H=memery,A=memery,I=memery;

	
	val1=vvalue(geo,NUM);val2=val1;
 	memcpy(I,mtemp=eyes(),mem);free(mtemp);
	memcpy(A,mtemp=zeros(),mem);free(mtemp);     
    memcpy(xu,xtemp=randperm(NUM),sizeof(int)*NUM);free(xtemp);
    
	while(nobest<2)
	{     
	   for(i=0;i<NUM;i++)
	   {	
	   	   ele=xu[i];	  
			   	   
	   	   memcpy(g1,vtemp=df(geo,NUM,ele),vem);free(vtemp);   
	   	   memcpy(x1,geo[ele],vem);
	   	    
	       if (x1[1]>sup||x1[1]<inf){I[1][1]=0;A[1][1]=1;}
	       else{I[1][1]=1;A[1][1]=0;}
	   	   memcpy(B,I,mem);
	   	   
	   	   breaknum=0;
	   	  do
	   	 {		
	   	      for(k=0;k<3;k++){pg[k]=x1[k]-Project(x1[k]-g1[k],k,sup,inf);}
	   	      epsilon=norm(pg);if (epsilon>(sup-inf)/2) epsilon=(sup-inf)/2;
	   	      
	   	  	 memcpy(H,mtemp=mpm(A,B),mem);free(mtemp);
             memcpy(dir,vtemp=mv2v(H,g1),vem);free(vtemp);		  
	   	  	 memcpy(x2,vtemp=local_findmin(geo,NUM,ele,dir,sup,inf),vem);free(vtemp);
	   	  	 
			 val2=vvalue(geo,NUM);
			 
			 memcpy(g2,vtemp=df(geo,NUM,ele),vem);free(vtemp);
			  
	   	     if (x2[1]>sup-epsilon||x2[1]<inf+epsilon){I[1][1]=0;A[1][1]=1;}
			 else{I[1][1]=1;A[1][1]=0;}
			  
			 memcpy(s,vtemp=v_v(x2,x1),vem);free(vtemp);memcpy(s,vtemp=mv2v(I,s),vem);free(vtemp);
			 memcpy(y,vtemp=v_v(g2,g1),vem);free(vtemp);memcpy(y,vtemp=mv2v(I,y),vem);free(vtemp);
 
			 rho=vv2n(s,y);
			 if (rho>0) {rho=1/rho;	memcpy(B,mtemp=matrixchange2(B,rho,y,s,I),mem);free(mtemp);}
			 else {memcpy(B,I,mem);}
			
			 memcpy(g1,g2,vem);
			 memcpy(x1,x2,vem);
             val1=val2;
             
             breaknum++;if (breaknum>100) break;
                   
         }while (epsilon>0.1);
		     if(val2<=vmin)
			    vmin=val2;
	   }  
	   
	   if (vlast-vmin<0.001) nobest=nobest+1;
	   breaknum2++;
	   vlast=vmin;
	   if (breaknum2>10) break;
	}  

    free(g1);free(g2);free(x1);free(x2);free(dir);free(s);free(y);free(B);free(xu);free(A);free(I);
    free(geo);
    return vlast;
    
} 
//=====================================TEMP=====================================

//=====================================L-BFGS===================================
void cellplus(double store[][3],double s[],int m)
{
	int i,j;
	for(i=0;i<m-1;i++)
	{
		for(j=0;j<3;j++)
		  store[i][j]=store[i+1][j];
	}
	for(j=0;j<3;j++)
		  store[m-1][j]=s[j];
}

//=============================================================

matrix op_lbfgs(double s[][3],double y[][3],int m)
{
	double rho[m];
	double temp;
	double V[m][3][3];
	double r_k;
	 
	matrix VM;VM=eyes();
	vector VS;	
	matrix H;H=zeros();
	
	int i,j,k;
	
	for(i=0;i<m;i++)
	{
		temp=y[i][0]*s[i][0]+y[i][1]*s[i][1]+y[i][2]*s[i][2];
		temp=1/temp;
		rho[i]=temp;
		for(j=0;j<3;j++)
		{
			for(k=0;k<3;k++)
			{
				V[i][j][k]=-temp*y[i][j]*s[i][k];
				if(j==k)
				{
					V[i][j][k]++;
				}
			}
        }
	}
	
	for(i=m-1;i>=0;i--)
	{
			free(VS);VS=vemery;
		VS=mv2v(VM,s[i]); 
		H=mpm(H,nm2m(rho[i],vv2m(VS,VS)));
		VM=mm2m(VM,V[i]);
	}
	r_k=vv2n(s[m-1],y[m-1])/vv2n(y[m-1],y[m-1]);
	H=mpm(H,nm2m(r_k,mm2m(VM,trans(VM))));
	
	free(VS);VS=NULL;
	free(VM);VM=NULL;
    return H;
} 

GEO WL_randmove1(GEO geo,double r,int NUM)
{
	GEO temp;temp=gemery;
	memcpy(temp,geo,gem);
	int i;
	for(i=0;i<NUM;i++)
	{
		shift(temp[i],r);
	}
	return temp; 
}
//=============================================================
GEO WL_randmove1_b(GEO geo,double r1,double r2,double r3,int NUM,double sup,double inf)
{
	GEO temp;temp=gemery;
	memcpy(temp,geo,gem);
	int i;
	for(i=0;i<NUM;i++)
	{
		temp[i][0]+=r1*(2*random(1)-1);
		if ((temp[i][1]<sup-r2)&&(temp[i][1]>inf+r2))
		    temp[i][1]+=r2*(2*random(1)-1);
		else if (temp[i][1]>sup-r2)
		    temp[i][1]=sup-r2+r2*(2*random(1)-1); 
		else
		    temp[i][1]=inf+r2+r2*(2*random(1)-1); 
		temp[i][2]+=r3*(2*random(1)-1);

	}
	return temp; 
}
//=============================================================
GEO WL_randmove2(GEO geo,double r,int NUM)
{
	GEO temp;temp=gemery;
	memcpy(temp,geo,gem);
	
	int ele;
	int i;
	ele=random(NUM);if (ele==NUM) ele=random(NUM);
	shift(temp[ele],r);
	return temp; 
}

//=============================================================
int WL_flatness(vector H,int len,double para)
{
	double ave;double temp;
	int out=0,i;
	double M[len];
	if ((temp=sum(H,len))>len*100)
	   {
	   	ave=temp/len;
	   	for (i=0;i<len;i++){M[i]=H[i]/ave;}
	   	if((getmin(M,len)>para)&&(getmax(M,len)<(2-para))){out=1;}
	   }
	return out;
} 


