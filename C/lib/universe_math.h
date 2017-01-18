#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <string> 
#include "matrix_math.h"
#define Pi 3.141592653

#define EYES {{1,0,0},{0,1,0},{0,0,1}}
#define ZEROS {{0,0,0},{0,0,0},{0,0,0}} 

#define SIZE 3 

#define memery (double(*)[SIZE])malloc(sizeof(double)*SIZE*SIZE)
#define vemery (double * )malloc(sizeof(double)*SIZE)
#define fullvemery (double * )malloc(sizeof(double)*SIZE*NUM)

typedef double (*matrix)[3];  
typedef double  *vector;  
typedef double (*GEO)[3];  

void gprint(GEO geo,int NUM)
{
	int i,j;
	
	for (i=0;i<NUM;i++)
	{
		for(j=0;j<3;j++)
		{
			printf("%lf ",geo[i][j]);
		}
		printf("%\n");
	}
	printf("%\n");
}

double random(double lim)
{
	double out;
    out=rand()/double(RAND_MAX);
    out*=lim;
    return(out);
}
double value(GEO geo,int r_geo)
{   
	int c_geo = 3;
	
	int row_num=r_geo*(r_geo-1)/2;	
	
	double dis[row_num],temp,scale;
	
	double a=0,b=0,out;
	int i,j,t,k=0;
    for(i=0;i<r_geo;i++)
    {
    	for(j=0;j<i;j++)
    	{
    		if(i!=j)
    		{   
    		    temp=0;
    		    for(t=0;t<c_geo;t++){temp+=pow((geo[i][t]-geo[j][t]),2);}
    			temp=sqrt(temp);
    	    	a=a+pow(temp,(-6));
    	    	b=b+pow(temp,(-12));
			} 
		}
	}
	out=-pow(a,2)/b;
	scale=pow(b/a,1.0/6);
	for(i=0;i<r_geo;i++){
		for(j=0;j<3;j++)
			geo[i][j]*=scale;	   
	}
	return(out);
} 

double vvalue(GEO geo,int r_geo)
{   

	int c_geo = 3;
	
	int row_num=r_geo*(r_geo-1)/2;	
	
	double dis[row_num],temp,scale;
	
	double a=0,b=0,out;
	int i,j,t,k=0;
    for(i=0;i<r_geo;i++)
    {
    	for(j=0;j<i;j++)
    	{
    		if(i!=j)
    		{   
    		    temp=0;
    		    for(t=0;t<c_geo;t++){temp=temp+pow((geo[i][t]-geo[j][t]),2);}
    			temp=sqrt(temp);
    	    	a=a+pow(temp,(-6));
    	    	b=b+pow(temp,(-12));
			} 
		}
	}
	out=b-2*a;
	
	return(out);
} 


void shift(vector line,double radius)
{   

	double theta,phi;
	
	theta=Pi*random(1);
	phi=2*Pi*random(1);

	line[0]+=radius*cos(theta)*cos(phi);
    line[1]+=radius*cos(theta)*sin(phi);
    line[2]+=radius*sin(theta);
} 

void fileprint(GEO geo,int r_geo,char name[])
{
	FILE *fpt;
	char filename[]="D:\\";
	char filename2[]=".txt"; 
	int i;
	strcat(filename,name);
	strcat(filename,filename2);
	fpt=fopen(filename, "w");     
    for (i=0;i<r_geo;i++)
    {
    	fprintf(fpt,"%lf %lf %lf\n",geo[i][0],geo[i][1],geo[i][2]);
	}
	fprintf(fpt,"the min energy is %lf\n",vvalue(geo,r_geo)); 
    fclose(fpt); 
}


GEO initialization(int NUM,int STR)
{
    int i,j,nn=0,k;
    int basis_L,last_N,outside_N;
    double center_L;
    double zero[NUM][3]; 
    GEO geo;geo=(double(*)[3])malloc(sizeof(double)*3*NUM);
	if(STR==1)
    {
    	for(i=0;i<NUM;i++)
    	{
    		geo[i][0]=geo[i][1]=geo[i][2]=0;
			shift(geo[i],1);
		}
	}
    else if (STR==2)
	{
		if(NUM<27)
		   printf("Num is too small, suggest to use sphere\n");
		else
		{
			basis_L= floor(pow(NUM,1.0/3));
			outside_N=floor((NUM-pow(basis_L,3))/6);
			last_N=NUM-pow(basis_L,3)-outside_N*6; 
			for(i=0;i<basis_L;i++){
			   	for(j=0;j<basis_L;j++){
			   		for(k=0;k<basis_L;k++){
			   			geo[nn][0]=i;
			   			geo[nn][1]=j;
			   			geo[nn][2]=k;
			   			nn++;
					   }
				   }
			   }
		    
			center_L=(basis_L+1)/2; 
			if(outside_N!=0)
			{
				for(i=0;i<outside_N;i++)
				{
					geo[nn][0]=2*center_L;                      geo[nn][1]=center_L+basis_L*(2*random(1)-1);geo[nn][2]=center_L+basis_L*(2*random(1)-1);nn++;
					geo[nn][0]=-1;                              geo[nn][1]=center_L+basis_L*(2*random(1)-1);geo[nn][2]=center_L+basis_L*(2*random(1)-1);nn++;
					geo[nn][0]=center_L+basis_L*(2*random(1)-1);geo[nn][1]=2*center_L;                      geo[nn][2]=center_L+basis_L*(2*random(1)-1);nn++;
					geo[nn][0]=center_L+basis_L*(2*random(1)-1);geo[nn][1]=-1;                              geo[nn][2]=center_L+basis_L*(2*random(1)-1);nn++;
					geo[nn][0]=center_L+basis_L*(2*random(1)-1);geo[nn][1]=center_L+basis_L*(2*random(1)-1);geo[nn][2]=2*center_L;                      nn++;
					geo[nn][0]=center_L+basis_L*(2*random(1)-1);geo[nn][1]=center_L+basis_L*(2*random(1)-1);geo[nn][2]=-1;                              nn++;
				}
		    }      
		    
		    if(last_N!=0)
		    {
		    	for(i=0;i<last_N;i++)
		    	{
		    		geo[nn][0]=center_L;geo[nn][1]=center_L;geo[nn][2]=center_L;
					shift(geo[nn],center_L);
					nn++;
				}
			}
			
		}
	} 
	else
	  printf("please choose the right para\n");	
	return geo;	
}

GEO initialization_b(int NUM,double sup,double inf)
{
    int i,j,nn=0,k;
    int basis_L,last_N,outside_N;
    double center_L;
    double zero[NUM][3]; 
    double box=sup-inf;
    GEO geo;
	geo=(double(*)[3])malloc(sizeof(double)*3*NUM);
    	for(i=0;i<NUM;i++)
    	{
    		geo[i][0]=2*box*(2*random(1)-1);
    		geo[i][1]=box*random(1)+inf;
    		geo[i][2]=2*box*(2*random(1)-1);
		}
	
	
	
	return geo;	
}
int *randperm(int NUM)
{   
    int i,j,k,p,*xu;
    xu=(int * )malloc(sizeof(int)*NUM);
    for(i=0;i<NUM;i++)
    {
    	xu[i]=i;
	}
	
	for(p=0;p<NUM;p++)
	{
	for(i=0;i<NUM;i++)
	{
		
		j=rand()%NUM;
		k=xu[i];
		xu[i]=xu[j];
		xu[j]=k;
	}
    }
	return xu;
} 


vector df(GEO geo,int NUM,int k)
{
	vector ddf;ddf=vemery;
	int i,j;
	double numx,numy,numz;
	double a=0,b=0,c[3]={0},e[3]={0},temp,temp1,temp2;

	for (i=0;i<NUM;i++)
	{
		
		if(i!=k)
		{
				numx=geo[i][0]-geo[k][0];numy=geo[i][1]-geo[k][1];numz=geo[i][2]-geo[k][2];
				temp=pow(numx,2)+pow(numy,2)+pow(numz,2);
				temp1=pow(temp,-4);temp2=pow(temp,-7);
				c[0]+=temp1*numx;c[1]+=temp1*numy;c[2]+=temp1*numz;
				e[0]+=temp2*numx;e[1]+=temp2*numy;e[2]+=temp2*numz;	
		}
		
	}

	ddf[0]=12*(c[0]-e[0]);
	ddf[1]=12*(c[1]-e[1]);
	ddf[2]=12*(c[2]-e[2]);
	return ddf;
}
vector df2(GEO geo,int NUM)
{
	vector ddf;ddf=fullvemery;
	int i,j,k;
	double numx,numy,numz;
	double a=0,b=0,c[3]={0},e[3]={0},temp,temp1,temp2;

	for (i=0;i<NUM;i++)
	{   
	double c[3]={0},e[3]={0};
		for (k=0;k<NUM;k++)
		 {
		   if(i!=k)
		  {
				numx=geo[k][0]-geo[i][0];numy=geo[k][1]-geo[i][1];numz=geo[k][2]-geo[i][2];
				temp=pow(numx,2)+pow(numy,2)+pow(numz,2);
				temp1=pow(temp,-4);temp2=pow(temp,-7);
				c[0]+=temp1*numx;c[1]+=temp1*numy;c[2]+=temp1*numz;
				e[0]+=temp2*numx;e[1]+=temp2*numy;e[2]+=temp2*numz;	
		  }
		}
	

	ddf[3*i+0]=12*(c[0]-e[0]);
	ddf[3*i+1]=12*(c[1]-e[1]);
	ddf[3*i+2]=12*(c[2]-e[2]);
	}
	return ddf;
}

double norm(double *dir)
{
	double temp=0;int i;
	for(i=0;i<3;i++)
	{
		temp+=dir[i]*dir[i];
	}
	return sqrt(temp);
}
vector normlize(vector dir)
{
	vector out;out=vemery;
	int i;
	double temp;
	temp=norm(dir);
	out[0]=dir[0]/temp;
	out[1]=dir[1]/temp;
	out[2]=dir[2]/temp;
	return out;	
}
int getminnum(double val[],int len)
{
	int minnum,i;
	double min;
	min=val[0];
	minnum=0;
	for(i=0;i<len;i++)
	{
		if (val[i]<=min)
          {
          	min=val[i];
          	minnum=i;
		  }
	}
	if((minnum==len-1)&&(val[minnum]==val[minnum-1]))
	  minnum--;
	  
	return minnum;
}


double getmin(double val[],int len)
{
	int minnum,i;
	double min;
	min=val[0];
	minnum=0;
	for(i=1;i<len;i++)
	{
		if (val[i]<min)
          {
          	min=val[i];
          	minnum=i;
		  }
	}
	return min;
}

double getmax(double val[],int len)
{
	int maxnum,i;
	double max;
	max=val[0];
	maxnum=0;
	for(i=1;i<len;i++)
	{
		if (val[i]>max)
          {
          	max=val[i];
          	maxnum=i;
		  }
	}
	return max;
}

vector findgeomax(GEO geo,int NUM)
{
	int n[3]={0},i,j;
	vector temp;temp=vemery;
	for (i=0;i<3;i++){temp[i]=geo[0][i];}
	for(j=0;j<3;j++){
		for(i=1;i<NUM;i++)
		{	
		  if (geo[i][j]>temp[j])
	      {
	      	temp[j]=geo[i][j];
	      	n[j]=i;
		  }
		}	
	}
	return temp;	
}

vector findgeomin(GEO geo,int NUM)
{
	int n[3]={0},i,j;
	vector temp;temp=vemery;
	for (i=0;i<3;i++){temp[i]=geo[0][i];}
   
	for(i=1;i<NUM;i++)
	{
		for(j=0;j<3;j++){
		 if (geo[i][j]<temp[j])
          {
          	temp[j]=geo[i][j];
          	n[j]=i;
		  }
		}	
	}
	return temp;	
}
double sum(vector v,int len)
{
	double out=0;
	int i;
	
	for(i=0;i<len;i++)
	{
		out+=v[i];
	}
	return out;
}

void vprint(vector v)
{
	int i;
	for(i=0;i<SIZE;i++)
	{
		printf("%lf ",v[i]);
	}
	printf("\n");
}

double Abs(double x)
{
	if(x<0)
	  x=-x;
	  
	return x;
}
double gaussian(double u)     //用Box_Muller算法产生高斯分布的随机数
{
   double r,t,z,x;
   double s1,s2;
   s1=(1.0+rand())/(RAND_MAX+1.0);
   s2=(1.0+rand())/(RAND_MAX+1.0);
   r=sqrt(-2*log(s2)/log(2.71828));
   t=2*3.14*s1;
   z=r*cos(t);
   x=u+z*0.2;
   if(x<0) x=0.1;
   if(x>1) x=1;
   return x;
}
