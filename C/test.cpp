#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>

#include "lib\oplocal.h"

#define NUM 15

#define divide 20
#define DIV 21
#define gemery (double(*)[3])malloc(sizeof(double)*3*NUM)

#define gem sizeof(double)*3*NUM
#define mem sizeof(double)*3*3
#define vem sizeof(double)*3
typedef double (*matrix)[3];  
typedef double  *vector;  
typedef double (*GEO)[3];  

int main(){
 	srand((unsigned)time(NULL));  
    
//	printf("GEO is list \n");gprint(geo,NUM);
    
    
	GEO S_old=gemery,S_GLM=gemery,S_new=gemery;
	
    double va;
    va=value(S_old,NUM);
    
    double E_new=va,E_old=va,E_GLM=0;
    double E_inf,E_sup;
    double E_delta;
    
    E_inf=-10;E_sup=0;E_delta=(E_sup-E_inf)/divide;
    
    int old_inf=divide/2;
    int old_sup=old_inf+1;
    int new_inf,new_sup;
    
    double lnOmiga_list[DIV]={0},H[DIV]={0},E_mem[1000]={0};
    double lnOmiga_old,lnOmiga_new;
    double lnf=1;

    int mark=1,ii=0,i;
    
    double unit;
    double scale;
    double r,sca_temp;
    int count=0;
    double E_opt=-52.2;
    
    vector lll=vemery;
    vector temp1;
    vector temp2;
    vector vtemp;
    GEO gtemp;
    
    memcpy(S_old,gtemp=initialization(NUM,1),gem);free(gtemp);
    while (E_GLM>E_opt)
    {

		   memcpy(lll,vtemp=v_v(temp1=findgeomax(S_old,NUM),temp2=findgeomin(S_old,NUM)),vem);free(vtemp);free(temp1);free(temp2);
	       unit=pow(lll[0]*lll[1]*lll[2]/NUM,1.0/3);     
	       sca_temp=exp(-0.001*count);count++;
		   r=unit*sca_temp;
		   if (sca_temp<0.1)
		      count=0;
		   memcpy(S_new,gtemp=WL_randmove1(S_old,r,NUM),gem);free(gtemp);
	       E_new=value(S_new,NUM);
	       if(E_new>=E_sup){}
	       else if(E_new<=E_inf)
	       {
	       	   
			  memcpy(S_GLM,S_new,gem);
			  op_local(S_GLM,NUM);
			  E_GLM=value(S_GLM,NUM);
			  if (E_GLM<-50){
			  
				 op_local(S_GLM,NUM);
				 E_GLM=value(S_GLM,NUM);
				 printf("scale=%lf global min=%lf\n",sca_temp,E_GLM);  }  	
		      if (E_GLM<E_opt)
			    goto here;	
	       }
	       else
	       {
		       	new_inf=floor((E_new-E_inf)/E_delta);
		       	new_sup=new_inf;
		       	scale=(E_new-E_inf)/E_delta-new_inf;
		       	
		       	lnOmiga_old = (1-scale)*lnOmiga_list[old_inf]+scale*lnOmiga_list[old_sup];
		        lnOmiga_new = (1-scale)*lnOmiga_list[new_inf]+scale*lnOmiga_list[new_sup];
			    
			   if (log(random(1))<lnOmiga_old-lnOmiga_new)
			   {
			   	memcpy(S_old,S_new,gem);
			   	E_old=E_new;
			   	old_sup=new_sup;
			   	old_inf=new_inf;
			   }
		   }
	       scale=(int)(scale+0.5)>(int)scale?(int)scale+1:(int)scale;
	       H[old_sup]+=1-scale;
	       H[old_inf]+=scale;
	       lnOmiga_list[old_inf]+=scale*lnf; 
	       lnOmiga_list[old_sup]+=(1-scale)*lnf;
	 }
	here:;
	printf("the global min is %lf\n",E_GLM);
	printf("the global min GEO is list\n");gprint(S_GLM,NUM);

}

