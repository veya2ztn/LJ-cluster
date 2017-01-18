#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>

#include "lib\oplocal.h"

#define divide 20
#define DIV 21
#define gemery (double(*)[3])malloc(sizeof(double)*3*NUM)

#define gem sizeof(double)*3*NUM
#define mem sizeof(double)*3*3
#define vem sizeof(double)*3
typedef double (*matrix)[3];  
typedef double  *vector;  
typedef double (*GEO)[3];  

int main() 
{
    int NUM=18;
    double sup=0.4;
    double inf=-0.4;
    double E_opt=-77.177043;
	
	GEO S_old=gemery,S_GLM=gemery,S_new=gemery,S_print=gemery;
    
    double E_delta;
    double E_inf;
	double E_sup;
    
    
    int old_inf=divide/2;
    int old_sup=old_inf+1;
    int new_inf,new_sup;
    
    double lnOmiga_list[DIV]={0},H[DIV]={0},E_mem[1000]={0};
    double lnOmiga_old,lnOmiga_new;
    double lnf=1;

    int ii=0,i;
    double unit;
    double scale;
    double r,sca_temp;
    int count=0,count2=0;
    
    
    int mark=0;
    vector lll=vemery;
    vector temp1;
    vector temp2;
    vector vtemp;
    GEO gtemp;
    
    memcpy(S_old,gtemp=initialization_b(NUM,sup,inf),gem);free(gtemp);
    
    double va;
    va=value(S_old,NUM);  
    gProject(S_old,NUM,sup,inf);
    va=vvalue(S_old,NUM);
    if (va>0.02*E_opt) va=0.02*E_opt;
    memcpy(S_print,S_old,gem);
    E_sup=0;E_inf=va;E_delta=(E_sup-E_inf)/divide;
    double E_new=va,E_old=va,E_GLM=va;
    double E_print=va; 
    
    int count3;
	count3=-E_opt/E_delta; 
	count3=count3/200;
	if(count3<10) count3=10;
	
	int timerecard=0;
	double timemax=0,timetemp=0;
	
	int timerecard2=0;
	double timeloop=0;
	
	clock_t start, finish,loopstart,loopend;
	loopstart=clock(); 
    
   int countlimit=1000;
    
    while (mark<10)
    {
			    count2++;
			    r=3*gaussian(1.0*count2/countlimit);  
	            memcpy(S_new,gtemp=WL_randmove1_b(S_old,r,(sup-inf)/3,r,NUM,sup,inf),gem);free(gtemp);  
	            E_new=value(S_new,NUM);
	            gProject(S_new,NUM,sup,inf);
                E_new=vvalue(S_new,NUM);          
                
		        if(E_new>=E_sup){}
		        else if(E_new<=E_inf)
			       {			  
					  memcpy(S_GLM,S_new,gem);
					  if(timerecard==100)
					  {
					  	 printf( "\n one BFGS_b algorithm need at most %0.2f seconds\n", timemax); 
					  	 timerecard++;
					  }
					  if(timerecard<100)
					  {
					  	timerecard++;
					  	start = clock(); 
					  	E_GLM = op_local_b_f(S_GLM,NUM,sup,inf);
					  	finish= clock(); 
					  	timetemp=(double)(finish - start) / CLOCKS_PER_SEC;
					  	if(timetemp>timemax)
					  	  timemax=timetemp;
					  }
					   else{E_GLM=op_local_b_f(S_GLM,NUM,sup,inf);}
					  if (pow((E_GLM-E_print),2)<0.001) mark++;
					  printf("->");
					  if (E_GLM<E_print)
					   {  
					      printf("\n");
						  E_print=E_GLM;
					      memcpy(S_print,S_GLM,gem);
					      printf("oh~oh~oh the min is %lf\n",E_print);
					      gprint(S_print,NUM);
					      printf("========================================\n");
					   }
					   
//					   if (mark>1000)
//					      goto here;
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
							   	memcpy(S_old,S_new,sizeof(double)*3*NUM);
							   	E_old=E_new;
							   	old_sup=new_sup;
							   	old_inf=new_inf;
						   }
				   }
	            scale=(int)(scale+0.5)>(int)scale?(int)scale+1:(int)scale;
	            lnOmiga_list[old_inf]+=scale*lnf; 
	            lnOmiga_list[old_sup]+=(1-scale)*lnf;
	            H[old_sup]+=1-scale;
	            H[old_inf]+=scale;
	            
	            if(count2>countlimit)
	               { 
                     loopend=clock();
	               	 count2=0;
	               	 timeloop=(double)(loopend - loopstart) / CLOCKS_PER_SEC;
	               	 printf("\n");
	               	 printf("one loop need around %0.2f seconds, at most 20 loop\n", timeloop); 
	               	 printf("\n");
	               	 if(timetemp<timeloop)
	               	 {
						 printf("the min energy now is %lf\n",E_print);
		               	 printf("the min-energy gem now is\n");
	                     gprint(S_print,NUM);
                     }
	               	 for(i=0;i<DIV;i++){lnOmiga_list[i]=H[i]=0;};lnf=1;
	               	 E_sup-=count3*E_delta;
	               	 E_inf-=count3*E_delta;
	               	 mark=0;
	               	 loopstart=clock();
				   }
				if(E_inf<E_opt/2)
				   break;
	}

    printf("the min energy is %lf\n",E_print);
    printf("the min-energy gem is\n");
    gprint(S_print,NUM);
    fileprint(S_print,NUM,"18-0.4");
}


