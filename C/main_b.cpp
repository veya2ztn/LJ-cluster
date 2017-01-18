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

void mainnn(int NUM,double sup,double inf,double opt) 
{
	
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
    memcpy(S_print,S_old,gem);
    E_sup=0;E_inf=va;E_delta=(E_sup-E_inf)/divide;
    double E_new=va,E_old=va,E_GLM=va;
    double E_print=va; 
    
    char str[4],str2[4];
    
    int count3;
	count3=-opt/E_delta; 
	count3=count3/200;
	if(count3<10) count3=10;
	
	int timerecard=0;
	double timemax=0,timetemp=0;
	
	int timerecard2=0;
	double timeloop=0;
	
	clock_t start, finish,loopstart,loopend;
	loopstart=clock(); 
	
    while (mark<10)
    {
                
			    count2++;
				
					        
	            memcpy(S_new,gtemp=WL_randmove1_b(S_old,1,(sup-inf)/3,1,NUM,sup,inf),gem);free(gtemp);  
	            E_new=value(S_new,NUM);
	            gProject(S_new,NUM,sup,inf);
                E_new=vvalue(S_new,NUM);

		        if(E_new>=E_sup){}
		        else if(E_new<=E_inf)
			       {			  
					  memcpy(S_GLM,S_new,gem);
					  if(timerecard==100)
					  {
					  	 printf( "one BFGS_b algorithm need at most %0.2f seconds\n", timemax);
						 timerecard++;   
					  }
					  if(timerecard<100)
					  {
					  	timerecard++;
					  	start = clock(); 
					  	E_GLM=op_local_b_f(S_GLM,NUM,sup,inf);
					  	finish = clock(); 
					  	timetemp=(double)(finish - start) / CLOCKS_PER_SEC;
					  	if(timetemp>timemax)
					  	  timemax=timetemp;
					  }
					  else{E_GLM=op_local_b(S_GLM,NUM,sup,inf);}
					  
					  if (pow((E_GLM-E_print),2)<0.001) mark++;
					  if (E_GLM<E_print) {E_print=E_GLM;memcpy(S_print,S_GLM,gem);}
					   
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
	            
	           if(count2>10000)
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
				if(E_inf<opt/2)
				   break;

	}
//    printf("\n");
    printf("the opt energy is %lf\n",opt);
    printf("the min energy is %lf\n",E_print);
    printf("the min-energy gem is\n");
    gprint(S_print,NUM);
    sprintf(str, "%d", NUM);
	
    fileprint(S_print,NUM,str);
}
    
int main()
{
    
	int num=18;
    double sup=1;
    double inf=-1;
    double E_opt; 
    double opt[50]={
	0,-1,-3,-6,-9,-12,-16.505384,-19.821489,-24.113360,-28.422532,
	-32.765970,-37.967600,-44.326801,-47.845157,-52.322627,
	-56.815742,-61.317995,-66.530949,-72.659782,-77.177043,
	-81.684571,-86.809782,-92.844472,-97.348815,-102.372663,
	-108.315616,-112.873584,-117.822402,-123.587371,-128.286571,
	-133.586422,-139.635524,-144.842719,-150.044528,-155.756643,
	-161.825363,-167.033672,-173.928427,-180.033185,-185.249839,
	-190.536277,-196.277534,-202.364664,-207.688728,-213.784862,
	-220.680330,-226.012256,-232.199529,-239.091864,-244.549926};
	
	while (num<=20)
	{   
	    printf("the num at now is %d\n",num);
		E_opt=opt[num-1];    
		mainnn(num,sup,inf,E_opt);
	    num++;       
	}
}
