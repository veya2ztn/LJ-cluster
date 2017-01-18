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
	int NUM=20;double sup=1;double inf=-1;double opt=-61.317995;
    GEO S_old=gemery,S_GLM=gemery,S_new=gemery,S_print=gemery;
    vector lll=vemery;
    vector temp1;
    vector temp2;
    vector vtemp;
    GEO gtemp;
    int ii=0,i;
    double E_opt=opt;

    
    double E_delta;
    double E_inf;
	double E_sup;
    
    
    int old_inf=divide/2;
    int old_sup=old_inf+1;
    int new_inf,new_sup;
    
    double lnOmiga_list[DIV]={0},H[DIV]={0},E_mem[1000]={0};
    double lnOmiga_old,lnOmiga_new;
    double lnf=1;

    
    double unit;
    double scale;
    double r,sca_temp;
    int count=0,count2=0;
    
    
    int mark=0;
    
    char str[4],str2[4];
    
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
	count3=count3/40;
	if(count3<10) count3=10;

	int timerecard=0;
	double timemax=0,timetemp=0;
	
	int timerecard2=0;
	double timeloop=0;
	
	clock_t start, finish,loopstart,loopend;
    loopstart=clock();
    
    int countlimit=10000;
    int loopcounts=0;
    
    int threebreak=1;
    double E_save=0;
    
    
    
    
    
    
    while (mark<10)
    {
			    
				count2++;
			    r=3*gaussian(1.0*count2/countlimit);  
	            memcpy(S_new,gtemp=WL_randmove1_b(S_old,r,(sup-inf)/3,r,NUM,sup,inf),gem);free(gtemp);  
	            E_new=value(S_new,NUM);
	            gProject(S_new,NUM,sup,inf);
                E_new=vvalue(S_new,NUM); 
				if (count2>10000)
				{
					loopend=clock();
					timeloop=(double)(loopend - loopstart) / CLOCKS_PER_SEC;
					printf("one loop need around %0.2f seconds\n", timeloop); 
					loopstart=clock();
					count2=0;
				}
			
}
            
}


