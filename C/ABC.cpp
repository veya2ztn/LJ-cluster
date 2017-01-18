#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>

#include "lib\oplocal.h"

#define gemery (double(*)[3])malloc(sizeof(double)*3*NUM)

#define SN 20
#define NUM 20
#define glimit 3
#define gmax 100

#define gem sizeof(double)*3*NUM
#define mem sizeof(double)*3*3
#define vem sizeof(double)*3
typedef double (*matrix)[3];  
typedef double  *vector;  
typedef double (*GEO)[3];  
main()
{
	int i,gnum=0,m,n;
	double sup=2,inf=-2,E_old,E_new;
	double G1store[SN][NUM][3];double Val1[SN];
	double G2store[SN][NUM][3];double Val2[SN];
	double Vdelta[SN];int limitturn[SN]={0};
	double E_GLM=0;
	int *xu;xu=(int * )malloc(sizeof(int)*SN);
	int *xtemp;
	GEO gtemp,gtemp1,gtemp2,gtemp3,gtemp4;
	
	double Vtemp[NUM][3];
	double tp1,tp2,tp3,p1,p2,p3,tpsum;
	
	double F,goodmintemp,bestmintemp;
	int goodminnum,bestminnum; 
//=====================initialization_begin====================
	for(i=0;i<SN;i++){memcpy(G1store[i],gtemp=initialization_b(NUM,sup,inf),gem);free(gtemp);}
	for(i=0;i<SN;i++){Val1[i]=op_local_b_veryfast(G1store[i],NUM,sup,inf);}
//======================initialization_end======================   
//==========================loop_start==========================
while(gnum<gmax)
{
//      ==================EMsearch_start=====================
    for(i=0;i<SN;i++){
    	E_old=Val1[i];
    	memcpy(xu,xtemp=randperm(SN),sizeof(int)*SN);free(xtemp);
    	gtemp1=G1store[xu[1]];gtemp2=G1store[xu[2]];gtemp3=G1store[xu[3]];
    	tp1=Val1[xu[1]];tp2=Val1[xu[2]];tp3=Val1[xu[3]];tpsum=tp1+tp2+tp3;
    	tp1=tp1/tpsum;tp2=tp2/tpsum;tp3=tp3/tpsum;
    	p1=1/3+tp2+tp3-2*tp1;p2=1/3+tp1+tp3-2*tp2;p3=1/3+tp1+tp2-2*tp3;
    	for(m=0;m<NUM;m++){
    		for(n=0;n<3;n++){
    	    	Vtemp[m][n]=p1*gtemp1[m][n]+p2*gtemp2[m][n]+p3*gtemp3[m][n];
			}
		}
		E_new=op_local_b_veryfast(Vtemp,NUM,sup,inf);
		if (E_new<E_old)
		   {
		   	memcpy(G2store[i],Vtemp,gem);
		    Val2[i]=E_new;
		   }
		else
		   {
		   	memcpy(G2store[i],G1store[i],gem);
		    Val2[i]=E_old;
		   } 
	}
	gnum++;
//      ==================EMsearch_end======================= 
//      ==================OLsearch_start=====================
    memcpy(xu,xtemp=randperm(SN),sizeof(int)*SN);free(xtemp);
//                  ====find the good one===
    goodmintemp=Val2[0];goodminnum=0;
    for(i=1;i<5;i++){if(Val2[xu[i]]<goodmintemp) {goodmintemp=Val2[xu[i]];goodminnum=xu[i];}}
//                  ====find the best one===
    bestmintemp=Val2[0];bestminnum=0;
    for(i=1;i<SN;i++){if(Val2[i]<bestmintemp) {bestmintemp=Val2[i];bestminnum=i;}}
    if (E_GLM>bestmintemp) E_GLM=bestmintemp;
//                  ====start to onlook ====
    memcpy(xu,xtemp=randperm(SN),sizeof(int)*SN);free(xtemp);
    for(i=1;i<5;i++){if(xu[i]==goodminnum) xu[i]=xu[6];} 
    gtemp1=G2store[xu[1]];gtemp2=G2store[xu[2]];gtemp3=G2store[xu[3]];gtemp4=G2store[xu[4]];
	F=random(1);
    if(random(1)<0.5)
    {
    	for(m=0;m<NUM;m++){
    		for(n=0;n<3;n++){
    	    	Vtemp[m][n]=G2store[goodminnum][m][n]+F*(gtemp1[m][n]+gtemp2[m][n]-gtemp3[m][n]-gtemp4[m][n]);
			}
		}
	}
	else
	{
		for(m=0;m<NUM;m++){
    		for(n=0;n<3;n++){
    	    	Vtemp[m][n]=G2store[bestminnum][m][n]+F*(gtemp1[m][n]+gtemp2[m][n]-gtemp3[m][n]-gtemp4[m][n]);
			}
		}
	}
	
    E_new=op_local_b_veryfast(Vtemp,NUM,sup,inf);
	
	if (E_new<goodmintemp)
	   {
	   	memcpy(G2store[goodminnum],Vtemp,gem);
	    Val2[goodminnum]=E_new;
	   }
	   
	if (E_GLM>E_new) E_GLM=E_new;
//      ==================OLsearch_end======================= 
//      ==================SCsearch_start=====================
    for(i=0;i<SN;i++)
	{
		
		if(Abs(Val2[i]-Val1[i])<0.1){
		   limitturn[i]++;
		   if(limitturn[i]>=glimit){
		   	   memcpy(G2store[i],gtemp=initialization_b(NUM,sup,inf),gem);free(gtemp);
		   	   Val2[i]=op_local_b_veryfast(G2store[i],NUM,sup,inf);
		   }
		}
		else{limitturn[i]=0;}
	}
        
    memcpy(G1store,G2store,sizeof(double)*3*NUM*SN); 
    memcpy(Val1,Val2,sizeof(double)*SN); 
//    vprint(Val2,SN);  
	printf("the GLM is at %d turn is %lf\n",gnum,E_GLM);
	
//      ==================SCsearch_end======================= 

}   
//==========================loop_end=============================
}
