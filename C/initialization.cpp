#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include "lib\universe_math.h"

#define NUM 40
#define ROW 15
#define STR 2

int main(){
    double geo[NUM][3]={0},temp;
    int i,j,nn=0,k;
    
    int basis_L,last_N,outside_N;
    double center_L;
    

    initialization(geo,NUM,STR);
//	if(STR==1)
//    {
//    	for(i=0;i<NUM;i++)
//    	{
//    		shift(geo[i],2);
//		}
//	}
//    else if (STR==2)
//	{
//		if(NUM<27)
//		   printf("Num is too small, suggest to use sphere\n");
//		else
//		{
//			basis_L= floor(pow(NUM,1.0/3));
//			outside_N=floor((NUM-pow(basis_L,3))/6);
//			last_N=NUM-pow(basis_L,3)-outside_N*6; 
//			for(i=0;i<basis_L;i++){
//			   	for(j=0;j<basis_L;j++){
//			   		for(k=0;k<basis_L;k++){
//			   			geo[nn][0]=i;
//			   			geo[nn][1]=j;
//			   			geo[nn][2]=k;
//			   			nn++;
//					   }
//				   }
//			   }
//		    
//			center_L=(basis_L+1)/2; 
//			if(outside_N!=0)
//			{
//				for(i=0;i<outside_N;i++)
//				{
//					geo[nn][0]=2*center_L;                      geo[nn][1]=center_L+basis_L*(2*random(1)-1);geo[nn][2]=center_L+basis_L*(2*random(1)-1);nn++;
//					geo[nn][0]=-1;                              geo[nn][1]=center_L+basis_L*(2*random(1)-1);geo[nn][2]=center_L+basis_L*(2*random(1)-1);nn++;
//					geo[nn][0]=center_L+basis_L*(2*random(1)-1);geo[nn][1]=2*center_L;                      geo[nn][2]=center_L+basis_L*(2*random(1)-1);nn++;
//					geo[nn][0]=center_L+basis_L*(2*random(1)-1);geo[nn][1]=-1;                              geo[nn][2]=center_L+basis_L*(2*random(1)-1);nn++;
//					geo[nn][0]=center_L+basis_L*(2*random(1)-1);geo[nn][1]=center_L+basis_L*(2*random(1)-1);geo[nn][2]=2*center_L;                      nn++;
//					geo[nn][0]=center_L+basis_L*(2*random(1)-1);geo[nn][1]=center_L+basis_L*(2*random(1)-1);geo[nn][2]=-1;                              nn++;
//				}
//		    }      
//		    
//		    if(last_N!=0)
//		    {
//		    	for(i=0;i<last_N;i++)
//		    	{
//		    		geo[nn][0]=center_L;geo[nn][1]=center_L;geo[nn][2]=center_L;
//					shift(geo[nn],center_L);
//					nn++;
//				}
//			}
//			
//		}
//	} 
//	else
//	  printf("please choose the right para\n");
	  	for(i=0;i<NUM;i++)
	  { 
	  	for(j=0;j<3;j++)
	  	  {
	  	  	printf("%lf ",*(*(geo+i)+j));
		  }
		printf("\n");  
	  }
	  printf("\n");  


}

