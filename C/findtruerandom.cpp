#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "lib\oplocal.h"
void srand ( unsigned int seed );    
int main(){
 	 
 	int i,j;
    int *xu;
 	for (i=0;i<6;i++)
 	{
     xu=randperm(10);
	 for (j=0;j<10;j++)
 	 {  
 		printf("%d ",xu[j]);
	 }
	 printf("\n");
     
	}
}
