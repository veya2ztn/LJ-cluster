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

int main{
	NUM=31;
	
    GEO S_old=gemery,S_GLM=gemery,S_new=gemery,S_print=gemery;
    
    double geo[][3]={}

}
