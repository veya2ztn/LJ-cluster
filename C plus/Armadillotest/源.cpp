#include <iostream>

#include "test.h"
#include<time.h>
using namespace std;
using namespace arma;

 
int main() {
	
	clock_t  start, finish;
	
	// for matrices with real elements

	mat A = randu<mat>(50, 50);
	mat B = A.t()*A;  // generate a symmetric matrix

	vec eigval;
	mat eigvec;

	

	start = clock();
	eig_sym(eigval, eigvec, B);
	finish = clock();

	cout << (double)(finish - start) / CLOCKS_PER_SEC << endl;

	// for matrices with complex elements

	cx_mat C = randu<cx_mat>(50, 50);
	cx_mat D = C.t()*C;

	vec eigval2;
	cx_mat eigvec2;
	start = clock();
	eig_sym(eigval2, eigvec2, D);
	finish = clock();
	cout << (double)(finish - start) / CLOCKS_PER_SEC << endl;




	cin.get();
	return 0;
}