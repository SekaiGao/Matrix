#include <iostream>
#include "matrix.h"
using namespace std;
int main()
{ /*
	基本运算
	*/
	cout.precision(3);
	cout << fixed;
	//初始化
	mat<ld> A(hilb(5)),U,S,V;
	A.SVD(U,S,V);
	cout << "U\n" << U << "\nS\n" << S << "\nV\n" << V << endl;
	return 0;
}