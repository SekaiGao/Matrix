#include<iostream>
#include"matrix.h"
using namespace std;
int main()
{	/*
	基本运算
	*/
	cout<<fixed;
	//初始化
	mat<ld>A(eye(5,5)),B(rand(5,5)),C{1,2,3,4},D{{1},{}},E(4,4),F({1,2},1);
	//改变大小
	A.Reshape(6,4);
	B.re(4,4);
	//赋值
	F=hilb(5);
	D={{0.1,0.2,0.3},{0.6,0.5,0.4},{1,0.5,0.9}};
	//销毁
	A.Destroy();
	//索引
	ld m_33=B[3][3];
	ld m_13=B(1,3);
	mat A_part=A({3,5},{2,4});
	//获取矩阵行列大小及分配内存大小
	int r=A.r(),c=A.c();
	int mr=A.max_row(),mc=A.max_col();
	auto[r0,c0]=C.size();
	//矩阵加减法
	B+=E;
	B=B+E;
	B-=eye(4,4);
	B=B-hilb(4);
	//转置
	mat<ld>Ct(C.t());
	//矩阵乘法
	auto C1=dot(C,{{0,r0},{0,c0}},Ct,{{0,c0},{0,r0}});
	auto C2=dot(C({0,r0},{0,c0}),Ct({0,c0},{0,r0}));
	Ct=C*Ct;
	C*=B;
	//数乘
	C*=Ct[0][0];
	Ct=1/Ct[0][0]*C;
	//幂
	mat<ld>Dp=D^100,Hs;
	Hs=hilb(4)^0.5;
	//行列式
	cout<<"|D|\n"<<D.Det()<<endl;
	//获取对角元
	mat<ld>K=Hs.diag();
	//行向量,列向量
	auto vec1(Dp.vrow(2)),vec2(Dp.vcol(1));
	//矩阵清0
	Dp.set(0);
	//迹
	D.tr();
	//平均值&求和
	D.mean();
	D.sum(1,2);
	//子矩阵
	A=Hs.cuts(1,2,1,2);
	//矩阵拼接
	mat<ld>A0=hilb(3).cat(eye(3,3));
	//排序
	mat<ld>B1=B.sort();
	//输出矩阵
	cout<<"A0\n"<<A0<<"\nB\n"<<B<<"\nC\n"
	<<C<<"\nHs\n"<<Hs<<"\nHs^2\n"<<Hs*Hs
	<<"\nvec2\n"<<vec2<<endl;

	
	/*
	矩阵分解
	*/
	A=hilb(3);
	cout<<"A\n"<<A<<endl;
	mat<ld>L,U,U1,Q,R,S,V,L1;
	A.LU(L,U);
	cout<<"L\n"<<L<<"U\n"<<U<<endl;
	A.QR(Q,R);
	cout<<"Q\n"<<Q<<"R\n"<<R<<endl;
	A.SVD(U1,S,V);
	cout<<"U1\n"<<U1<<"S\n"<<S<<"V\n"<<V<<endl;
	A.usf_qr(Q,R);
	cout<<"Q\n"<<Q<<"R\n"<<R<<endl;
	A.Cholesky(L1);
	cout<<"L1\n"<<L1<<endl;

	/*
	逆
	*/
	mat<ld>Ainv(A.Inv()),Apinv(A.Pinv());
	cout<<Ainv<<endl;//逆
	cout<<Apinv<<endl;//伪逆
	A0=A0.rref();//初等行变换
	cout<<A0.cuts(0,2,3,5)<<endl;

	/*
	解方程Ax=b
	*/	
	mat<ld>b{{1,1,1}};
	cout<<A.gaussj(b).t()<<endl;
	cout<<b.lusolve(L,U).t()<<endl;
	cout<<Ainv*b.t()<<endl;
	cout<<Apinv*b.t()<<endl;
	
	/*
	奇异值,特征值
	*/
	cout.precision(16);
	cout<<A.singular().t()<<endl;
	cout<<A.eig().t()<<endl;

	/*
	其他方法参考 matrix.h
	*/

	return 0;
}