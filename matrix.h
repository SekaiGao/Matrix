/*
* This file is just a matrix(linear algebra) library written for fun.
* It may contain many bugs.
* If you want to compute faster and more accurately,
* you should choose: Eigen/Lapack(C/C++), numpy/scipy/torch(python), matlab
* Copyright (C) 2022 Sekai Gao
* All rights reserved.
*/
#ifndef MATRIX_H_
#define MATRIX_H_
#include<cmath>
#include<time.h>
#include<complex>
#include<initializer_list>
typedef long double ld;
typedef std::complex<long double> cld;
const long double pi=3.14159265358979323846;

template<class T=ld>
class mat
{
	//当前行列
	int row=0;
	int col=0;
	//已分配内存所占行列
	int max_r=0;
	int max_c=0;
	//索引方式:arr[i][j]
	T**arr=nullptr;
public:

	mat(){}

	mat(int row,int col):row(row),col(col),max_r(row),max_c(col)//初始化行列
	{//没有实现稀疏矩的三元组存储即与稠密矩阵的转换
	//向量可以表示为mat(1,N)或者mat(N,1),考虑局部性,行向量更快
	//矩阵初始化会自动清0
		if(row<0||col<0)
		{std::cerr<<"\nInit failure.\n";
			exit(EXIT_FAILURE);
		}
		if(!(arr=(T**)malloc(row*sizeof(T*))))
		{std::cerr<<"\nMalloc failure.\n";
			exit(EXIT_FAILURE);
		}
		for(int i=0;i<row;++i)
		if(!(arr[i]=(T*)calloc(col,sizeof(T))))//calloc会自己清0
		{std::cerr<<"\nCalloc failure.\n";
			exit(EXIT_FAILURE);
		}
	}

	mat(const std::initializer_list<int>&li,bool axis=0)//初始化列表形式
	{
	/*
	初始化方式: mat A{...} A({...},1)
	赋值方式: A={...}
	*/
	if(axis)
	{
		row=li.size();
		col=1;
	}
	else
	{
		row=1;
		col=li.size();
	}
	if(!(arr=(T**)malloc(row*sizeof(T*))))
	{std::cerr<<"\nMalloc failure.\n";
		exit(EXIT_FAILURE);
	}
	max_r=row;
	max_c=col;
	int i=0;
	if(axis)
	for(auto&x:li)
	{
		if(!(arr[i]=(T*)calloc(col,sizeof(T))))//calloc会自己清0
		{std::cerr<<"\nCalloc failure.\n";
			exit(EXIT_FAILURE);
		}
		arr[i++][0]=x;
	}
	else 
	{
		if(!(arr[0]=(T*)calloc(col,sizeof(T))))//calloc会自己清0
		{std::cerr<<"\nCalloc failure.\n";
			exit(EXIT_FAILURE);
		}
		i=0;
		for(auto&x:li)
		arr[0][i++]=x;
	}
	}

	mat(const std::initializer_list<std::initializer_list<T>>&li):row(li.size())//初始化列表形式
	{
	/*
	初始化方式: mat A{{},{},...} A({{},{},...})
	赋值方式: A={{},{},...}
	*/
	if(!(arr=(T**)malloc(row*sizeof(T*))))
	{std::cerr<<"\nMalloc failure.\n";
		exit(EXIT_FAILURE);
	}
	max_r=row;
	col=0;
	for(auto&x:li)//col为最大列数
    if(int(x.size())>col)
    col=x.size();
	max_c=col;
	int i=0;
	for(auto&x:li)
	{
		if(!(arr[i]=(T*)calloc(col,sizeof(T))))//calloc会自己清0
		{std::cerr<<"\nCalloc failure.\n";
			exit(EXIT_FAILURE);
		}
		int j=0;
		for(const T&e:x)
		arr[i][j++]=e;
		++i;
	}
	}

	mat(const mat&a)//()初始化
	{//A(B)
		re(a.row,a.col);
		for(int i=0;i<a.row;++i)
		for(int j=0;j<a.col;++j)
		arr[i][j]=a.arr[i][j];
		row=a.row,col=a.col;
	}

	bool re(int row,int col)//改变矩阵大小
	{//若新矩阵大小在已分配内存容量内,不重新分配内存
		if(row<=max_r&&col<=max_c)
		{
		this->row=row,this->col=col;
		return 1;
		}
		return Reshape(row,col);
	}

	bool Reshape(int row,int col)//重新分配内存
	{//重新分配内存会把原来的数据copy过去,但是多出来的列不是0(realloc不置0),行是0(calloc分配))
		if((row<=0||col<=0))//允许缩小内存
		return 1;
		if(this->row==0)
		{if(!(arr=(T**)malloc(row*sizeof(T*))))
		{std::cerr<<"\nMalloc failure.\n";
		return 0;
		}}
		else
		if(!(arr=(T**)realloc(arr,row*sizeof(T*))))//不可以realloc未分配内存的
		{std::cerr<<"\nRelloc failure.\n";
		return 0;
		}
		int min=this->row<row?this->row:row;
		for(int i=0;i<min;++i)
		if(!(arr[i]=(T*)realloc(arr[i],col*sizeof(T))))
		{std::cerr<<"\nRelloc failure.\n";
			return 0;
		}
		for(int i=min;i<row;++i)
		if(!(arr[i]=(T*)calloc(col,sizeof(T))))
		{std::cerr<<"\nCalloc failure.\n";
			return 0;
		}
		this->row=row,this->col=col;
		max_r=row,max_c=col;
		return 1;
	}

	void Destroy()//释放内存
	{
		if(arr)
		{//std::cout<<"~\n";
		for(int i=0;i<max_r;++i)
		free(arr[i]);
		free(arr);
		}
		row=0,col=0;
		max_r=0,max_c=0;
	}

	~mat()//析构
	{
		if(arr)
		{//std::cout<<'~'<<row<<' '<<col<<'\n';
		//#pragma omp parallel for
		for(int i=0;i<max_r;++i)
		free(arr[i]);
		free(arr);
		}
		row=0,col=0;
		max_r=0,max_c=0;
	}

	mat operator()(const std::initializer_list<int>&rSlice,const std::initializer_list<int>cSlice)//高级索引
	{//包头不包尾
		int r0,r1,c0,c1;
		if(rSlice.size()==0)
		r0=0,r1=row;
		else if(rSlice.size()==1)
		r0=*rSlice.begin(),r1=r0+1;
		else
		r0=*rSlice.begin(),r1=*(rSlice.begin()+1);
		if(cSlice.size()==0)
		c0=0,c1=col;
		else if(cSlice.size()==1)
		c0=*cSlice.begin(),c1=c0+1;
		else
		c0=*cSlice.begin(),c1=*(cSlice.begin()+1);
		if(r0<0||r0>row-1||r1<r0||r1>row)
		{
		std::cerr<<"\nIndex out of range.\n";
		exit(1);
		}
		if(c0<0||c0>col-1||c1<c0||c1>col)
		{
		std::cerr<<"\nIndex out of range.\n";
		exit(1);
		}
		int rt=r1-r0,ct=c1-c0;
		mat tem(rt,ct);
		for(int i=0;i<rt;++i)
		for(int j=0;j<ct;++j)
		tem[i][j]=arr[i+r0][j+c0];
		return tem;
	}

	T*&operator[](int i)//索引(返回引用,可更改返回的元素)
	{//或许可以设计python的负值索引
	return arr[i];
	}

	T&operator()(int i,int j)//检测边界
	{
	if(i<0||i>row-1||j<0||j>col-1)
	{
		std::cerr<<"\nIndex out of range.\n";
		exit(1);
	}
	return arr[i][j];
	}

	friend mat operator+(const mat&a,const mat&b)//重载 +
	{
	if(a.row!=b.row||a.col!=b.col)
	{
		std::cerr<<"\nAdd failure.\n";
		exit(1);
	}
	mat tem(a.row,a.col);//局部对象出函数就析构
	for(int i=0;i<a.row;++i)
	for(int j=0;j<a.col;++j)
	tem.arr[i][j]=a.arr[i][j]+b.arr[i][j];
	return tem;
	}

	void operator+=(const mat&b)
	{//由于没有局部对象的生成析构,+=事实上比+快
	if(row!=b.row||col!=b.col)
	{
		std::cerr<<"\nAdd failure.\n";
		exit(1);
	}
	for(int i=0;i<row;++i)
	for(int j=0;j<col;++j)
	arr[i][j]+=b.arr[i][j];
	}

	friend mat operator-(const mat&a,const mat&b)
	{
	if(a.row!=b.row||a.col!=b.col)
	{
		std::cerr<<"\nSub failure.\n";
		exit(1);
	}
	mat tem(a.row,a.col);
	for(int i=0;i<a.row;++i)
	for(int j=0;j<a.col;++j)
	tem.arr[i][j]=a.arr[i][j]-b.arr[i][j];
	return tem;
	}

	void operator-=(const mat&b)
	{
	if(row!=b.row||col!=b.col)
	{//比-快
		std::cerr<<"\nSub failure.\n";
		exit(1);
	}
	for(int i=0;i<row;++i)
	for(int j=0;j<col;++j)
	arr[i][j]-=b.arr[i][j];
	}

	friend mat operator*(const mat&a,const mat&b)//用AVX优化?
	{//矩阵乘法
	if(a.col!=b.row)
	{
		std::cerr<<"\nMulti failure.\n";
		exit(1);
	}
	mat tem(a.row,b.col);
	T s;
	//#pragma omp parallel for//50阶以上开多线程并行
	for(int i=0;i<a.row;++i)
	//#pragma omp parallel for
	for(int j=0;j<a.col;++j)
	{s=a.arr[i][j];
	for(int k=0;k<b.col;++k)
	tem.arr[i][k]+=s*b.arr[j][k];//利用空间局部性
	}
	return tem;
	}

	void operator*=(const mat&a)
	{//与*速度差不多
	if(col!=a.row)
	{
		std::cerr<<"\nMulti failure.\n";
		exit(1);
	}
	T tem[row][a.col],s;
	for(int i=0;i<row;++i)
	for(int j=0;j<a.col;++j)
	tem[i][j]=0;
	//#pragma omp parallel for//50阶以上开
	for(int i=0;i<row;++i)
	//#pragma omp parallel for
	for(int j=0;j<col;++j)
	{s=arr[i][j];
	for(int k=0;k<a.col;++k)
	tem[i][k]+=s*a.arr[j][k];
	}
	re(row,a.col);
	for(int i=0;i<row;++i)
	for(int j=0;j<col;++j)
	arr[i][j]=tem[i][j];
	}
	
	friend mat operator*(const T&a,const mat&b)
	{//数乘
	mat tem(b.row,b.col);
	for(int i=0;i<b.row;++i)
	for(int j=0;j<b.col;++j)
	tem.arr[i][j]=a*b.arr[i][j];
	return tem;
	}

	void operator*=(const T&a)
	{
	for(int i=0;i<row;++i)
	for(int j=0;j<col;++j)
	arr[i][j]*=a;
	}

	mat&operator=(const mat&a)//解决赋值时内存分配(深拷贝)
	{
		if(this==&a)
		return *this;
		re(a.row,a.col);//局部对象析构一次
		for(int i=0;i<a.row;++i)
		for(int j=0;j<a.col;++j)
		arr[i][j]=a.arr[i][j];
		return *this;
	}

	mat operator^(int x)//幂(快速幂))
	{
		if(col!=row)
		{
		std::cerr<<"\n mismatch.\n";
		exit(1);
		}
		if(x<1)
		{
		std::cerr<<"\n exp failure.\n";
		exit(1);
		}
		mat tem(row,row),p(*this);
		for(int i=0;i<row;++i)
		tem.arr[i][i]=1;
		while(x)
		{
			if(x&1)
			tem*=p;
			p*=p;
			x>>=1;
		}
		return tem;
	}

	mat operator^(double x)//幂(可分数幂)
	{//如果陷入死循环,大概率是有复特征值
		if(col!=row)
		{
		std::cerr<<"\n mismatch.\n";
		exit(1);
		}
		mat Q,E;
		usf_qr(Q,E);
		for(int i=0;i<col;++i)
		E[i][i]=pow(E[i][i],x);
		return Q*E*Q.t();//特征分解:A^n=Q*(E^n)*Q.T
	}

	friend mat dot(const mat&a,const mat&b)
	{
		return a*b;
	}

	friend mat dot(const mat&a,const std::initializer_list<std::initializer_list<int>>&rangeA,const mat&b,const std::initializer_list<std::initializer_list<int>>&rangeB)
	{//矩阵索引点乘(避免多次初始化与析构)
		auto*b1=rangeA.begin();
		int ra0=*(b1->begin()),ra1=*(b1->begin()+1);
		int ca0=*(b1+1)->begin(),ca1=*((b1+1)->begin()+1);
		b1=rangeB.begin();
		int rb0=*(b1->begin()),rb1=*(b1->begin()+1);
		int cb0=*(b1+1)->begin(),cb1=*((b1+1)->begin()+1);
		if(ra0<0||ra0>a.row-1||ra1<ra0||ra1>a.row)
		{
		std::cerr<<"\nIndex out of range.\n";
		exit(1);
		}
		if(ca0<0||ca0>a.col-1||ca1<ca0||ca1>a.col)
		{
		std::cerr<<"\nIndex out of range.\n";
		exit(1);
		}
		if(rb0<0||rb0>b.row-1||rb1<rb0||rb1>b.row)
		{
		std::cerr<<"\nIndex out of range.\n";
		exit(1);
		}
		if(cb0<0||cb0>b.col-1||cb1<cb0||cb1>b.col)
		{
		std::cerr<<"\nIndex out of range.\n";
		exit(1);
		}
		int r=ra1-ra0,ca=ca1-ca0,rb=rb1-rb0,c=cb1-cb0;
		if(ca!=rb)
		{
		std::cerr<<"\nMulti failure.\n";
		exit(1);
		}
		mat tem(r,c);
		T s;
		//#pragma omp parallel for//50阶以上开多线程并行
		for(int i=0;i<r;++i)
		//#pragma omp parallel for
		for(int j=0;j<ca;++j)
		{s=a.arr[i+ra0][j+ca0];
		for(int k=0;k<c;++k)
		tem.arr[i][k]+=s*b.arr[j+rb0][k+cb0];//利用空间局部性
		}
		return tem;
	}

	mat t()//转置
	{
		mat tem(col,row);
		for(int i=0;i<col;++i)
		for(int j=0;j<row;++j)
		tem.arr[i][j]=arr[j][i];
		return tem;
	}
	
	mat cut(int c)//代数余子式
	{//返回i行对应代数余子式
		mat tem(row-1,col-1);
		for(int i=0;i<tem.row;++i)
		for(int j=0;j<tem.col;++j)
		{
		tem.arr[i][j]=arr[i+1][j<c?j:j+1];
		}
		return tem;
	}
	
	mat cuts(int r1,int r2,int c1,int c2)//子矩阵
	{//返回row1——row2,col1——col2构成的子矩阵
		mat tem(r2-r1+1,c2-c1+1);
		for(int i=0;i<tem.row;++i)
		for(int j=0;j<tem.col;++j)
		{
		tem.arr[i][j]=arr[i+r1][j+c1];
		}
		return tem;
	}
	
	friend std::ostream&operator<<(std::ostream&os,const mat&s)//输出矩阵
	{
		if(s.row==0||s.col==0)
		os<<"null\n";
		for(int i=0;i<s.row;++i)
		{
		for(int j=0;j<s.col;++j)
		os<<s.arr[i][j]<<'\t';
		os<<std::endl;
		}
		return os;
	}

	mat diag()//对角元
	{//返会对角元组成向量
		int min=row<col?row:col;
		mat tem(1,min);
		for(int i=0;i<min;++i)
		tem[0][i]=arr[i][i];
		return tem;
	}

	long double det()//行列式 O(n!)
	{
		if(row!=col)
		{std::cerr<<"\ndet\nThis matrix has no det.\n";
		exit(1);
		}
		if(row==1)
		return arr[0][0];
		else if(row==2)
		return arr[0][0]*arr[1][1]-arr[0][1]*arr[1][0];
		else 
		{double sum=0;
		for(int i=0;i<col;++i)
		sum+=arr[0][i]*cut(i).det()*pow(-1,i);
		return sum;
		}
	}

	T tr()//迹
	{//等于所有特征值之和
		T sum=0;
		int m=row>col?col:row;
		for(int i=0;i<m;++i)
		sum+=arr[i][i];
		return sum;
	}
	
	void swap_row(int iloc,int jloc)//行交换
	{//交换矩阵A的iloc,jloc两行
	if(iloc==jloc)
	return;
	long double tem=0;
	for(int i=0;i<col;++i)
	{
		tem=arr[iloc][i];
		arr[iloc][i]=arr[jloc][i];
		arr[jloc][i]=tem;
	}
	}

	void set(T o=0)//矩阵置0
	{//堆上分配的内存似乎无法用memset清0
		for(int i=0;i<row;++i)
		for(int j=0;j<col;++j)
		arr[i][j]=0;
	}

	mat cat(const mat&B)//矩阵拼接
	{//把当前对象与B拼接起来
		if(row!=B.row)
		{
		std::cerr<<"\ncat\nunmatch.\n";
		exit(1);
		}
		mat tem(*this);
		tem.re(row,col+B.col);
		for(int i=0;i<row;++i)
		for(int j=0;j<B.col;++j)
		tem[i][j+col]=B.arr[i][j];
		return tem;
	}

	mat gaussj(mat&b)//gauss-jordan消元法
	{
	struct po{
	int i=0;
	int j=0;
	long double val=0;
	}pos;
	int ed[col+1]={0};
	mat A,x(1,col);
	A.re(row,col+1);
	int n=col;
	for(int i=0;i<n;++i)
	for(int j=0;j<n;++j)
	A[i][j]=arr[i][j];
	for(int i=0;i<n;++i)
	A[i][n]=b[0][i];//增广
	for(int k=0;k<n;++k)
	{pos.val=0;
	for(int i=0;i<n;++i)//选取最大元素所在行
	for(int j=0;j<n;++j)
	if(fabs(A[i][j])>fabs(pos.val)&&ed[i]!=1)
	pos.i=i,pos.j=j,pos.val=A[i][j];
	ed[pos.i]=1;//标记已化1的列
	
	for(int i=0;i<=n;++i)//将主元化为1
	A[pos.i][i]/=pos.val;
	for(int i=0;i<n;++i)//其他行减选取行
	{
		if(i==pos.i)
		continue;
		long double tem=A[i][pos.j];
		for(int j=0;j<=n;++j)
		A[i][j]-=tem*A[pos.i][j];
	}
	}
	for(int i=0;i<n;++i)
	{
	for(int j=0;j<n;++j)
	{if(A[i][j]==1)
	{
	x[0][j]=A[i][n];
	continue;
	}}
	}
	return x;
	/*
	当b不在Col A中时
	Ax=b最小二乘解:
	1. A'Ax=A'b
	2. A=QR, Rx=Q'b
	3. x=pinv(A)b
	*/
	}

	mat rref()//简化行阶梯形矩阵
	{//基于gaussj消元法
	struct po{
	int i=0;
	int j=0;
	long double val=0;
	}pos;
	int m=row>col?col:row;
	int ed[row]={0};
	mat A(*this),x(1,col);
	
	for(int k=0;k<m;++k)
	{pos.val=0;
	for(int i=0;i<row;++i)//选取最大元素所在行
	for(int j=0;j<m;++j)
	if(fabs(A[i][j])>fabs(pos.val)&&ed[i]!=1)
	pos.i=i,pos.j=j,pos.val=A[i][j];
	ed[pos.i]=1;//标记已化1的列
	
	for(int i=0;i<col;++i)//将主元化为1
	A[pos.i][i]/=pos.val;
	for(int i=0;i<row;++i)//其他行减选取行
	{
		if(i==pos.i)
		continue;
		long double tem=A[i][pos.j];
		for(int j=0;j<col;++j)
		A[i][j]-=tem*A[pos.i][j];
	}
	}
	for(int i=0;i<row;++i)
	for(int j=0;j<m;++j)
	if(A[i][j]==1)
	A.swap_row(i,j);
	return A;
	}

	int PLU(mat&L,mat&U,mat&P)//PLU分解
	{//PA=LU(列选主元LU)),数值稳定性比LU好
	if(row!=col)
	{std::cerr<<"\nPLU\nThis matrix is not a square.\n";
	exit(1);
	}
	struct pos//列选主元时存储行位置i与该处值val
	{
	int i=0;
	T val=0;
	};
	mat A;
	A=*this;
	L.re(row,row),U.re(row,row),P.re(row,row);
	int i,j,k,sw=1;
	for(i=0;i<row;++i)
	for(j=0;j<row;++j)
	P.arr[i][j]=(i==j)?1:0;
	for(k=0;k<row;++k)
	{pos po;
	po.i=k,po.val=A[k][k];
	for(i=k;i<row;++i)
	{
		if(fabs(po.val)<fabs(A[i][k]))//找最大列元素
		po.i=i,po.val=A[i][k];//存储列主元
	}
	A.swap_row(k,po.i);//把列主元行和k行交换
	P.swap_row(k,po.i);
	if(k!=po.i)
	sw*=-1;
	for(i=k+1;i<row;++i)//开始LU分解
	{
		A[i][k]/=po.val;
		for(j=k+1;j<row;++j)
		A[i][j]-=A[k][j]*A[i][k];
	}
	}
	for(i=0;i<row;++i)
	{
	for(j=0;j<row;++j)//将LU分解结果存入L,U中
	{
		if(i>j)
		L.arr[i][j]=A[i][j],U.arr[i][j]=0;
		else
		U.arr[i][j]=A[i][j],L.arr[i][j]=0;
	}
	L.arr[i][i]=1;//L的对称轴元素置1

	}
	return sw;
	}

	void LU(mat&L,mat&U)//LU分解
	{
	if(row!=col)
	{std::cerr<<"\nLU\nThis matrix is not a square.\n";
	exit(1);
	}
	mat A;
	A=*this;
	L.re(row,row),U.re(row,row);
	int i,j,k;
	for(k=0;k<row;++k)
	{
	for(i=k+1;i<row;++i)//开始LU分解
	{
		A[i][k]/=A[k][k];
		for(j=k+1;j<row;++j)
		A[i][j]-=A[k][j]*A[i][k];
	}
	}
	for(i=0;i<row;++i)
	{
	for(j=0;j<row;++j)//将LU分解结果存入L,U中
	{
		if(i>j)
		L.arr[i][j]=A[i][j],U.arr[i][j]=0;
		else
		U.arr[i][j]=A[i][j],L.arr[i][j]=0;
	}
	L.arr[i][i]=1;//L的对称轴元素置1
	}
	}

	long double Det()//行列式 O(n^3)
	{//(基于PLU分解)
	mat l,u,p;
	int sign=PLU(l,u,p),it=u.row;
	T sum=sign;
	for(int i=0;i<it;++i)
	if(std::isnan(u[i][i]))
	{sum=0;
	break;}
	else
	sum*=u[i][i];
	return sum;
	}

	mat Inv()//逆
	{//基于PLU分解
	mat L,U,P,l(row,row),u(row,row);
	PLU(L,U,P);
	l.set(),u.set();
	for (int i=0;i<row;i++)//求U的逆u
	{
	u[i][i]=1/U[i][i];
	for(int k=i-1;k>=0;k--)
	{
		T s=0;
		for(int j=k+1;j<=i;j++)
		s=s+U[k][j]*u[j][i];
		u[k][i]=-s/U[k][k];
	}
	}
	for(int i=0;i<row;i++)//求L的逆l
	{
	l[i][i]=1;
	for(int k=i+1;k<row;k++)
	{
	for(int j=i;j<=k-1;j++)
	l[k][i]=l[k][i]-L[k][j]*l[j][i];
	}
	}
	L=u*l;
	P=L*P;
	return P;
	}
	
	mat lusolve(mat&L,mat&U)//LU分解求解方程组Ax=b
	{
	int i,k;
	mat x(1,col),y(1,col),b(*this);
	for(i=0;i<col;++i)//解Ly=b
	{
		y[0][i]=b[0][i];
		for(k=i+1;k<col;++k)
		b[0][k]-=L[k][i]*y[0][i];
	}
	for(i=col-1;i>=0;--i)//解Ux=y
	{
		x[0][i]=y[0][i]/U[i][i];
		for(k=i-1;k>=0;--k)
		y[0][k]-=U[k][i]*x[0][i];
	}
	return x;
	}

	mat lusolve(mat&L,mat&U,mat&P)//PLU分解求解方程组Ax=b
	{
	int i,j,k;
	mat Pb(1,col),x(1,col),y(1,col);
	for(i=0;i<col;++i)//求P*b
	for(j=0;j<col;++j)
	Pb[0][i]+=arr[0][j]*P[i][j];
	for(i=0;i<col;++i)//解Ly=Pb
	{
		y[0][i]=Pb[0][i];
		for(k=i+1;k<col;++k)
		Pb[0][k]-=L[k][i]*y[0][i];
	}
	for(i=col-1;i>=0;--i)//解Ux=y
	{
		x[0][i]=y[0][i]/U[i][i];
		for(k=i-1;k>=0;--k)
		y[0][k]-=U[k][i]*x[0][i];
	}
	return x;
	}
	
	bool Cholesky(mat&L)//A=LL'
	{//将任意n阶对称正定矩阵A分解成下三角矩阵L
		if(row!=col)
		{std::cerr<<"\nCholesky\nThis matrix is not a square.\n";
		return 0;
		}
		L.re(row,row);
		for(int k=0;k<row;k++)  
    	{  
        T sum=0;  
		//L[k][k]=sqrt(A[k][k]-sum(L[k][i]^2)) {0<=i<k}
        for(int i=0;i<k;i++)  
        sum+=L[k][i]*L[k][i];  
        sum=arr[k][k]-sum;  
        L[k][k]=sqrt(sum>0?sum:0);
		//L[i][k]=(A[i][k]-sum(L[i][j]*L[k][j]))/L[k][k] {0<=j<k}
        for(int i=k+1;i<row;i++)  
        {  
            sum=0;  
            for(int j=0;j<k;j++)  
            sum+=L[i][j]*L[k][j];  
            L[i][k]=(arr[i][k]-sum)/L[k][k];  
        }
        for(int j=0;j<k;j++)  
        L[j][k]=0;  
    	}//如果对称轴上有0则矩阵不正定
		return 1;
	}

	bool Cholesky(mat&L,mat&D)//A=LDL'
	{//精度更高的Cholesky分解
	//D是对角都是正数的对角阵(如果不全是正的说明不正定)
		if(row!=col)
		{std::cerr<<"\nCholesky\nThis matrix is not a square.\n";
		return 0;
		}
		mat A=(*this);
		for(int k=0;k<row;k++)  
    	{  
        for(int i=0;i<k;i++)  
        A[k][k]-=A[i][i]*A[k][i]*A[k][i];  
        for(int j=k+1;j<row;j++)   
		{
		for(int i=0;i<k;i++)  
		A[j][k]-=A[j][i]*A[i][i]*A[k][i];  
        A[j][k]/=A[k][k];  
		}
    	}  
    	L.re(row,row),D.re(row,row); 
		L.set(0),D.set(0); 
    	for(int i=0;i<row;i++)  
    	{ 
        D[i][i]=A[i][i];  
        L[i][i]=1;  
    	}  
    	for(int i=0;i<row;i++)   
        for(int j=0;j<i;j++)  
        L[i][j]=A[i][j];  
		return 1;
	}

	void max(int&iloc,int&jloc,T&m)//获得最大元素及其索引
	{
		m=arr[0][0];
		iloc=0,jloc=0;
		for(int i=0;i<row;++i)
		for(int j=0;j<col;++j)
		if(m<arr[i][j])
		m=arr[i][j],iloc=i,jloc=j;
	}
	
	void min(int&iloc,int&jloc,T&m)//获得最小元素及其索引
	{
		m=arr[0][0];
		iloc=0,jloc=0;
		for(int i=0;i<row;++i)
		for(int j=0;j<col;++j)
		if(m>arr[i][j])
		m=arr[i][j],iloc=i,jloc=j;
	}

	T sum(bool axis=0,int i=0)//求第i行(列)之和
	{
		T sum=0;
		if(axis==0)
		{
			for(int j=0;j<col;++j)
			sum+=arr[i][j];
			return sum;
		}
		for(int j=0;j<row;++j)
		sum+=arr[j][i];
		return sum;
	}

	T mean(bool axis=0,int i=0)//求第i行(列)之和
	{
		T sum=0;
		if(axis==0)
		{
			for(int j=0;j<col;++j)
			sum+=arr[i][j];
			return sum/col;
		}
		for(int j=0;j<row;++j)
		sum+=arr[j][i];
		return sum/row;
	}

	T inf_norm()//向量无穷范数
	{
		T m=fabs(arr[0][0]),m1;
		for(int i=0;i<row;++i)
		for(int j=0;j<col;++j)
		{
		m1=fabs(arr[i][j]);
		if(m<m1)
		m=m1;
		}
		return m;
	}
	
	T norm(int p=2)//向量p范数
	{
		T sum=0;
		for(int i=0;i<row;++i)
		for(int j=0;j<col;++j)
		sum+=fabs(pow(arr[i][j],p));
		sum=pow(sum,1/(long double)(p));
		return sum;
	} 
	
	T dnorm(int p=2)//对角元p范数
	{
		T sum=0;
		for(int i=0;i<row;++i)
		sum+=fabs(pow(arr[i][i],p));
		sum=pow(sum,1/(long double)(p));
		return sum;
	} 
	
	T cond()//条件数
	{//可能会是无穷大
		mat v=singular();
		return v[0][0]/v[0][v.col-1];
	}
	
	int rank()//秩
	{
		mat v=singular();
		int ti=0;
		for(int i=0;i<v.col;++i)
		if(v[0][i]>1e-14)
		ti++;
		return ti;
	}
	
	mat vcol(int i0)//列向量
	{//第i0列
		if(i0>=col)
		{
			std::cerr<<"\nIndex out of range.\n";
			exit(1);
		}
		mat c(row,1);
		for(int i=0;i<row;++i)
		c.arr[i][0]=arr[i][i0];
		return c;
	}

	mat vrow(int i0)//行向量
	{
		if(i0>=row)
		{
			std::cerr<<"\nIndex out of range.\n";
			exit(1);
		}
		mat c(1,col);
		for(int i=0;i<col;++i)
		c.arr[0][i]=arr[i0][i];
		return c;
	}

	mat Schmidt()//施密特正交化
	{//已归一化
		mat S(*this);
		ld sum;
		for(int k=0;k<col;++k)
		{
		for(int t=0;t<k;++t)
		{
		sum=0;
		for(int j=0;j<row;++j)
		sum+=arr[j][k]*S[j][t];
		for(int i=0;i<row;++i)
		S[i][k]-=sum*S[i][t];
		}
		sum=0;
		for(int i=0;i<row;++i)
		sum+=pow(S[i][k],2);
		sum=sqrt(sum);
		for(int i=0;i<row;++i)
		S[i][k]/=sum;
		}//归一化后也就是QR分解的Q
		return S;
	}	

	void normalize()//归一化
	{//对当前对象归一化
		ld sum=0;
		if(row==1)
		{
		for(int i=0;i<col;++i)
		sum+=pow(arr[0][i],2);
		sum=sqrt(sum);
		for(int i=0;i<col;++i)
		arr[0][i]=arr[0][i]/sum;
		}
		else 
		{
		for(int k=0;k<col;++k)
		{sum=0;
		for(int i=0;i<row;++i)
		sum+=pow(arr[i][k],2);
		sum=sqrt(sum);
		for(int i=0;i<row;++i)
		arr[i][k]=arr[i][k]/sum;
		}
		}
	}

	mat jacobi(long double acc=0.01)//jacobi正交变换
	{//精度(过高会陷入死循环)
	T theta,c,s,max;
	mat eigen(1,row),X(row,row),ma;
	ma=(*this);
	for(int i=0,ti=ma.row;i<ti;++i)//I
	X[i][i]=1;
	while(true)
	{
	max=fabs(ma[1][0]);
	int iloc=1,jloc=0;
	for(int i=2,ti=ma.row;i<ti;++i)
	for(int j=0;j<i;++j)
	if(max<fabs(ma[i][j]))
	max=fabs(ma[i][j]),iloc=i,jloc=j;//选主元
	if(max<acc)//终止标志
	break;
	if(ma[iloc][iloc]==ma[jloc][jloc])
	theta=pi/4;
	else
	theta=atan(2*ma[iloc][jloc]/(ma[iloc][iloc]-ma[jloc][jloc]))/2;
	c=cos(theta),s=sin(theta);
	X[iloc][iloc]=c,X[jloc][jloc]=c;//构造Jacobi旋转矩阵
	X[iloc][jloc]=s,X[jloc][iloc]=-s;
	ma=X.t()*ma*X;//正交变换
	X[iloc][iloc]=1,X[jloc][jloc]=1;//复原为I,用于下次构造
	X[iloc][jloc]=0,X[jloc][iloc]=0;
	}
	for(int i=0;i<ma.row;++i)
	eigen[0][i]=ma[i][i];//存储特征值
	for(int i=0;i<ma.row;++i)//降序排序
	for(int j=0;j<ma.row-i-1;++j)
	if(eigen[0][j]<eigen[0][j+1])
	c=eigen[0][j],eigen[0][j]=eigen[0][j+1],eigen[0][j+1]=c;
	return eigen;
	}

	long double powerit()//幂法
	{//瑞利商作近似特征值近似
	mat ma;
	ma=(*this);
	mat u,x(1,ma.row);
	long double norm;
	for(int i=0;i<ma.row;++i)
	x[0][i]=1;//初始向量全1,不会变成特征向量
	while(true)
	{
		u=x,norm=fabs(x[0][0]);
		for(int i=1;i<ma.row;++i)
		if(fabs(x[0][i])>norm)//用瑞利商收敛更快
		norm=fabs(x[0][i]);//求无穷范数
		for(int i=0;i<ma.row;++i)
		x[0][i]/=norm;//归一化
		x=(ma*x.t()).t();
		if((x-u).norm(2)<1e-14)
		break;
	}
	return norm;
	}

	mat rqi(ld&lam)//反幂法
	{//利用上次计算出的近似特征值作为lam迭代
	mat u,x(1,row),L,U,P;
	long double norm,ln=-1;
	for(int i=0;i<row;++i)
	x[0][i]=1;//初始向量
	while(fabs(norm-ln)<1e-14)
	{
		u=x,norm=u.norm(2);
		for(int i=0;i<row;++i)
		u[0][i]/=norm;
		lam=(u*((*this)*u.t()))[0][0];//瑞利商
		for(int i=0;i<row;++i)
		arr[i][i]-=lam;
		PLU(L,U,P);
		x=u.lusolve(L,U,P);
		for(int i=0;i<row;++i)
		arr[i][i]+=lam;//复原
		ln=norm;
	}
	u=x,norm=u.norm(2);
	for(int i=0;i<row;++i)
	u[0][i]/=norm;
	lam=(u*((*this)*u.t()))[0][0];
	return u;
	}

	mat rit(ld&lam,ld s)//逆向幂迭代
	{//迭代s附近的主特征值
	mat u,x(1,row),I(row,row),L,U,P;
	long double norm,ln=-1;
	I.set(0);
	for(int i=0;i<row;++i)
	x[0][i]=1,I[i][i]=s;
	((*this)-I).PLU(L,U,P);
	while(true)
	{
		u=x,norm=u.norm(2);
		for(int i=0;i<row;++i)
		u[0][i]/=norm;
		x=u.lusolve(L,U,P);
		lam=(u*x.t())[0][0];
		if((u-x).norm(2)<1e-14)
		break;
		//lam=(x*(*this)*x.t())[0][0];//瑞利商
	}
	lam=1/lam+s;
	u=x,norm=u.norm(2);
	for(int i=0;i<row;++i)
	u[0][i]/=norm;
	return u;
	}

	long double prs(bool op)//最大位移幂法(适合求最大特征值)
	{//最大绝对分量作特征值近似
	mat u,x(1,row),ma;
	ma=(*this);
	long double norm,n=ma.norm(2);
	n=(op==true)?n:-n;
	for(int i=0;i<ma.row;++i)
	x[0][i]=1,ma[i][i]+=n;//初始向量全1
	int i=0;
	while(true)
	{
		i++;
		u=x,norm=fabs(x[0][0]);
		for(int i=1;i<ma.row;++i)
		if(fabs(x[0][i])>norm)
		norm=fabs(x[0][i]);//求无穷范数
		for(int i=0;i<ma.row;++i)
		x[0][i]/=norm;//归一化
		x=(ma*x.t()).t();
		if((x-u).norm(2)<1e-14||i>250000)
		break;
	}
	n=(op==true)?n:-n;
	return (n-norm)*pow(-1,op);
	}

	long double iprs(bool op)//最大位移反幂法(适合求最小特征值)
	{
	mat u,x(1,row),L,U,P,ma;
	ma=(*this);
	long double norm,n=ma.norm(2);
	n=(op==true)?n:-n;
	for(int i=0;i<ma.row;++i)
	ma[i][i]+=n,x[0][i]=1;//初始向量全1,不会变成特征向
	ma.PLU(L,U,P);//LU分解
	int i=0;
	while(true)
	{
		i++;
		u=x,norm=fabs(x[0][0]);
		for(int i=1;i<ma.row;++i)
		if(fabs(x[0][i])>norm)
		norm=fabs(x[0][i]);//求无穷范数
		for(int i=0;i<ma.row;++i)
		x[0][i]/=norm;//归一化
		x=x.lusolve(L,U,P);//求解
		if((x-u).norm(2)<1e-14||i>250000)
		break;
	}
	n=(op==true)?n:-n;
	return (n-1/norm)*pow(-1,op);
	}

	long double invpowerit()//反幂法
	{
	mat ma;
	ma=(*this);
	mat u,x(1,ma.row),L,U,P;
	long double norm;
	for(int i=0;i<ma.row;++i)
	x[0][i]=1;//初始向量全1
	ma.PLU(L,U,P);//LU分解
	int i=0;
	while(true)
	{
		i++;
		u=x,norm=fabs(x[0][0]);
		for(int i=1;i<ma.row;++i)
		if(fabs(x[0][i])>norm)
		norm=fabs(x[0][i]);//求无穷范数
		for(int i=0;i<ma.row;++i)
		x[0][i]/=norm;//归一化
		x=x.lusolve(L,U,P);//求解
		if((x-u).norm(2)<1e-14||i>250000)
		break;
	}
	return 1/norm;
	}

	bool QR(mat&Q,mat&R)//正交分解
	{//QR分解,返回Q矩阵与R,基于householder 
	int r=row,c=col,t=r>c?c:r-1;
	R=(*this);
	long double sum;
	mat x,H(r,r),v;
	Q.re(r,r);
	for(int j=0;j<r;++j)
	for(int k=j+1;k<r;++k)
	Q[j][k]=0,Q[k][j]=0;
	for(int i=0;i<r;++i)//Q化单位矩阵
	Q[i][i]=1;
	for(int i=0;i<t;++i)
	{
		x.re(1,r-i);
		for(int j=0;j<r-i;++j)
		x[0][j]=-R[j+i][i];
		x[0][0]+=x.norm(2);//w-x
		sum=(x*x.t())[0][0];//V'*V
		v=x.t()*x;//V*V'
		v*=-2/sum;//  V*V'/V'V
		for(int j=0;j<r-i;++j)//H=I-P
		v[j][j]+=1;
		for(int j=0;j<i;++j)
		for(int k=j+1;k<r;++k)
		H[j][k]=0,H[k][j]=0;
		for(int j=0;j<i;++j)
		H[j][j]=1;//化单位阵
		for(int j=i;j<r;++j)
		for(int k=i;k<r;++k)
		H[j][k]=v[j-i][k-i];
		Q*=H;//H1*H2*H3...
		R=H*R;
	}
	return 1;
	}

	bool usf_qr(mat&Q0,mat&L)//特征分解(非移动幂迭代)
	{//特征分解:A=Q0*L*(Q0)^-1
	mat lam(1,row),Q(row,row),R,Qb;
	long double Ro=sqrt(row);
	for(int i=0;i<row;++i)
	Q[i][i]=1;//I
	Qb=Q,R=(*this);
	while(true)
	{
		(R*Q).QR(Q,R);
		Qb*=Q;
		if(fabs(Q.dnorm(2)-Ro)<1e-16)
		break;
	}
	Q0=Qb;//特征向量
	L=R*Q;//对角阵(对于一些不可对角化的矩阵,可能不会收敛到对角矩阵)
	return 1;
	}

	void hessen(mat&Q,mat&R)//化上hessenberg矩阵
	{//返回Q,R,Q=Hn*...H2*H1,R为hessenberg矩阵
	//A=Q*R*Q.t()
	int r=row;
	R=*(this);
	long double sum;
	mat x,H(r,r),v;
	Q.re(r,r);
	for(int j=0;j<r;++j)
	for(int k=j+1;k<r;++k)
	Q[j][k]=0,Q[k][j]=0;
	for(int i=0;i<r;++i)//Q化单位矩阵
	Q[i][i]=1;
	for(int i=0;i<r-2;++i)
	{
		x.re(1,r-i-1);
		for(int j=0;j<r-i-1;++j)
		x[0][j]=-R[j+i+1][i];
		x[0][0]+=x.norm(2);//w-x
		sum=(x*x.t())[0][0];//Vt*V
		v=x.t()*x;//V*Vt
		v*=-2/sum;
		for(int j=0;j<r-i-1;++j)//H=I-P
		v[j][j]+=1;
		for(int j=0;j<i+1;++j)
		for(int k=j+1;k<r;++k)
		H[j][k]=0,H[k][j]=0;
		for(int j=0;j<i+1;++j)
		H[j][j]=1;
		for(int j=i+1;j<r;++j)
		for(int k=i+1;k<r;++k)
		H[j][k]=v[j-i-1][k-i-1];
		Q*=H;//H1*H2*H3...
		R=H*R*H;
	}
	}
	
	mat<cld> eig()//计算所有(复)特征值(移动幂迭代),不返回特征向量
	{
	mat A;
	A=(*this);
	mat<cld> lam(1,row);
	mat Q(row,row),R;
	int n=row-1,k;
	long double sum,tem;
	while(n>0)
	{	k=0;
		sum=fabs(A[n][0]);
		for(int i=1;i<n;++i)
		if(sum<fabs(A[n][i]))
		sum=fabs(A[n][i]);
		while(sum>1e-14&&k<500)//避免无法截断
		{k++;
		tem=A[n][n];
		for(int i=0;i<A.row;++i)
		A[i][i]-=tem;
		A.QR(Q,R);
		A=R*Q;
		for(int i=0;i<A.row;++i)
		A[i][i]+=tem;
		}
		if(k<500)
		{lam[0][n]=A[n][n];
		n--;
		A=A.cuts(0,n,0,n);
		}
		else 
		{
			tem=pow((A[n-1][n-1]-A[n][n]),2)+4*A[n][n-1]*A[n-1][n];
			lam[0][n]=(A[n-1][n-1]+A[n][n]+sqrt(cld(tem)))/(long double)2.0;
			lam[0][n-1]=(A[n-1][n-1]+A[n][n]-sqrt(cld(tem)))/(long double)2.0;			
			n-=2;
			A=A.cuts(0,n,0,n);
		}
	}
	if(n>-1)
	lam[0][0]=A[0][0];
	return lam;
	}
	
	mat krylov()//krylov子空间
	{
	int r=row;
	mat x0(1,r),tem,K(r,r),I(r,r),Q,R;
	for(int j=0;j<r;++j)
	for(int k=j+1;k<r;++k)
	I[j][k]=0,I[k][j]=0;
	for(int i=0;i<r;++i)//I化单位矩阵
	I[i][i]=1,x0[0][i]=1;
	
	for(int i=0;i<r;++i)
	{
		tem=I*x0.t();
		for(int j=0;j<r;++j)
		K[j][i]=tem[j][0];
		I*=(*this);
	}
	K.QR(Q,R);
	K=Q.t()*(*this)*Q;//化三对角,与A有相同特征值
	return K;
	}

	mat Lanczos()//Lanczos三对角化
	{
	int r=row;
	long double beta=0,alpha=0;
	mat q0(r,1),q1(r,1),aq(r,1),x(r,1),u,K(r,r);
	for(int j=0;j<r;++j)
	q0[j][0]=0,x[j][0]=1,q1[j][0]=1/sqrt(r);
	for(int j=0;j<r;++j)
	for(int k=j+1;k<r;++k)
	K[j][k]=0,K[k][j]=0;
	for(int i=0;i<r;++i)//K化单位矩阵
	K[i][i]=1;
	for(int i0=0;i0<r;++i0)
	{//Aqi=βi−1qi−1+aiqi+βiqi+1
		u=(*this)*q1;
		alpha=(q1.t()*u)[0][0];
		for(int i=0;i<r;++i)
		q0[i][0]*=beta,aq[i][0]=q1[i][0]*alpha;
		u=u-q0-aq;
		beta=u.norm(2);
		q0=q1;
		for(int i=0;i<r;++i)
		q1[i][0]=u[i][0]/beta;
		if(beta==0)
		break;
		K[i0][i0]=alpha;
		if(i0!=r-1)
		K[i0][i0+1]=beta,K[i0+1][i0]=beta;
	}
	return K;//lanczos三对角,与A有相同特征值
	}

	bool svd(mat&U,mat&S,mat&V)//奇异值分解(基于特征分解)
	{//计算A^T*A,误差较大
		mat A=(*this).t()*(*this),D,Q;
		//int min=row>col?col:row;
		A.hessen(V,A);//化上hessenberg矩阵
		A.usf_qr(Q,D);
		V*=Q;
		S.re(row,col),U.re(row,row);
		for(int i=0;i<row;++i)
		for(int j=0;j<col;++j)
		S[i][j]=(i==j)?sqrt(fabs(D[i][i])):0;
		for(int i=0;i<row;++i)
		for(int j=0;j<row;++j)
		{U[i][j]=0;
		for(int k=0;k<col;++k)
		U[i][j]+=arr[i][k]*V[k][j];
		U[i][j]/=S[j][j];
		}
		return 1;
	}

	void bidiag(mat&U,mat&S,mat&V)//双对角化
	{//A=USV'
	//矩阵双对角化(基于householder)
	int r=row,c=col;
	int t0=r>c?c:r-1,t1=r<c-1?r+1:(r<c?c:r-1);
	int max=t0>t1?t0:t1;
	S=(*this);
	long double sum;
	mat x,H0(r,r),H1(c,c),v;
	V.re(c,c),U.re(r,r);//V化单位阵
	for(int i=0;i<c;++i)
	for(int j=0;j<c;++j)
	V[i][j]=i==j?1:0;
	for(int i=0;i<r;++i)//U化单位阵
	for(int j=0;j<r;++j)
	U[i][j]=i==j?1:0;
	for(int i=0;i<max;++i)
	{
		if(i<t0)
		{//列向量应用householder变换
		x.re(1,r-i);
		for(int j=0;j<r-i;++j)
		x[0][j]=-S[j+i][i];
		x[0][0]+=x.norm(2);//w-x
		sum=(x*x.t())[0][0];//Vt*V
		v=x.t()*x;//V*Vt
		v*=-2/sum;
		for(int j=0;j<r-i;++j)//H=I-P
		v[j][j]+=1;
		for(int j=0;j<i;++j)
		for(int k=j+1;k<r;++k)
		H0[j][k]=0,H0[k][j]=0;
		for(int j=0;j<i;++j)
		H0[j][j]=1;
		for(int j=i;j<r;++j)
		for(int k=i;k<r;++k)
		H0[j][k]=v[j-i][k-i];
		U*=H0;//U1*U2*U3...
		S=H0*S;
		}
		if(i<t1-1)
		{//行向量应用householder变换
		x.re(1,c-i-1);
		for(int j=0;j<c-i-1;++j)
		x[0][j]=-S[i][i+j+1];//行索引会超
		x[0][0]+=x.norm(2);//w-x
		sum=(x*x.t())[0][0];//Vt*V
		v=x.t()*x;//V*Vt
		v*=-2/sum;
		for(int j=0;j<c-i-1;++j)//H=I-P
		v[j][j]+=1;
		for(int j=0;j<i+1;++j)
		for(int k=j+1;k<c;++k)
		H1[j][k]=0,H1[k][j]=0;
		for(int j=0;j<i+1;++j)
		H1[j][j]=1;
		for(int j=i+1;j<c;++j)
		for(int k=i+1;k<c;++k)
		H1[j][k]=v[j-i-1][k-i-1];
		V=H1*V;//...V1*V2*V3
		S*=H1;
		}
	}
	V=V.t();
	}	

	bool SVD(mat&U,mat&S,mat&V)//奇异值分解
	{//Golub-Kahan算法,精度较高
		int p,fl;
		if(row>col)//r>c时会出现问题,应转置
		p=col,fl=1;
		else p=row,fl=0;
		mat A=fl?t():(*this),D,Q,Q1;
		if(fl)A.bidiag(V,D,U);
		else A.bidiag(U,D,V);//A=UDV',此处U,V有对左右奇异向量的良好近似
		(D.t()*D).usf_qr(Q,S);//特征分解D'D=QRQ'
		(D*Q).QR(Q1,S);//DQ=Q1S => D=Q1SQ'
		if(fl)V*=Q1,U*=Q,S=S.t();
		else U*=Q1,V*=Q;//A=UQ1SQ'V' => A=USV'
		
		for(int i=0;i<p;++i)//奇异值符号修正
		if(S[i][i]<0)
		{
		S[i][i]*=-1;
		for(int j=0;j<row;++j)
		V[j][i]*=-1;//修正右奇异向量
		}
		return 1;
	}

	mat singular()//奇异值
	{//返回奇异值
	int r=row,c=col;
	mat S=r>c?t():(*this),D,Q,Q1;
	int t0=r>c?c:r-1,t1=r<c-1?r+1:(r<c?c:r-1);
	int max=t0>t1?t0:t1;
	long double sum;
	mat x,H0(r,r),H1(c,c),v;
	for(int i=0;i<max;++i)
	{
		if(i<t0)
		{//列向量应用householder变换
		x.re(1,r-i);
		for(int j=0;j<r-i;++j)
		x[0][j]=-S[j+i][i];
		x[0][0]+=x.norm(2);//w-x
		sum=(x*x.t())[0][0];//Vt*V
		v=x.t()*x;//V*Vt
		v*=-2/sum;
		for(int j=0;j<r-i;++j)//H=I-P
		v[j][j]+=1;
		for(int j=0;j<i;++j)
		for(int k=j+1;k<r;++k)
		H0[j][k]=0,H0[k][j]=0;
		for(int j=0;j<i;++j)
		H0[j][j]=1;
		for(int j=i;j<r;++j)
		for(int k=i;k<r;++k)
		H0[j][k]=v[j-i][k-i];
		S=H0*S;
		}
		if(i<t1-1)
		{//行向量应用householder变换
		x.re(1,c-i-1);
		for(int j=0;j<c-i-1;++j)
		x[0][j]=-S[i][i+j+1];
		x[0][0]+=x.norm(2);//w-x
		sum=(x*x.t())[0][0];//Vt*V
		v=x.t()*x;//V*Vt
		v*=-2/sum;
		for(int j=0;j<c-i-1;++j)//H=I-P
		v[j][j]+=1;
		for(int j=0;j<i+1;++j)
		for(int k=j+1;k<c;++k)
		H1[j][k]=0,H1[k][j]=0;
		for(int j=0;j<i+1;++j)
		H1[j][j]=1;
		for(int j=i+1;j<c;++j)
		for(int k=i+1;k<c;++k)
		H1[j][k]=v[j-i-1][k-i-1];
		S*=H1;
		}
	}
	(S.t()*S).usf_qr(Q,H0);//特征分解
	S*=Q;
	mat H(r,r);
	int t=r>c?c:r-1;
	for(int i=0;i<t;++i)
	{
		x.re(1,r-i);
		for(int j=0;j<r-i;++j)
		x[0][j]=-S[j+i][i];
		x[0][0]+=x.norm(2);//w-x
		sum=(x*x.t())[0][0];//V'*V
		v=x.t()*x;//V*V'
		v*=-2/sum;//  V*V'/V'V
		for(int j=0;j<r-i;++j)//H=I-P
		v[j][j]+=1;
		for(int j=0;j<i;++j)
		for(int k=j+1;k<r;++k)
		H[j][k]=0,H[k][j]=0;
		for(int j=0;j<i;++j)
		H[j][j]=1;//化单位阵
		for(int j=i;j<r;++j)
		for(int k=i;k<r;++k)
		H[j][k]=v[j-i][k-i];
		S=H*S;
	}
	t=r>c?c:r;
	v.re(1,t);
	for(int i=0;i<t;++i)
	v[0][i]=fabs(S[i][i]);
	return v;
	}	

	mat Pinv(ld tol=1e-16)//伪逆
	{//Moore-Penrose逆
		mat U,S,V;
		SVD(U,S,V);
		int t=row>col?col:row,ti;
		for(ti=0;ti<t;++ti)
		if(S[ti][ti]>tol)
		S[ti][ti]=1/S[ti][ti];//S^-1
		else break;
		ti--;//非0奇异值个数
		S=S.cuts(0,ti,0,ti);
		U=U.cuts(0,row-1,0,ti);
		V=V.cuts(0,row-1,0,ti);
		return V*S*U.t();
	}

	bool sorted(int iloc=0)//Shell排序
	{
		const int Sedgewick[]=//最快分组序列
    	{0,1,5,19,41,109,209,505,929,2161};
		int l=0,h=9,p,ed,d=0;
		T st;
		while(l<=h)//二分查找位置
		{
		p=(h+l)/2;
		if(Sedgewick[p]>col)h=p-1;
		else l=p+1;
		}
		for(;h>=0;d=Sedgewick[h--])
		for(int i=d;i<col;++i)
		{
		for(ed=i-d,st=arr[iloc][i];ed>=0&&st<arr[iloc][ed];ed-=d)
		arr[iloc][ed+d]=arr[iloc][ed];
		arr[iloc][ed+d]=st;
		}
		return 1;
	}

	mat sort(int iloc=0)//Shell排序
	{
		mat vec(1,col);
		for(int i=0;i<col;++i)
		vec[iloc][i]=arr[iloc][i];
		const int Sedgewick[]=//最快分组序列
    	{0,1,5,19,41,109,209,505,929,2161};
		int l=0,h=9,p,ed,d=0;
		T st;
		while(l<=h)//二分查找位置
		{
		p=(h+l)/2;
		if(Sedgewick[p]>col)h=p-1;
		else l=p+1;
		}
		for(;h>=0;d=Sedgewick[h--])
		for(int i=d;i<col;++i)
		{
		for(ed=i-d,st=vec[iloc][i];ed>=0&&st<vec[iloc][ed];ed-=d)
		vec[iloc][ed+d]=vec[iloc][ed];
		vec[iloc][ed+d]=st;
		}
		return vec;
	}

	inline int r()//行数
	{
		return row;
	}
	
	inline int c()//列数
	{
		return col;
	}

	std::pair<int,int>size()//矩阵大小
	{
		return{row,col};
	}
	
	inline int max_row()//已分配的行数
	{
		return max_r;
	}

	inline int max_col()//已分配的列数
	{
		return max_c;
	}

	std::pair<int,int>max_size()//矩阵最大大小
	{
		return{max_r,max_c};
	}

};

mat<long double>hilb(int n)//希尔伯特矩阵
{//一个病态矩阵
	mat<long double>tem(n,n);
	for(int i=0;i<n;++i)
	for(int j=0;j<n;++j)
	tem[i][j]=1/(long double)(i+j+1);
	return tem;
}

template<class T=ld>
mat<T>eye(int m,int n)//单位矩阵
{
	mat<long double>tem(m,n);
	for(int i=0;i<m;++i)
	for(int j=0;j<n;++j)
	tem[i][j]=i==j?1:0;
	return tem;
}

template<class T=ld>
mat<T>ones(int m,int n)//0矩阵
{
	mat<long double>tem(m,n);
	for(int i=0;i<m;++i)
	for(int j=0;j<n;++j)
	tem[i][j]=1;
	return tem;
}

template<class T=ld>
mat<T>zeros(int m,int n)//1矩阵
{
	mat<long double>tem(m,n);//calloc自动清零
	return tem;
}

template<class T=ld>
mat<T>rand(int m,int n)//随机矩阵
{
	long long Time=time(NULL),a=16807,m0=~(1<<31);
	long long x0=Time>>3;
	mat<long double>tem(m,n);
	for(int i=0;i<m;++i)
	for(int j=0;j<n;++j)
	{
		x0=(a*x0)%m0;
		tem[i][j]=ld(x0)/m0;
	}
	return tem;
}

#endif