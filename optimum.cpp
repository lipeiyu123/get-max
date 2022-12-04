/*********************************************************************************************************************
 *                            			All Rights Reserved 
 * @file       	  optimum
 * @author     	  LPY（91Mr.Lee）
 * @version    	  1.0.0
 * @date       	  2022/12/04
 ********************************************************************************************************************/


#include <iostream>
#include <math.h>

#define DEBUG_ON 1
#define DEBUG_OFF 0

#define IFDEBUG	DEBUG_OFF

/*
	double Fun(double x , double param[] ,int order);
	x 		自变量
	param[]		函数各次项系数
	order		函数次数
	return		函数值
	...
	dot_Fun、dot2_Fun 分别返回该函数的一阶、二阶导数值
*/

double Fun(double x , double param[] ,int order){
	 
	int length = order + 1;
	double y = 0;
	for(int i = 0 ; i <length ; i++){
		
		y = y + param[i] * pow(x,i);
	}
	return y;
}

double dot_Fun(double x , double param[] ,int order){
	 
	int length = order + 1;
	double y = 0;
	for(int i = 1 ; i <length ; i++){
		
		y = y + i * param[i] * pow(x,i-1);
	}
	return y;
}


double dot2_Fun(double x , double param[] ,int order){
	 
	int length = order + 1;
	double y = 0;
	for(int i = 2 ; i <length ; i++){
		
		y = y + i * (i-1) * param[i] * pow(x,i-2);
	}
	return y;
}

/*
double newton_max(double start_point , 
		  double accuracy ,
 		  double param[],
		  int order,
		  double(*f)(double x , double param[] ,int order) ,
		  double(*dot_f)(double x , double param[] ,int order) )

	start_point     搜索起始点
	accuracy	精度
	param[]		函数各次项系数
	order		函数次数
	return		函数最大值
	f		要搜索的函数
	dot_f		要搜索函数的一阶导数
	...
	secant_method   相同
*/

double newton_max(double start_point , 
		  double accuracy ,
 		  double param[],
		  int order,
		  double(*f)(double x , double param[] ,int order) ,
		  double(*dot_f)(double x , double param[] ,int order) ){

	double xk = start_point + (2 * accuracy);
	double xk_n = start_point ;
	
	#if	IFDEBUG 
		int i = 1;
		std::cout << "newton method" << std::endl;
	#endif

	while(	fabs( xk_n - xk ) > accuracy ){
		
		xk = xk_n;
		xk_n = xk - ( f(xk,param,order) / dot_f(xk,param,order) );		
		


	#if	IFDEBUG
		std::cout << "times " << i << "  x = " << xk_n << std::endl;
		i++;
	#endif				
	}
	
	return xk_n;
}

double secant_method(double start_point , 
		     double accuracy ,
 		     double param[],
		     int order,  
		     double(*dot_f)(double x , double param[] ,int order) ){

	double xk = start_point + (2 * accuracy);
	double xk_n = start_point ;
	double xk_last = 0;
	
	#if	IFDEBUG 
		int i = 1;
		std::cout << "secant method" << std::endl;
	#endif

	while(	fabs( xk_n - xk ) > accuracy ){
		
		xk_last = xk;
		xk = xk_n;
		xk_n = xk - (( xk - xk_last )/( dot_f(xk,param,order) - dot_f(xk_last,param,order) )) * dot_f(xk,param,order);		
		


	#if	IFDEBUG
		std::cout << "times " << i << "  x = " << xk_n << std::endl;
		i++;
	#endif				
	}
	
	return xk_n;

}



int main(){

	double (*f)(double,double * ,int);
	double (*dot_f)(double,double * ,int);
	double (*dot2_f)(double,double * ,int);

	double f1_para[6] = {38.88,26.64,-14.12,-10.66,1.1,1};		//f1 param
	int f1_order = 5;						//f1 order

	double f2_para[6] = {16,16,-8,-8,1,1};				
	int f2_order = 5;

	double temp = 0;

	f = Fun;
	dot_f = dot_Fun;
	dot2_f = dot2_Fun;	

	temp = newton_max(0,0.001,f1_para,f1_order,dot_f,dot2_f);		//newton method f1
	std::cout << "x = " << temp <<std::endl;				//output x
	std::cout << "f1(x) = "<< f(temp , f1_para , f1_order) << std::endl;	//test f(X)

	temp = secant_method(0,0.001,f1_para,f1_order,dot_f);			//secant method f1
	std::cout << "x = " << temp <<std::endl;				//output x
	std::cout << "f1(x) = "<< f(temp , f1_para , f1_order) << std::endl;	//test f(X)


	temp = newton_max(0.5,0.005,f2_para,f2_order,dot_f,dot2_f);		//newton method f2
	std::cout << "x = " << temp <<std::endl;
	std::cout << "f2(x) = "<< f(temp , f2_para , f2_order) << std::endl;

	temp = secant_method(0.5,0.005,f2_para,f2_order,dot_f);			//secant method f2
	std::cout << "x = " << temp <<std::endl;
	std::cout << "f2(x) = "<< f(temp , f2_para , f2_order) << std::endl;


	return 0;
}































