#include <iostream>
#include "Constant.h"
#include "Functions.h"
#include <string>

using namespace std;

double(*Arr_StateValue)[Ny] = new double[Nx][Ny];
double(*Arr_StateValueNew)[Ny] = new double[Nx][Ny];
State(*Arr_State)[Ny] = new State[Nx][Ny];


int main()
{
	Func_Init();
	int filenameindex = 0;
	for (int i = 0; i < 105; i++)
	{
		cout << "The" << i + 1 << "-th recursion is completed" << endl;
		Func_RecursionMT();
	}
	Func_SavaData("solution.dat", "solution.csv");
}

