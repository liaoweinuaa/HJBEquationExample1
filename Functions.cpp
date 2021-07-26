#include "Constant.h"
#include <iostream>
#include <string>
#include <mutex>
#include <thread>
#include "Functions.h"
#include <fstream>
using namespace std;

//Arraies to hold the values of the solution at the grid points
extern double(*Arr_StateValue)[Ny];
extern double(*Arr_StateValueNew)[Ny];

//The states represented by the grid points
extern State(*Arr_State)[Ny];



//Endpoint cost function
double Func_EndPoint(State s)
{
	return 0;
}

void Func_Init()
{
	memset(Arr_StateValue, 0, sizeof(double) * Nx * Ny);
	memset(Arr_StateValueNew, 0, sizeof(double) * Nx * Ny);
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{

			Arr_State[i][j].x = xmin + gridx * i;
			Arr_State[i][j].y = ymin + gridy * j;
			Arr_StateValue[i][j] = Func_EndPoint(Arr_State[i][j]);

		}
	}


}


//Target set
bool Func_IsTargetSet(State s)
{
	if (s.x <= 0.2 && s.x >= -0.2 && s.y <= 0.2 && s.y >= -0.2)
	{
		return true;
	}
	return false;
}


//Dynamic system, s(t+dt)=s(t)+f(s,u)dt+ 0.5 * df/ds * f(s,u)*dt*dt
State Func_Transition(State s, double u)
{
	if (Func_IsTargetSet(s))
	{
		return s;
	}

	double x = s.x;
	double y = s.y;


	double dotx = y + 1 * x * x;
	double doty = -1 * x + y * y * y + u;


	double ax = dotx * (2 * x) + doty;
	double ay = -dotx + 3 * y * y * doty;


	double newx = x + dotx * dt + 0.5 * ax * dt * dt;
	double newy = y + doty * dt + 0.5 * ay * dt * dt;

	
	if (newx >= xmax)
	{
		newx = xmax - 0.0001;
	}

	if (newx <= xmin)
	{
		newx = xmin + 0.0001;
	}

	if (newy >= ymax)
	{
		newy = ymax - 0.0001;
	}

	if (newy <= ymin)
	{
		newy = ymin + 0.0001;
	}

	return State{
		newx,newy
	};
}



//Bilinear interpolation see https://en.wikipedia.org/wiki/Bilinear_interpolation
double Func_InterPolation(State s0)
{
	if (Func_IsTargetSet(s0))
	{
		return Func_EndPoint(s0);
	}
	double ialpha_double = ((s0.x - xmin) / gridx);
	double itheta_double = ((s0.y - ymin) / gridy);

	int ialpha = int(ialpha_double);
	int itheta = int(itheta_double);

	double v00 = Arr_StateValue[ialpha][itheta];
	double v01 = Arr_StateValue[ialpha][itheta + 1];
	double v10 = Arr_StateValue[ialpha + 1][itheta];
	double v11 = Arr_StateValue[ialpha + 1][itheta + 1];


	double dalpha0 = ialpha_double - ialpha;
	double dalpha1 = 1 - dalpha0;
	double dtheta0 = itheta_double - itheta;
	double dtheta1 = 1 - dtheta0;


	double V00 = dalpha0 * dtheta0;
	double V01 = dalpha0 * dtheta1;
	double V10 = dalpha1 * dtheta0;
	double V11 = dalpha1 * dtheta1;


	return (v00 * V11 + v01 * V10 + v10 * V01 + v11 * V00);
}



//Running cost function
double Func_RunningCost(State s0, State s1, double u)
{
	if (Func_IsTargetSet(s0))
	{
		return 0;
	}

	return 1;
}


//Bellman's principle of optimality
double Func_ValueUnderOptInput(State s0)
{
	double minvalue = 1000000;
	for (int i = 0; i < Nu; i++)
	{
		double u = umin + i * gridu;

		State snew = Func_Transition(s0, u);

		double snewvalue = Func_InterPolation(snew);
		double value = Func_RunningCost(s0, snew, u) * dt + snewvalue;
		if (value < minvalue)
		{
			minvalue = value;
		}

	}
	return minvalue;
}

//Single-threaded recursion
mutex mt;
void Func_RecursionST(int* current_index)
{
	int loc_index;
	while (true)
	{
		mt.lock();
		loc_index = (*current_index);
		(*current_index)++;
		mt.unlock();
		if (loc_index >= Ntotal)
		{
			break;
		}

		int ix = loc_index / (Ny);
		int iy = (loc_index - ix * (Ny));


		Arr_StateValueNew[ix][iy] = Func_ValueUnderOptInput(Arr_State[ix][iy]);
	}
}

//Ten-thread recursion
void Func_RecursionMT()
{
	int current_index = 0;
	thread t0(Func_RecursionST, &current_index);
	thread t1(Func_RecursionST, &current_index);
	thread t2(Func_RecursionST, &current_index);
	thread t3(Func_RecursionST, &current_index);
	thread t4(Func_RecursionST, &current_index);
	thread t5(Func_RecursionST, &current_index);
	thread t6(Func_RecursionST, &current_index);
	thread t7(Func_RecursionST, &current_index);
	thread t8(Func_RecursionST, &current_index);
	thread t9(Func_RecursionST, &current_index);



	t0.join();
	t1.join();
	t2.join();
	t3.join();
	t4.join();
	t5.join();
	t6.join();
	t7.join();
	t8.join();
	t9.join();
	memcpy(Arr_StateValue, Arr_StateValueNew, sizeof(double) * Nx * Ny);
}


//Save data
void Func_SavaData(string filename0, string filename1)
{
	//Save data in binary form
	ofstream ofs(filename0, ios::binary | ios::out);
	ofs.write((const char*)Arr_StateValue, sizeof(double) * Nx * Ny);
	ofs.close();

	//Save data in text form
	ofstream ofs1;
	ofs1.open(filename1);
	for (int i = 0; i < Nx; i += savestep)
	{
		for (int j = 0; j < Ny; j += savestep)
		{
			ofs1 << Arr_StateValue[i][j] << endl;
		}
	}
	ofs1.close();
}