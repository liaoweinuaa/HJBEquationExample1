#pragma once
//Computational domain
const double xmin = -1;
const double xmax = 1;
const double ymin = -1;
const double ymax = 1;

//Number of grid points
const int Nx = 201;
const int Ny = 201;



const int savestep = 1;

const int Ntotal = Nx * Ny;
const double gridx = (xmax - xmin) / (Nx - 1);
const double gridy = (ymax - ymin) / (Ny - 1);

//Admissible control inputs
const double umin = -1;
const double umax = 1;

const int Nu = 31;
const double gridu = (umax - umin) / (Nu - 1);


//Time step length
const double dt = 0.02;


//Structure of system state
struct State
{
	double x;
	double y;
};


