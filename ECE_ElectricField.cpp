/*
Author: Hyeonjae Park
Class: ECE 4122 (A)
Last Date Modified: Sep 27, 2023
Description:

Calculation of electric field & Setting variable for electric field

*/

#include <cmath>

#include "ECE_ElectricField.h"

void ECE_ElectricField::computeFieldAt(double rx, double ry, double rz)
{

	double r; // Set radius

	r = sqrt(pow((rx - x), 2) + pow((ry - y), 2) + pow((rz - z), 2));
	
	Ex = (rx - x) / (pow(r, 3));
	Ey = (ry - y) / (pow(r, 3));
	Ez = (rz - z) / (pow(r, 3));
}

void ECE_ElectricField::getElectricField(double& E_x, double& E_y, double& E_z)
{
	E_x = Ex;
	E_y = Ey;
	E_z = Ez;
}