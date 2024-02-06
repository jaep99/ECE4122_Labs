/*
Author: Hyeonjae Park
Class: ECE 4122 (A)
Last Date Modified: Sep 27, 2023
Description:

Header file for the electric field

*/

#pragma once
#include "ECE_PointCharge.h"

class ECE_ElectricField : public PointCharge
{
    protected:
	    double Ex; // Electric field in the x-direction
	    double Ey; // Electric field in the y-direction
	    double Ez; // Electric field in the z-direction

    public:
        void computeFieldAt(double, double, double);
	    void getElectricField(double&, double&, double&);
};