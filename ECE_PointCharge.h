/*
Author: Hyeonjae Park
Class: ECE 4122 (A)
Last Date Modified: Sep 27, 2023
Description:

Header file for the point charge

*/

#pragma once
#pragma once

class PointCharge
{
    protected:
	    double x; // x-coordinate
	    double y; // y-coordinate
	    double z; // z-coordinate
	    double q; // charge of the point

    public:
	    void setLocation(double, double, double);
	    void setCharge(double);
};