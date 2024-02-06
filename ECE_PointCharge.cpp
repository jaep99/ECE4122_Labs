/*
Author: Hyeonjae Park
Class: ECE 4122 (A)
Last Date Modified: Sep 27, 2023
Description:

Setting variable for point charge

*/

#include "ECE_PointCharge.h"


void PointCharge::setLocation(double x_i, double y_i, double z_i)
{
	x = x_i;
	y = y_i;
	z = z_i;
}

void PointCharge::setCharge(double q_i)
{
	q = q_i;
}