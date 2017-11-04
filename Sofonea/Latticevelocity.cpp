#include "Latticevelocity.h"

int* setLatticeVelocityX()
{
	int *v_x = new int[BASIS];
	v_x[0] = 0;  v_x[1] = 1;  v_x[2] = 0;
	v_x[3] = -1; v_x[4] = 0;  v_x[5] = 1;
	v_x[6] = -1; v_x[7] = -1; v_x[8] = 1;

	return v_x;

}

int* setLatticeVelocityY()
{
	int* v_y = new int [BASIS];
	v_y[0] = 0; v_y[1] = 0;  v_y[2] = 1;
	v_y[3] = 0; v_y[4] = -1; v_y[5] = 1;
	v_y[6] = 1; v_y[7] = -1; v_y[8] = -1;
	return v_y;
}