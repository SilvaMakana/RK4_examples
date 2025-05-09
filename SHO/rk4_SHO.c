#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* Making a RK4 routine for solving a second order ODE that can be cast into two first order ODE */
/* yʻʻ + k y = 0*/


/* This is the differential equation: yʻ = z */
double diff_equ(double x, double y, double z)
{
	double first_der;
	first_der = z;
	return(first_der);
}

/* This is the other ODE coupled to the first: zʻ = -ky */
double diff_equ1(double k, double x, double y, double z)
{
	double second_der;
	second_der = -k*k * y;
	return(second_der);
}

/* This takes a step size (dx), an inital set of data y(x0)==y0 and z(x0)==z0 */
void rk4(double k, double dx, double x0, int n, double y0, double z0, double* y, double* z)
{
	int i;
	/*double *y, *z;*/
	double *k1, *k2, *k3, *k4;
	double *l1, *l2, *l3, *l4;
	double *x;
	double diff_equ(double, double, double);
	double diff_equ1(double, double, double, double);

	y[0] = y0;
	z[0] = z0;
	x = (double *)malloc(sizeof(double) * n);
	x[0] = x0;

	/*k1 = dx * diff_equ(x0, y0, z0);
	l1 = dx * diff_equ1(x0, y0, z0);

	k2 = dx * diff_equ(x0 + dx / 2, y0 + k1 / 2, z0 + l1 / 2);
	l2 = dx * diff_equ1(x0 + dx / 2, y0 + k1 / 2, z0 + l1 / 2);

	k3 = dx * diff_equ(x0 + dx / 2, y0 + k2 / 2, z0 + l2 / 2);
	l3 = dx * diff_equ1(x0 + dx / 2, y0 + k2 / 2, z0 + l2 / 2);

	k4 = dx * diff_equ(x0 + dx, y0 + k3, z0 + l3);
	l4 = dx * diff_equ1(x0 + dx, y0 + k3, z0 + l3);*/

	/*return (y0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6);*/
	
	k1 = (double *)malloc(sizeof(double) * n);
	k2 = (double *)malloc(sizeof(double) * n);
	k3 = (double *)malloc(sizeof(double) * n);
	k4 = (double *)malloc(sizeof(double) * n);
	l1 = (double *)malloc(sizeof(double) * n);
	l2 = (double *)malloc(sizeof(double) * n);
	l3 = (double *)malloc(sizeof(double) * n);
	l4 = (double *)malloc(sizeof(double) * n);
	for (i = 1; i < n; i++)
		{
			k1[i-1] = dx * diff_equ(x[i-1], y[i-1], z[i-1]);
			l1[i-1] = dx * diff_equ1(k, x[i-1], y[i-1], z[i-1]);

			k2[i-1] = dx * diff_equ(x[i-1] + dx / 2, y[i-1] + k1[i-1] / 2, z[i-1] + l1[i-1] / 2);
			l2[i-1] = dx * diff_equ1(k, x[i-1] + dx / 2, y[i-1] + k1[i-1] / 2, z[i-1] + l1[i-1] / 2);

			k3[i-1] = dx * diff_equ(x[i-1] + dx / 2, y[i-1] + k2[i-1] / 2, z[i-1] + l2[i-1] / 2);
			l3[i-1] = dx * diff_equ1(k, x[i-1] + dx / 2, y[i-1] + k2[i-1] / 2, z[i-1] + l2[i-1] / 2);

			k4[i-1] = dx * diff_equ(x[i-1] + dx, y[i-1] + k3[i-1], z[i-1] + l3[i-1]);
			l4[i-1] = dx * diff_equ1(k, x[i-1] + dx, y[i-1] + k3[i-1], z[i-1] + l3[i-1]);

			y[i] = y[i-1] + (k1[i-1] + 2 * k2[i-1] + 2 * k3[i-1] + k4[i-1]) / 6;
			z[i] = z[i-1] + (l1[i-1] + 2 * l2[i-1] + 2 * l3[i-1] + l4[i-1]) / 6;
			/*printf("%i\t %g\t \n", i, x[i]);*/
		}

	free((char*)x);
	free((char*)k1);
 	free((char*)l1);
 	free((char*)k2);
 	free((char*)l2);
 	free((char*)k3);
 	free((char*)l3);
 	free((char*)k4);
 	free((char*)l4);
}

/* Here we execute it */
int main()
{
	double h=0.1, k, x, x0, x1, y0, y2, z0, *y, *z;
	int i, n;
	printf("Starting point x0, initial position y0, initial derivative z0, spring constant k, and final x point x1: ", x0, y0, z0, k, x1);
	scanf("%lf %lf %lf %lf %lf", &x0, &y0, &z0, &k, &x1);
	n = 1 + (x1 - x0)/h;
	printf("n=%d\n",n);
	y = (double *)malloc(sizeof(double) * n);
	z = (double *)malloc(sizeof(double) * n);
	/*printf("y=%08lX z=%08lX\n",(unsigned long)y,(unsigned long)z);*/

	rk4(k, h, x0, n, y0, z0, y, z);
	/*for (y[0] = 1, i = 1; i < n; i++){y[i] = rk4(h, x0 + h * (i - 1), y[i-1]);}*/
		
 
	printf("x\ty\tdy_dt\ty_soln\trel. err. of y[i]/y_soln-1\n------------\n");
	for (i = 0; i < n; i++) {
		x = x0 + h * i;
		/*y2 = pow(x * x / 4 + 1, 2);*/
		// y2 = exp(-3*x) + 2*exp(2*x);
        y2 = cos(k * x);
		printf("%g\t%g\t%g\t%g\t%g\n", x, y[i], z[i], y2, y[i]/y2 - 1);
	}
 
 	free((char*)y);
 	free((char*)z);
	return 0;
}
