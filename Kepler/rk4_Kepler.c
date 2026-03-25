#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
/* Making a RK4 routine for solving a second order ODE that can be cast into two first order ODE */
/* Vector equation, so should have 2n - first order ODEs (n = 2 for motion in the plane) */
/* y'' = -G M / y^2*/


/* This is the differential equation: r' = vr and phi' = vphi */
/* FLAG == 0 is for r and FLAG == 1 is for phi */
double first_diff_equ(double t, double r, double phi, double vr, double vphi, int FLAG)
{
	double first_der_r, first_der_phi;
	first_der_r = vr;
    first_der_phi = vphi;

    if (FLAG == 0)
    {
        return(first_der_r);
    }
    
	else if (FLAG == 1)
    {
        return(first_der_phi);
    }

    return NAN;
    
}

/* This is the other ODE coupled to the first: vr' = -GM/r^2 + r^2 * vphi and vphi' = -2vr * vphi / r */
double second_diff_equ(double t, double r, double phi, double vr, double vphi, int FLAG)
{
	double second_der_r, second_der_phi;
    double G = 1.0, M = 1.0; // GM in SI units
	second_der_r = - (G * M) / (r * r) + r * vphi * vphi;
    second_der_phi = - 2 * vr * vphi / r;

    if (FLAG == 0)
    {
        return(second_der_r);
    }
    
	else if (FLAG == 1)
    {
        return(second_der_phi);
    }

    if (r <= 0) 
    {
        return NAN;
    }
    if (r < 1e-12) 
    {
        r = 1e-12; // or stop integration
    }

    return NAN;

}

/* This takes a step size (dx), an inital set of data (r(t0) == r0, phi(t0) == phi0) and (vr(t0) == vr0, vphi(t0) == vphi0) */
void rk4(double dt, double t0, double r0, double phi0, int n, double vr0, double vphi0, double* r, double* phi, double* vr, double* vphi)
{
	int i, FLAG_r = 0, FLAG_phi = 1;
	/*double *y, *z;*/
	double *k1, *k2, *k3, *k4; // For r; FLAG = 0
    double *k11, *k22, *k33, *k44; // For vr
	double *l1, *l2, *l3, *l4; //For phi
    double *l11, *l22, *l33, *l44; //For vphi
	double *t;
	double first_diff_equ(double, double, double, double, double, int);
	double second_diff_equ(double, double, double, double, double, int);

	r[0] = r0;
    phi[0] = phi0;
	vr[0] = vr0;
    vphi[0] = vphi0;
	t = (double *)malloc(sizeof(double) * n);
	t[0] = t0;

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

    k11 = (double *)malloc(sizeof(double) * n);
	k22 = (double *)malloc(sizeof(double) * n);
	k33 = (double *)malloc(sizeof(double) * n);
	k44 = (double *)malloc(sizeof(double) * n);
    
	l1 = (double *)malloc(sizeof(double) * n);
	l2 = (double *)malloc(sizeof(double) * n);
	l3 = (double *)malloc(sizeof(double) * n);
	l4 = (double *)malloc(sizeof(double) * n);

    l11 = (double *)malloc(sizeof(double) * n);
	l22 = (double *)malloc(sizeof(double) * n);
	l33 = (double *)malloc(sizeof(double) * n);
	l44 = (double *)malloc(sizeof(double) * n);
	for (i = 1; i < n; i++)
		{
			k1[i-1] = dt * first_diff_equ(t[i-1], r[i-1], phi[i-1], vr[i-1], vphi[i-1], FLAG_r);
            k11[i-1] = dt * second_diff_equ(t[i-1], r[i-1], phi[i-1], vr[i-1], vphi[i-1], FLAG_r);
			l1[i-1] = dt * first_diff_equ(t[i-1], r[i-1], phi[i-1], vr[i-1], vphi[i-1], FLAG_phi);
            l11[i-1] = dt * second_diff_equ(t[i-1], r[i-1], phi[i-1], vr[i-1], vphi[i-1], FLAG_phi);

			k2[i-1] = dt * first_diff_equ(t[i-1] + dt / 2, r[i-1] + k1[i-1] / 2, phi[i-1] + l1[i-1] / 2, vr[i-1] + k11[i-1] / 2, vphi[i-1] + l11[i-1] / 2, FLAG_r);
			k22[i-1] = dt * second_diff_equ(t[i-1] + dt / 2, r[i-1] + k1[i-1] / 2, phi[i-1] + l1[i-1] / 2, vr[i-1] + k11[i-1] / 2, vphi[i-1] + l11[i-1] / 2, FLAG_r);
            l2[i-1] = dt * first_diff_equ(t[i-1] + dt / 2, r[i-1] + k1[i-1] / 2, phi[i-1] + l1[i-1] / 2, vr[i-1] + k11[i-1] / 2, vphi[i-1] + l11[i-1] / 2, FLAG_phi);
            l22[i-1] = dt * second_diff_equ(t[i-1] + dt / 2, r[i-1] + k1[i-1] / 2, phi[i-1] + l1[i-1] / 2, vr[i-1] + k11[i-1] / 2, vphi[i-1] + l11[i-1] / 2, FLAG_phi);

            k3[i-1] = dt * first_diff_equ(t[i-1] + dt / 2, r[i-1] + k2[i-1] / 2, phi[i-1] + l2[i-1] / 2, vr[i-1] + k22[i-1] / 2, vphi[i-1] + l22[i-1] / 2, FLAG_r);
			k33[i-1] = dt * second_diff_equ(t[i-1] + dt / 2, r[i-1] + k2[i-1] / 2, phi[i-1] + l2[i-1] / 2, vr[i-1] + k22[i-1] / 2, vphi[i-1] + l22[i-1] / 2, FLAG_r);
            l3[i-1] = dt * first_diff_equ(t[i-1] + dt / 2, r[i-1] + k2[i-1] / 2, phi[i-1] + l2[i-1] / 2, vr[i-1] + k22[i-1] / 2, vphi[i-1] + l22[i-1] / 2, FLAG_phi);
            l33[i-1] = dt * second_diff_equ(t[i-1] + dt / 2, r[i-1] + k2[i-1] / 2, phi[i-1] + l2[i-1] / 2, vr[i-1] + k22[i-1] / 2, vphi[i-1] + l22[i-1] / 2, FLAG_phi);

            k4[i-1] = dt * first_diff_equ(t[i-1] + dt, r[i-1] + k3[i-1], phi[i-1] + l3[i-1], vr[i-1] + k33[i-1], vphi[i-1] + l33[i-1], FLAG_r);
			k44[i-1] = dt * second_diff_equ(t[i-1] + dt, r[i-1] + k3[i-1], phi[i-1] + l3[i-1], vr[i-1] + k33[i-1], vphi[i-1] + l33[i-1], FLAG_r);
            l4[i-1] = dt * first_diff_equ(t[i-1] + dt, r[i-1] + k3[i-1], phi[i-1] + l3[i-1], vr[i-1] + k33[i-1], vphi[i-1] + l33[i-1], FLAG_phi);
            l44[i-1] = dt * second_diff_equ(t[i-1] + dt, r[i-1] + k3[i-1], phi[i-1] + l3[i-1], vr[i-1] + k33[i-1], vphi[i-1] + l33[i-1], FLAG_phi);

            t[i] = t[i-1] + dt;
            r[i] = r[i-1] + (k1[i-1] + 2 * k2[i-1] + 2 * k3[i-1] + k4[i-1]) / 6;
            phi[i] = phi[i-1] + (l1[i-1] + 2 * l2[i-1] + 2 * l3[i-1] + l4[i-1]) / 6;
            vr[i] = vr[i-1] + (k11[i-1] + 2 * k22[i-1] + 2 * k33[i-1] + k44[i-1]) / 6;
            vphi[i] = vphi[i-1] + (l11[i-1] + 2 * l22[i-1] + 2 * l33[i-1] + l44[i-1]) / 6;

			/*printf("%i\t %g\t \n", i, x[i]);*/
		}

	free((char*)t);
	free((char*)k1);
    free((char*)k11);
 	free((char*)l1);
    free((char*)l11);
 	free((char*)k2);
    free((char*)k22);
 	free((char*)l2);
    free((char*)l22);
 	free((char*)k3);
    free((char*)k33);
 	free((char*)l3);
    free((char*)l33);
 	free((char*)k4);
    free((char*)k44);
 	free((char*)l4);
    free((char*)l44);
}

/* Here we execute it */
int main()
{
	double h=0.01, t, t0, t1, r0, phi0, vr0, vr02, vphi0, *r, *phi, *vr, *vphi;
    double x, y;
    double L_ini, L_t, Delta_L;
    double E_ini, E_t, Delta_E;
	double eccen_ini, eccen_ini2, eccen_t, eccen_t2, L_ini2, Delta_eccen, eccen_tol = 1.0e-12;
	int i, n, sign_vr;
	char FLAG_INI[10];

	printf("Input velocities or energy: \n");
	scanf("%9s", FLAG_INI); // The acceptable flags are VEL (for specifying velocities or ENERGY (for specifying specific energy and angular momentum))

	if (strcmp(FLAG_INI, "VEL") == 0)
	{
		printf("Starting time t0, initial position (r0, phi0), initial derivative (vr0, vphi0), and final t point t1: \n");
		scanf("%lf %lf %lf %lf %lf %lf", &t0, &r0, &phi0, &vr0, &vphi0, &t1);
    	printf("Read: t0=%g r0=%g phi0=%g vr0=%g vphi0=%g t1=%g\n", t0, r0, phi0, vr0, vphi0, t1);
	}
	
	else if (strcmp(FLAG_INI, "ENERGY") == 0)
	{
		/* ---TODO: This block uses specific energy and eccentricity as inputs, but it is not useful for describing parabolic orbits--- */

		// printf("Starting time t0, initial position (r0, phi0), initial specific energy and eccentricity (E, e), and final t point t1: \n");
		// scanf("%lf %lf %lf %lf %lf %d %lf", &t0, &r0, &phi0, &E_ini, &eccen_ini, &sign_vr, &t1);
    	// printf("Read: t0 = %g r0 = %g phi0 = %g E_ini = %g eccen_ini = %g sign of radial velocity = %d t1 = %g\n", t0, r0, phi0, E_ini, eccen_ini, sign_vr, t1);

		// L_ini2 = (eccen_ini * eccen_ini - 1) / (2 * E_ini);

		// if (L_ini2 < 0) // Check if the square of ang. momentum is positive
		// {
		// 	printf("L_ini^2 = %g \n", L_ini2);
		// 	printf("Negative L_ini^2... BAD combination of E and eccen... try again \n");
		// 	return(1);
		// }
		
		// else
		// {
		// 	L_ini = sqrt(L_ini2);
		// }

		/* --- */

		printf("Starting time t0, initial position (r0, phi0), initial specific energy and ang. momentum (E, L), and final t point t1: \n");
		scanf("%lf %lf %lf %lf %lf %d %lf", &t0, &r0, &phi0, &E_ini, &L_ini, &sign_vr, &t1);
    	printf("Read: t0 = %g r0 = %g phi0 = %g E_ini = %g L_ini = %g sign of radial velocity = %d t1 = %g\n", t0, r0, phi0, E_ini, L_ini, sign_vr, t1);

		vphi0 =  L_ini / (r0 * r0);
		vr02 = 2 * (E_ini + 1 / r0) - vphi0 * r0 * vphi0 * r0;

		if (vr02 < 0) // Check if the square of vr0 is positive
		{
			printf("vr0^2 = %g \n", vr02);
			printf("Negative vr0^2... BAD combination of E and L... try again \n");
			return(1);
		}
		
		else
		{
			if (sign_vr != 1 && sign_vr != -1)
				{
    				printf("sign_vr must be either +1 or -1 \n");
    				return(1);
				}
			vr0 = sign_vr * sqrt(vr02);
		}
		printf("Calculated vr0 and vphi0: t0=%g r0=%g phi0=%g vr0=%g vphi0=%g t1=%g\n", t0, r0, phi0, vr0, vphi0, t1);
	}

	else
	{
		printf("INVALID FLAG! Choose either VEL or ENERGY \n");
		return(2);
	}

	n = 1 + (t1 - t0)/h;
	printf("Number of time steps n = %d\n",n);

	r = (double *)malloc(sizeof(double) * n);
    phi = (double *)malloc(sizeof(double) * n);
	vr = (double *)malloc(sizeof(double) * n);
    vphi = (double *)malloc(sizeof(double) * n);
	/*printf("y=%08lX z=%08lX\n",(unsigned long)y,(unsigned long)z);*/

	rk4(h, t0, r0, phi0, n, vr0, vphi0, r, phi, vr, vphi);
	/*for (y[0] = 1, i = 1; i < n; i++){y[i] = rk4(h, x0 + h * (i - 1), y[i-1]);}*/

   L_ini = vphi0 * r0 * r0;
   E_ini = 0.5 * (vr0 * vr0 + vphi0 * r0 * vphi0 * r0) - 1. / r0;
   eccen_ini2 = 1 + 2 * E_ini * L_ini * L_ini;

   if (fabs(eccen_ini2) < eccen_tol) // For numerical purposes, check if the square of the eccentricity is close to 0
   {
	printf("eccentricity^2 = %g \n", eccen_ini2);
	printf("Enforcing e = 0 \n");
	eccen_ini2 = 0.0;
   }

   eccen_ini = sqrt(eccen_ini2);
   printf("Initial energy, angular momentum, and eccentricity are E_ini = %lg, L_ini = %lg, and eccen_ini = %lg \n", E_ini, L_ini, eccen_ini);
		
 
	printf("t \t r \t phi \t vr \t vphi \t x(t) \t y(t) \t E(t) \t L(t) \t (E_ini - E(t)) / E_ini \t (L_ini - L(t)) / L_ini \t (eccen_ini - eccen_t) / eccen_ini \n------------\n");
	for (i = 0; i < n; i++) {
		t = t0 + h * i;
        L_t = vphi[i] * r[i] * r[i];
        E_t = 0.5 * (vr[i] * vr[i] + vphi[i] * r[i] * vphi[i] * r[i]) - 1. / r[i];
		// eccen_t = sqrt(1 + 2 * E_t * L_t * L_t);
		// Delta_E = (E_ini - E_t) / E_ini;
        // Delta_L = (L_ini - L_t) / L_ini;

		if (fabs(E_ini) < 1e-14) Delta_E = 0.0;
		else Delta_E = (E_ini - E_t) / E_ini;

		if (fabs(L_ini) < 1e-14) Delta_L = 0.0;
		else Delta_L = (L_ini - L_t) / L_ini;
        

		eccen_t2 = 1 + 2 * E_t * L_t * L_t;
		if (fabs(eccen_t2) < eccen_tol) // For numerical purposes, check if the square of the eccentricity is close to 0
   			{
				// printf("eccentricity^2 = %g \n", eccen_ini2);
				// printf("Enforcing e = 0 \n");
				eccen_t2 = 0.0;
   			}
		else if (eccen_t2 < 0.0) // For numerical purposes, check if the square of the eccentricity is negative
			{
    			printf("ERROR: eccen_t^2 became negative at i=%d, t=%g, eccen_t2=%g\n", i, t, eccen_t2);
    			return(3);
			}
   		eccen_t = sqrt(eccen_t2);
		
		if (fabs(eccen_ini) < eccen_tol)
    		Delta_eccen = 0.0;
		else
    		Delta_eccen = (eccen_ini - eccen_t) / eccen_ini;

        x = r[i] * cos(phi[i]);
        y = r[i] * sin(phi[i]);
		// printf("%lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \n", t, r[i], phi[i], vr[i], vphi[i], E_t, L_t, Delta_E, Delta_L);
        printf("%.15e \t %.15e \t %.15e \t %.15e \t %.15e \t %.15e \t %.15e \t %.15e \t %.15e \t %.15e \t %.15e \t %.15e \n", t, r[i], phi[i], vr[i], vphi[i], x, y, E_t, L_t, Delta_E, Delta_L, Delta_eccen);
    }
 
 	free((char*)r);
 	free((char*)phi);
    free((char*)vr);
 	free((char*)vphi);
	return 0;
}
