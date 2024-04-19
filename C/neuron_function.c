#include <stdio.h>
#include <stdlib.h> /* for the rand function */
#include <math.h>
#include "nrutil.h"
#define MAX_FILENAME_LENGTH 100

const double PIOVERTWO = 1.57079632679489661923;
const double PI = 3.14159265358979323846;
const double TWOPI = 6.28318520717958647692;

void FullModel(double x, double y[], double dydx[]);
void rk4(double[], double[], int, double, double, double[],
		 void (*derivs)(double, double[], double[]));
void run_HH_model(int seconds, double K_o, double pH, int current);  
double initial_K(double pH);      

/* Global variables */
double alpha_n, alpha_m, alpha_h, beta_n, beta_m, beta_h, m_inf,
	Kin, Naout, E_k, E_na, E_cl, E_ca, Ina, Ik, Icl, Itildepump, Itildeglia, Itildediffusion,
	Cm, g_na, g_naL, g_k, g_kL, g_ahp, g_clL, g_ca, gamma_val, beta, tau, phi,
	rho, glia, epsilon, kbath, I_ext;

int main(void)
{
    int seconds;
    double pH, K_o;
    double pH_vals[] = {7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7};

	printf("Integrate for how many seconds? ");
	scanf("%ld", &seconds);

    system("mkdir results");

    int length = sizeof(pH_vals) / sizeof(pH_vals[0]);

    for (int i = 0; i < length; i++) {
        K_o = initial_K(pH_vals[i]);
        run_HH_model(seconds, K_o, pH_vals[i], 200);
    }

    return 0;
}

double initial_K(double pH){
    double K_o;
    double K_o_initial = 4.0;
    double pH_initial = 7.4;

    K_o = K_o_initial + (pH_initial - pH)*6.0;

    return K_o;
}

void run_HH_model(int seconds, double K_o, double pH, int current){

    int i, j, totalsteps, skip;
	double *y, *derivs, time, timestep;
	FILE *fp;

	I_ext = current;
	
	/*
	y[1]=V, the membrane voltage
	y[2]=n, gating variable
	y[3]=h, gating variable
	y[4]=[K]_o, the extracellular potassium concentration
	y[5]=[Na]_i, the intracellular sodium concentration
	y[6]=[Ca], the calcium concentration
	*/
	y = dvector(1, 6);
	derivs = dvector(1, 6);

	/* set initial conditions */
	y[1] = -50;
	y[2] = 0.08553;
	y[3] = 0.96859;
	y[4] = K_o;
	y[5] = 15.5;
	y[6] = 0.0;

	/* set parameters */

	gamma_val = 0.044494542;
	tau = 1000.0;
	beta = 7.0;

	rho = 1.25;
	glia = 200.0 / 3.0;
	epsilon = 4.0 / 3.0;
	kbath = K_o;

	E_cl = 26.64 * log(6.0 / 130.0);
	E_ca = 120.0;

	Cm = 1.0;
	g_na = 100.0;
	g_naL = 0.0175;
	g_k = 40.0;
	g_kL = 0.05;
	g_ahp = 0.01;
	g_clL = 0.05;
	g_ca = 0.1;
	phi = 3.0;

	/* integration parameters */
	time = 0.0;
	timestep = 0.01;
	skip = 20;
	totalsteps = (int)(1000 * seconds / timestep);

	/* preintegrate to remove transients */
	for (j = 1; j <= 5000; j++)
	{
		FullModel(time, y, derivs);
		rk4(y, derivs, 6, time, timestep, y, FullModel);
		time = time + timestep;
	}

	time = 0.0;

    char filename[MAX_FILENAME_LENGTH];
    snprintf(filename, MAX_FILENAME_LENGTH, "results/datafile_%.1f.dat", pH);

	fp = fopen(filename, "w");
	fprintf(fp, "%lf %lf %lf %lf\n", time / 1000.0, y[1], y[4], y[5]);
	
	for (i = 1; i <= (totalsteps / skip); i++)
	{
		if (i == 30) {
			I_ext = 0;
		}

		for (j = 1; j <= skip; j++)
		{
			FullModel(time, y, derivs);
			rk4(y, derivs, 6, time, timestep, y, FullModel);
			time = time + timestep;
		}

		fprintf(fp, "%lf %lf %lf %lf\n", time / 1000.0, y[1], y[4], y[5]);
		fflush(fp);
	}

	fclose(fp);

	printf("Done!\n");
}
	
/*****/

void FullModel(double x, double y[], double dydx[])
{
	alpha_n = 0.01 * (y[1] + 34.0) / (1.0 - exp(-0.1 * (y[1] + 34.0)));
	beta_n = 0.125 * exp(-(y[1] + 44.0) / 80.0);
	alpha_m = 0.1 * (y[1] + 30.0) / (1.0 - exp(-0.1 * (y[1] + 30.0)));
	beta_m = 4.0 * exp(-(y[1] + 55.0) / 18.0);
	alpha_h = 0.07 * exp(-(y[1] + 44.0) / 20.0);
	beta_h = 1.0 / (1.0 + exp(-0.1 * (y[1] + 14.0)));

	m_inf = alpha_m / (alpha_m + beta_m);

	Kin = 158.0 - y[5];
	Naout = 144.0 - beta * (y[5] - 18.0);

	E_k = 26.64 * log((y[4] / Kin));
	E_na = 26.64 * log((Naout / y[5]));

	Ina = g_na * (m_inf * m_inf * m_inf) * y[3] * (y[1] - E_na) + g_naL * (y[1] - E_na);
	Ik = (g_k * y[2] * y[2] * y[2] * y[2] + g_ahp * y[6] / (1 + y[6])) * (y[1] - E_k) + g_kL * (y[1] - E_k);
	Icl = g_clL * (y[1] - E_cl);

	Itildepump = (rho / (1.0 + exp((25.0 - y[5]) / 3.0))) * (1 / (1 + exp(5.5 - y[4])));
	Itildeglia = (glia / (1.0 + exp((18.0 - y[4]) / 2.5)));
	Itildediffusion = epsilon * (y[4] - kbath);

	/*
	y[1]=V, the membrane voltage
	y[2]=n, gating variable
	y[3]=h, gating variable
	y[4]=[K]_o, the extracellular potassium concentration
	y[5]=[Na]_i, the intracellular sodium concentration
	y[6]=[Ca], the calcium concentration
	*/
	dydx[1] = (1.0 / Cm) * (-Ina - Ik - Icl + I_ext);
	dydx[2] = phi * (alpha_n * (1 - y[2]) - beta_n * y[2]);
	dydx[3] = phi * (alpha_h * (1 - y[3]) - beta_h * y[3]);
	dydx[4] = (1 / tau) * (gamma_val * beta * Ik - 2.0 * beta * Itildepump - Itildeglia - Itildediffusion);
	dydx[5] = (1 / tau) * (-gamma_val * Ina - 3.0 * Itildepump);
	dydx[6] = -0.002 * g_ca * (y[1] - E_ca) / (1 + exp(-(y[1] + 25.0) / 2.5)) - y[6] / 80.0;
}

/*****/

/* note #undef's at end of file */
#define NRANSI
#include "nrutil.h"

void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
		 void (*derivs)(double, double[], double[]))
{
	int i;
	double xh, hh, h6, *dym, *dyt, *yt;

	dym = dvector(1, n);
	dyt = dvector(1, n);
	yt = dvector(1, n);
	hh = h * 0.5;
	h6 = h / 6.0;
	xh = x + hh;
	for (i = 1; i <= n; i++)
		yt[i] = y[i] + hh * dydx[i];
	(*derivs)(xh, yt, dyt);
	for (i = 1; i <= n; i++)
		yt[i] = y[i] + hh * dyt[i];
	(*derivs)(xh, yt, dym);
	for (i = 1; i <= n; i++)
	{
		yt[i] = y[i] + h * dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x + h, yt, dyt);
	for (i = 1; i <= n; i++)
		yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
	free_dvector(yt, 1, n);
	free_dvector(dyt, 1, n);
	free_dvector(dym, 1, n);
}
#undef NRANSI
