#include <iostream>
#include <cmath> //Ӣ��

using namespace std;
const int Q = 9;
const int L = 100;
const int H = 40;
const int NX = 100;
const int NY = 40;
const double rho0 = 1.0;
const double Fi = 4e-4;
double w[Q] = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
int e[Q][2] = {{0, 0}, {0, 1}, {1, 0}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
double rho[NX + 1][NY + 1], u[NX + 1][NY + 1][2], f[NX + 1][NY + 1][Q], u0[NX + 1][NY + 1][2], F[NX + 1][NY + 1][Q];
double U0 = 0;
int i, j, k, dx, dy, tan_u, ip, jp;
double c, dt, Lx, Ly, niu, Re;

double feq(int i, int j, int k)
{
	double eu = e[k][0] * u[i][j][0] + e[k][1] * u[i][j][1];
	double u2 = u[i][j][0] * u[i][j][0] + u[i][j][1] * u[i][j][1];
	double fe = rho[i][j] * w[k];
	double f = fe * (1.0 + 3 * eu + 4.5 * eu * eu - 1.5 * u2);
	return f;
}

void init()
{
	dx = 1;
	dy = 1;
	c = dx / dt;
	Lx = dx * NX;
	Ly = dy * NY;
	niu = 0.1;
	tan_u = 0.5 + 3 * niu;

	for (i = 0; i < NX; i++)
		for (j = 0; j < NY; j++)
		{
			rho[i][j] = rho0;
			u[i][j][0] = 0;
			u[i][j][1] = 0;
			for (k = 0; k < Q; k++)
				f[i][j][k] = feq(i, j, k);
		}
}

double Fi(int i, int j, int k)
{
}

void collision()
{
	for (i = 1; i < NX; i++)
		for (j = 1; j < NY; j++)
		{

			for (k = 0; k < Q; k++)
			{
				ip = i - e[k][0];
				jp = j - e[k][1];
				f[i][j][k] = f[ip][jp][k] - (f[ip][jp][k] - feq(ip, jp, k)) / tan_u + Fi[i][j][k];
			}
		}
}

void calculatemacroscopic()
{
	for (i = 0; i < NX; i++)
		for (j = 0; j < NY; j++)
		{
			u0[i][j][0] = u[i][j][0];
			u0[i][j][1] = u[i][j][1];
			rho[i][j] = 0;
			u[i][j][0] = 0;
			u[i][j][1] = 0;

			for (k = 0; k < Q; k++)
			{
				F[i][j][k] = f[i][j][k];
				rho[i][j] += F[i][j][k];
				u[i][j][0] += e[k][0] * F[i][j][k];
				u[i][j][1] += e[k][1] * F[i][j][k];
			}
			u[i][j][0] /= rho[i][j];
			u[i][j][1] /= rho[i][j];
		}
}

void border()
{
	for (i = 0; i < NX; i++) // ���±߽�
	{
		rho[i][0] = rho[i][1];
		rho[i][NY] = rho[i][NY - 1];
		f[i][NY][k] = feq(i, NY - 1, k) + f[i][NY - 1][k] - feq(i, NY - 1, k);
		f[i][0][k] = feq(i, 0, k) + f[i][1][k] - feq(i, 1, k);
	}
}
