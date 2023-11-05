#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;
const int Q = 9;
const int NX = 128;
const int NY = 128;
const double U = 0.1;
const double rho0 = 1.0;

int e[Q][2] = {{0, 0}, {0, 1}, {1, 0}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
double w[Q] = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};
double rho[NX + 1][NY + 1], u[NX + 1][NY + 1][2], u0[NX + 1][NY + 1][2], f[NX + 1][NY + 1][Q], F[NX + 1][NY + 1][Q];
int i, j, k, dx, dy, Lx, Ly, dt, c, ip, jp, n;
double re, niu, tan_f, error;

double feq(int i, int j, int k)
{
	double feq, f;
	double eu = e[k][0] * u[i][j][0] + e[k][1] * u[i][j][1];
	double eu2 = eu * eu;
	double u2 = u[i][j][0] * u[i][j][0] + u[i][j][1] * u[i][j][1];
	f = 1 + 3 * eu + 4.5 * eu2 - 1.5 * u2;
	feq = w[k] * rho[i][j] * f;
	return feq;
}

void init()
{
	dx = 1;
	dy = 1;
	Lx = dx * NX;
	Ly = dy * NY;
	dt = dx;
	c = dt / dx;
	re = 1000;
	niu = U * Lx / re;
	tan_f = 3 * niu + 0.5;

	for (i = 0; i <= NX; i++)

		for (j = 0; j <= NY; j++)
		{
			rho[i][j] = rho0;
			u[i][j][0] = 0;
			u[i][j][1] = 0;
			u[i][NY][0] = U;
			for (k = 0; k < Q; k++)
				f[i][j][k] = feq(i, j, k);
		}
}

void evolution()
{
	for (i = 1; i < NX; i++)
		for (j = 1; j < NY; j++)
			for (k = 0; k < Q; k++)
			{
				ip = i - e[k][0];
				jp = j - e[k][1];
				F[i][j][k] = f[ip][jp][k] + 1.0 / tan_f * (feq(ip, jp, k) - f[ip][jp][k]);
			}

	for (i = 1; i < NX; i++) // 宏观量

		for (j = 1; j < NY; j++)
		{
			u0[i][j][0] = u[i][j][0];
			u0[i][j][1] = u[i][j][1];
			rho[i][j] = 0;
			u[i][j][0] = 0;
			u[i][j][1] = 0;
			for (k = 0; k < Q; k++)
			{
				f[i][j][k] = F[i][j][k];
				rho[i][j] += f[i][j][k];
				u[i][j][0] += e[k][0] * f[i][j][k];
				u[i][j][1] += e[k][1] * f[i][j][k];
			}
			u[i][j][0] /= rho[i][j];
			u[i][j][1] /= rho[i][j];
		}
	for (j = 1; j < NY; j++) // bianjie
	{
		for (k = 0; k < Q; k++)
		{
			rho[NX][j] = rho[NX - 1][j];
			f[NX][j][k] = feq(NX, j, k) + f[NX - 1][j][k] - feq(NX - 1, j, k);
			rho[0][j] = rho[1][j];
			f[0][j][k] = feq(0, j, k) + f[1][j][k] - feq(1, j, k);
		}
	}
	for (i = 0; i <= NX; i++)
	{
		for (k = 0; k < Q; k++)
		{
			rho[i][NY] = rho[i][NY - 1];
			f[i][NY][k] = feq(i, NY - 1, k) + f[i][NY - 1][k] - feq(i, NY - 1, k);
			rho[i][0] = rho[i][1];
			f[i][0][k] = feq(i, 0, k) + f[i][1][k] - feq(i, 1, k);
		}
	}
}

void output(int m)
{
	ostringstream name;
	name << "cavity" << m << ".dat";
	ofstream out(name.str().c_str());
	out << "Title = \"LBM Lid Driven Flow \"\n"
		<< "VARIABLES = \"x\",\"y\",\"ρ\",\"U\",\"V\"\n"
		<< "ZONE T = \"BOX\", I = " << NX + 1 << ", J = " << NY + 1 << ", F = POINT" << endl;
	for (j = 0; j <= NY; j++)
		for (i = 0; i <= NX; i++)
		{
			out << double(i) / Lx << " "
				<< double(j) / Ly << " " << rho[i][j] << " " << u[i][j][0] << " " << u[i][j][1] << endl;
		}
}

void Error()
{
	double temp1, temp2;
	temp1 = 0;
	temp2 = 0;
	for (i = 1; i < NX; i++)
		for (j = 1; j < NY; j++)
		{
			temp1 += ((u[i][j][0] - u0[i][j][0]) * (u[i][j][0] - u0[i][j][0]) + (u[i][j][1] - u0[i][j][1]) * (u[i][j][1] - u0[i][j][1]));
			temp2 += (u[i][j][0] * u[i][j][0] + u[i][i][1] * u[i][j][1]);
		}
	temp1 = sqrt(temp1);
	temp2 = sqrt(temp2);
	error = temp1 / (temp2 + 1e-30);
}

int main()
{
	using namespace std;
	init();
	for (n = 0;; n++)
	{
		evolution();
		if (n % 100 == 0)
		{
			if (n % 100 == 0)
			{
				Error();
				cout << "The" << n << "th computation result;" << endl
					 << "The u,v of point (NX/2,NY/2)is:" << setprecision(6) << u[NX / 2][NY / 2][0] << "," << u[NX / 2][NY / 2][1] << endl;
				cout << "The max relative error of uv is:" << setiosflags(ios::scientific) << error << endl;
				if (n >= 1000)
				{
					if (n % 1000 == 0)
						output(n);
					if (error < 1.0e-6)
						break;
				}
			}
		}
	}
	return 0;
}