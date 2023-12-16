#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#define Ly 200.0
#define M 201
#define M1 (M+1)
#define N 201
#define N1 (N+1)
#define imin 1
#define imax (N-1)
#define jmin 1
#define jmax (M-1)
#define CFL 0.9
#define Q 9
#define sgn(x) ((x)==0.0?0.0:((x)<0.0?-1.0:1.0))
#define NUM_THREADS 4
const double U0 = 0.02;
const double W = 4.0;
//std::string path = "D:\C++file\diag_GKS\data";
std::string path = "./";
std::string indexC = "_1orderT";
const double sigma = 0.01;
const double indexT = 2.0;
const double lambda = 8.0 / 15.0;//0.02/(3.0*sigma); //lambda ~ U0*dx/beta
//const double lambda=0.02/(3.0*sigma); //lambda ~ U0*dx/beta
//const double lambda=1.0*U0/(3.0*sigma);//3.0*Pe  //Pe=U0*W/(lambda*beta);    lambda ~ U0*dx/beta
double tau, cs2, eta; //mobility and interface function, theta=4*phi*(1-phi)/W
double dx, dy, dt, rdx, rdy, beta, kappa, PI;
double phi[M1][N1], phi0[M1][N1], mu[M1][N1], ux[M1][N1], uy[M1][N1], Fx[M1][N1], Fy[M1][N1]; // Fx/Fy: interface flux in x/y direction
double phi_x[M1][N1], phi_y[M1][N1], phi_xx[M1][N1], phi_xy[M1][N1], phi_yy[M1][N1], D2mu[M1][N1];
double ux_x[M1][N1], ux_y[M1][N1], uy_x[M1][N1], uy_y[M1][N1], D2ux[M1][N1], D2uy[M1][N1];
double e[Q][2], tp[Q];

void gks_ini(void);
void boundary0(void);
void boundary(void);
void flux_x(void);
void flux_y(void);
void Evol(void);
void datadeal(void);
void contour(int);

int main()
{
    int m = 0;
    double t0 = Ly / U0;//period
    double T0 = indexT * t0;//period number
    double t = 0.0;//time


    omp_set_num_threads(NUM_THREADS);
    gks_ini();
    datadeal();
    contour(m);
    printf("mmax=%d \n", (int)(T0 / dt + 0.5));

    while (t < T0)
    {
        boundary();
        Evol();
        m++; t += dt;
        if (m % 1000 == 0)
        {
            printf("m=%d t=%e phicur=%e \n", m, t / t0, phi[M / 2][N / 2]);
        }
        /* if(m % 100000 ==0)
         {
             datadeal();
             contour(m);
         }*/
    }
    datadeal();
    contour(m);
}

void gks_ini()
{
    int i, j, k, ik, jk, I1, I2, J1, J2;
    double x, y, z, d, xc, yc, R, Pe, We;
    double D2Phi;

    PI = 4.0 * atan(1.0);
    dx = dy = Ly / (jmax - jmin + 1); rdx = 1.0 / dx; rdy = 1.0 / dy;
    dt = CFL * dx;
    //W=4.0*dx;//width of interface

  //  W=4*eps*Ly;
  //  beta=0.1;
  //  kappa=beta*W*W/8.0;
  /*
    kappa=0.1*W*W;
    beta=8*kappa/(W*W);
    sigma=sqrt(2*beta*kappa)/6;
  */
  //sigma=0.01;
    beta = 12.0 * sigma / W;
    kappa = 1.5 * sigma * W;

    //  lambda=1.0e-4;
    //  Pe=2000.0;
    //  lambda=Ly*U0/(beta*Pe);

      // lambda=coef*U0*dx/beta;                      //change for parameter 2022.10.20
      // cs2=1.0/3; eta=1.0; tau=lambda/(cs2*eta);
      // Pe=Ly*U0/(beta*lambda);

      //cs2=1.0/3.0; eta=1.0; tau=lambda/(cs2*eta);
    cs2 = 1.0 / 3.0; tau = 1.0; eta = lambda / (cs2 * tau);
    Pe = U0 * Ly / (lambda * beta);
    We = 1.0 * U0 * U0 * Ly / sigma;

    printf("Pe=%e We=%e lambda=%e sigma=%e kappa=%e beta=%e  dt/tau=%e tau=%e\n", Pe, We, lambda, sigma, kappa, beta, dt / tau, tau);

    //init field
    xc = yc = Ly / 2.0; R = 0.4 * Ly;
    J1 = (jmax - jmin + 1) / 4; J2 = 3 * (jmax - jmin + 1) / 4.0; R = Ly / 4.0;
#pragma omp parallel for private(i,x,y,d,z)
    for (j = jmin - 1; j <= jmax + 1; j++) for (i = imin - 1; i <= imax + 1; i++)
    {
        phi[j][i] = 0.0;
        x = (i - imin) * dx; y = (j - jmin) * dx;
        d = (x - xc) * (x - xc) + (y - yc) * (y - yc);  d = sqrt(d);
        z = d - R;
        phi[j][i] = 0.0;
        //Diag
        if ((i < 93) || (i > 109) || (j < 101)) {
            phi[j][i] = 0.5 + 0.5 * tanh(-2 * z / W); //Disk
        }
        ux[j][i] = - U0 * PI * (y / Ly - 0.5);
        uy[j][i] = U0 * PI * (x / Ly - 0.5);

        /*
        //Disk
        if((i<93)||(i>109)||(j<101)){
          phi[j][i]=0.5+0.5*tanh(-2*z/W); //Disk
        }
        ux[j][i]=-U0*PI*(y/Ly-0.5);
        uy[j][i]=U0*PI*(x/Ly-0.5);
        */
    }

    //init D2phi
    for (j = jmin; j <= jmax; j++) for (i = imin; i <= imax; i++)
        phi0[j][i] = phi[j][i];

    e[0][0] = 0.0; e[1][0] = 1.0; e[2][0] = 0.0; e[3][0] = -1.0; e[4][0] = 0.0; e[5][0] = 1.0; e[6][0] = -1.0; e[7][0] = -1.0; e[8][0] = 1.0;
    e[0][1] = 0.0; e[1][1] = 0.0; e[2][1] = 1.0; e[3][1] = 0.0; e[4][1] = -1.0; e[5][1] = 1.0; e[6][1] = 1.0; e[7][1] = -1.0; e[8][1] = -1.0;
    tp[0] = 4.0 / 9; tp[1] = tp[2] = tp[3] = tp[4] = 1.0 / 9; tp[5] = tp[6] = tp[7] = tp[8] = 1.0 / 36;

#pragma omp parallel for private(i,k,D2Phi,jk,ik)
    for (j = jmin; j <= jmax; j++) for (i = imin; i <= imax; i++)
    {
        D2Phi = 0.0;
        for (k = 0; k < Q; k++)
        {
            jk = j + e[k][1]; ik = i + e[k][0];
            D2Phi += tp[k] * (phi[jk][ik] - phi[j][i]);
        }
        //D2Phi -= phi[j][i];
        D2Phi *= (6 * rdx * rdx);
        mu[j][i] = 4 * beta * phi[j][i] * (phi[j][i] - 1.0) * (phi[j][i] - 0.5) - kappa * D2Phi;
    }

    boundary0();

    for (j = jmin; j <= jmax; j++) for (i = imin; i <= imax; i++)
    {
        ux_x[j][i] = ux_y[j][i] = uy_x[j][i] = uy_y[j][i] = 0.0;
        D2ux[j][i] = D2uy[j][i] = 0.0;
        phi_x[j][i] = (phi[j][i + 1] - phi[j][i - 1]) / (2 * dx);
        phi_y[j][i] = (phi[j + 1][i] - phi[j - 1][i]) / (2 * dx);
        phi_xx[j][i] = (phi[j][i + 1] - 2 * phi[j][i] + phi[j][i - 1]) / (dx * dx);
        phi_yy[j][i] = (phi[j + 1][i] - 2 * phi[j][i] + phi[j - 1][i]) / (dx * dx);
        phi_xy[j][i] = (phi[j + 1][i + 1] + phi[j - 1][i - 1] - phi[j + 1][i - 1] - phi[j - 1][i + 1]) / (4 * dx * dx);

        D2mu[j][i] = (mu[j][i + 1] - 4 * mu[j][i] + mu[j][i - 1] + mu[j + 1][i] + mu[j - 1][i]) / (dx * dx);
    }
}

void boundary0()
{
    int i, j, k, i1, i2, j1, j2;
    i1 = imin - 1; i2 = imax + 1;
    j1 = jmin - 1; j2 = jmax + 1;

    for (j = jmin; j <= jmax; j++)
    {
        phi[j][i1] = phi[j][imax];
        phi[j][i2] = phi[j][imin];
        mu[j][i1] = mu[j][imax];
        mu[j][i2] = mu[j][imin];
    }
    for (i = imin - 1; i <= imax + 1; i++)  //the corner cells are included
    {
        phi[j1][i] = phi[jmax][i];
        phi[j2][i] = phi[jmin][i];
        mu[j1][i] = mu[jmax][i];
        mu[j2][i] = mu[jmin][i];
    }
}

void boundary()
{
    int i, j, k, i1, i2, j1, j2;
    i1 = imin - 1; i2 = imax + 1;
    j1 = jmin - 1; j2 = jmax + 1;

#pragma omp parallel for
    for (j = jmin; j <= jmax; j++)
    {
        phi[j][i1] = phi[j][imax];
        phi[j][i2] = phi[j][imin];
        mu[j][i1] = mu[j][imax];
        mu[j][i2] = mu[j][imin];
        phi_x[j][i1] = phi_x[j][imax];
        phi_x[j][i2] = phi_x[j][imin];
        phi_y[j][i1] = phi_y[j][imax];
        phi_y[j][i2] = phi_y[j][imin];
        phi_xx[j][i1] = phi_xx[j][imax];
        phi_xx[j][i2] = phi_xx[j][imin];
        phi_yy[j][i1] = phi_yy[j][imax];
        phi_yy[j][i2] = phi_yy[j][imin];
        phi_xy[j][i1] = phi_xy[j][imax];
        phi_xy[j][i2] = phi_xy[j][imin];
        D2mu[j][i1] = D2mu[j][imax];
        D2mu[j][i2] = D2mu[j][imin];
    }
#pragma omp parallel for
    for (i = imin - 1; i <= imax + 1; i++)  //the corner cells are included
    {
        phi[j1][i] = phi[jmax][i];
        phi[j2][i] = phi[jmin][i];
        mu[j1][i] = mu[jmax][i];
        mu[j2][i] = mu[jmin][i];
        phi_x[j1][i] = phi_x[jmax][i];
        phi_x[j2][i] = phi_x[jmin][i];
        phi_y[j1][i] = phi_y[jmax][i];
        phi_y[j2][i] = phi_y[jmin][i];
        phi_xx[j1][i] = phi_xx[jmax][i];
        phi_xx[j2][i] = phi_xx[jmin][i];
        phi_yy[j1][i] = phi_yy[jmax][i];
        phi_yy[j2][i] = phi_yy[jmin][i];
        phi_xy[j1][i] = phi_xy[jmax][i];
        phi_xy[j2][i] = phi_xy[jmin][i];
        D2mu[j1][i] = D2mu[jmax][i];
        D2mu[j2][i] = D2mu[jmin][i];
    }
}


void flux_x()
{
    int i, j;
    double ux_b, uy_b, phi_bh, phi_b, phi_xb, phi_yb;  // normal to phase interface;
    double mu_x, Ax;
    double ux_xb, ux_yb, uy_xb, uy_yb, D2ux_b, D2mu_b, phi_yyb, phi_xyb;

#pragma omp parallel for private(i,mu_x, Ax,\
        ux_b,uy_b, phi_bh, phi_b, phi_xb, phi_yb,\
        ux_xb, ux_yb,uy_xb,uy_yb,D2ux_b,D2mu_b,\
        phi_yyb,phi_xyb)
    for (j = jmin; j <= jmax; j++) for (i = imin; i <= imax + 1; i++)
    {
        phi_b = 0.5 * (phi[j][i] + phi[j][i - 1]);
        ux_b = 0.5 * (ux[j][i] + ux[j][i - 1]);
        uy_b = 0.5 * (uy[j][i] + uy[j][i - 1]);
        phi_xb = (phi[j][i] - phi[j][i - 1]) * rdx;
        //    phi_yb=0.25*(phi[j+1][i-1]+phi[j+1][i]-phi[j-1][i-1]-phi[j-1][i])*rdy;  //dy at (i-1/2, j)
        phi_yb = 0.5 * (phi_y[j][i] + phi_y[j][i - 1]);  //dy at (i-1/2, j)
        mu_x = (mu[j][i] - mu[j][i - 1]) * rdx;
        D2mu_b = 0.5 * (D2mu[j][i] + D2mu[j][i - 1]); //LW part
        //  phi_bh = phi_b -0.5*dt*(ux_b*phi_xb+uy_b*phi_yb-lambda*D2mu_b); 2nd-CE
        phi_bh = phi_b - 0.5 * dt * (ux_b * phi_xb + uy_b * phi_yb); //1st-CE
        /*
            ux_xb=(ux[j][i]-ux[j][i-1])*rdx;
            ux_yb=0.5*(ux_y[j][i]+ux_y[j][i-1]);
            uy_xb=(uy[j][i]-uy[j][i-1])*rdx;
            uy_yb=0.5*(uy_y[j][i]+uy_y[j][i-1]);
            D2ux_b=0.5*(D2ux[j][i]+D2ux[j][i-1]);
            phi_yyb=0.5*(phi_yy[j][i]+phi_yy[j][i-1]);
            phi_xyb=0.5*(phi_xy[j][i]+phi_xy[j][i-1]);
            Ax=phi_xb*ux_xb + phi_yb*(2*ux_yb - uy_xb) + phi_b*D2ux_b + ux_b*phi_yyb-uy_b*phi_xyb;
            Fx[j][i]=phi_bh*ux_b - lambda*mu_x + lambda*tau*Ax; // 2nd-CE
        */
        Fx[j][i] = phi_bh * ux_b - lambda * mu_x; // 1st-CE
    }
}

void flux_y()
{
    int i, j;
    double ux_b, uy_b, phi_bh, phi_b, phi_xb, phi_yb;  // normal to phase interface;
    double mu_y, Ay;
    double ux_xb, ux_yb, uy_xb, uy_yb, D2uy_b, D2mu_b, phi_xxb, phi_xyb;
#pragma omp parallel for private(j,mu_y, Ay,\
        ux_b,uy_b, phi_bh,phi_b,phi_xb, phi_yb,\
        ux_xb, ux_yb,uy_xb,uy_yb,D2uy_b,D2mu_b,\
        phi_xxb,phi_xyb)
    for (i = imin; i <= imax; i++) for (j = jmin; j <= jmax + 1; j++)
    {
        phi_b = 0.5 * (phi[j][i] + phi[j - 1][i]);
        ux_b = 0.5 * (ux[j][i] + ux[j - 1][i]);
        uy_b = 0.5 * (uy[j][i] + uy[j - 1][i]);
        phi_yb = (phi[j][i] - phi[j - 1][i]) * rdy;
        phi_xb = 0.5 * (phi_x[j][i] + phi_x[j - 1][i]);  //dx at (i, j-1/2)
        //  phi_xb=0.25*(phi[j-1][i+1]+phi[j][i+1]-phi[j-1][i-1]-phi[j][i-1])*rdx;  //dx at (i, j-1/2)
        mu_y = (mu[j][i] - mu[j - 1][i]) * rdy;
        D2mu_b = 0.5 * (D2mu[j][i] + D2mu[j - 1][i]); //LW part
        //  phi_bh = phi_b -0.5*dt*(ux_b*phi_xb+uy_b*phi_yb-lambda*D2mu_b); //2nd-CE
        phi_bh = phi_b - 0.5 * dt * (ux_b * phi_xb + uy_b * phi_yb);  // 1st-CE
        /*
            uy_yb=(uy[j][i]-uy[j-1][i])*rdy;
            uy_xb=0.5*(uy_y[j][i]+uy_x[j-1][i]);
            ux_yb=(ux[j][i]-ux[j-1][i])*rdy;
            ux_xb=0.5*(ux_x[j][i]+ux_x[j-1][i]);
            D2uy_b=0.5*(D2uy[j][i]+D2uy[j-1][i]);
            phi_xxb=0.5*(phi_xx[j][i]+phi_xx[j-1][i]);
            phi_xyb=0.5*(phi_xy[j][i]+phi_xy[j-1][i]);
            Ay=phi_yb*uy_yb + phi_xb*(2*uy_xb - ux_yb) + phi_b*D2uy_b + uy_b*phi_xxb-ux_b*phi_xyb;
            Fy[j][i]=phi_bh*uy_b - lambda*mu_y + lambda*tau*Ay; //2nd-CE
        */
        Fy[j][i] = phi_bh * uy_b - lambda * mu_y; //1st-CE
    }
}

void Evol()
{
    int i, j, k, ik, jk;
    double D2Phi;
    flux_x();
    flux_y();
    //macro data update with flux 
#pragma omp parallel for private(i)
    for (j = jmin; j <= jmax; j++) for (i = imin; i <= imax; i++)
    {
        phi[j][i] += ((Fx[j][i] - Fx[j][i + 1]) * rdx + (Fy[j][i] - Fy[j + 1][i]) * rdy) * dt;
    }
#pragma omp parallel for private(i,D2Phi,k,jk,ik)
    for (j = jmin; j <= jmax; j++) for (i = imin; i <= imax; i++)
    {
        D2Phi = 0.0;
        for (k = 0; k < Q; k++)
        {
            jk = j + e[k][1]; ik = i + e[k][0];
            D2Phi += tp[k] * (phi[jk][ik] - phi[j][i]);
        }
        //D2Phi -= phi[j][i];
        D2Phi *= (6 * rdx * rdx);
        mu[j][i] = 4 * beta * phi[j][i] * (phi[j][i] - 1) * (phi[j][i] - 0.5) - kappa * D2Phi;
    }

    boundary0();

#pragma omp parallel for private(i)
    for (j = jmin; j <= jmax; j++) for (i = imin; i <= imax; i++)
    {
        phi_x[j][i] = (phi[j][i + 1] - phi[j][i - 1]) / (2 * dx);
        phi_y[j][i] = (phi[j + 1][i] - phi[j - 1][i]) / (2 * dx);
        phi_xx[j][i] = (phi[j][i + 1] - 2 * phi[j][i] + phi[j][i - 1]) / (dx * dx);
        phi_yy[j][i] = (phi[j + 1][i] - 2 * phi[j][i] + phi[j - 1][i]) / (dx * dx);
        phi_xy[j][i] = (phi[j + 1][i + 1] + phi[j - 1][i - 1] - phi[j + 1][i - 1] - phi[j - 1][i + 1]) / (4 * dx * dx);

        D2mu[j][i] = (mu[j][i + 1] - 4 * mu[j][i] + mu[j][i - 1] + mu[j + 1][i] + mu[j - 1][i]) / (dx * dx);
    }
}

void  datadeal()
{
    int i, j;
    std::string str1 = path + "phi0" + indexC + std::to_string((int)indexT) + ".dat";
    std::string str2 = path + "phi" + indexC + std::to_string((int)indexT) + ".dat";
    std::string str3 = path + "ux" + indexC + std::to_string((int)indexT) + ".dat";
    std::string str4 = path + "uy" + indexC + std::to_string((int)indexT) + ".dat";
    std::ofstream fout;
    fout.open(str1, std::ios::out);
    for (j = jmin; j <= jmax; j++)
    {
        for (i = imin; i <= imax; i++) fout << phi0[j][i] << " ";
        fout << std::endl;
    }
    fout.close();

    fout.open(str2, std::ios::out);
    for (j = jmin; j <= jmax; j++)
    {
        for (i = imin; i <= imax; i++) fout << phi[j][i] << " ";
        fout << std::endl;
    }
    fout.close();

    fout.open(str3, std::ios::out);
    for (j = jmin; j <= jmax; j++)
    {
        for (i = imin; i <= imax; i++) fout << ux[j][i] << " ";
        fout << std::endl;
    }
    fout.close();

    fout.open(str4, std::ios::out);
    for (j = jmin; j <= jmax; j++)
    {
        for (i = imin; i <= imax; i++) fout << uy[j][i] << " ";
        fout << std::endl;
    }
    fout.close();
}

void contour(int n)
{
    std::string str1 = path + "contour" + indexC + std::to_string((int)indexT) + std::to_string((int)n) + ".dat";
    std::ofstream fout;
    fout.open(str1, std::ios::out);
    fout << "TITLE = Vector\n";
    fout << "Variables = X,Y,PHI,U,V \n";
    fout << "Zone " << "I = " << N - 1 << ",J = " << M - 1 << ",F = POINT" << std::endl;
    for (int j = jmin; j <= jmax; j++)
    {
        for (int i = imin; i <= imax; i++)
        {
            fout << i << "\t" << j << "\t" << phi[j][i] << "\t" << ux[j][i] << "\t" << uy[j][i] << std::endl;;
        }
    }
    fout.close();
}