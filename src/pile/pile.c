/*
Overview:
    This program stacks three cylinders using the control parameter A (see Fig. S1 of the paper)
    and simulates their time evolution under gravity and wall motions.

Output:
    - Time evolution of cylinder–wall/cylinder–cylinder contact forces (for visualization).
    - Final equilibrium configuration written to:
        data/init_deltat_A<val>.txt  (tangential displacements)
        data/init_p_A<val>.txt       (positions/velocities/angles)
      These can be reused as initial conditions for the compression protocol.
    - Discrete snapshots for animation (see 'animation.c' to generate GIFs).
*/

#include <stdio.h>
#include <math.h>
#include <assert.h>

#define PI atan(1.0) * 4 // π

/*=========Material and Model Parameters (dimensionless)==========*/
#define N 3        // Number of cylinders
#define wall 3     // Number of walls
#define r 1.0      // Cylinder radius
#define e 0.3      // Coefficient of restitution
#define mu 0.7     // Friction coefficient between cylinders
#define kn 1e+2    // Normal elastic constant
#define kt 4e+1    // Tangential elastic constant
#define g 1.0      // Gravitational acceleration
#define vfall 0.1  // Initial falling velocity of the released cylinder
#define hight 0.1  // Release hight of the second cylinder (1 + sqrt(3) + hight)
#define A 1.0e-3   // Indentation depth
#define vwall 1e-3 // Velocity of the moving wall

/*============Time-Stepping Parameters=============*/
#define dt 1e-3  // Time step size
#define data 100 // Number of frames for GIF animation
#define tend 10  // Total simulation time

// Set three arrays' element at index n to (x, y, z), respectively.
void init(int n, double *a, double *b, double *c, double x, double y, double z)
{
    a[n] = x;
    b[n] = y;
    c[n] = z;
}

// Initialize an n-element vector 'o' with zeros.
void vec_init0(int n, double o[n])
{
    int i = 0;
    for (i = 0; i < n; i++)
    {
        o[i] = 0;
    }
}

// Initialize an m×n matrix 'o' with zeros.
void matrix_init0(int m, int n, double o[m][n])
{
    int i = 0, j = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            o[i][j] = 0;
        }
    }
}

int main(void)
{
    FILE *fp, *fp1, *fp2, *fp3, *fp4;
    int count = 1, i, j, k, igraph, jgraph, kgraph;
    double t;     // Simulation time
    int it;       // Integer time step index
    double tau;   // First-contact time (when cylinder 2 first makes contact)
    int flag = 0; // Set to 1 after cylinder reaches the cylinders 0 and 1

    /*==============cylinder state==============*/
    double x[N], xnew[N], y[N], ynew[N], theta[N], thetanew[N]; // Positions (x,y) and angles
    double u[N], unew[N], v[N], vnew[N], omega[N], omeganew[N]; // Velocities (vx,vy) and angular vel.

    double fx[N], fy[N]; // Net forces on each cylinder.

    /*===========Cylinder-cylinder interactions============*/
    double l[N][N];                                // Pair distances.
    double lx[N][N], ly[N][N], du[N][N], dv[N][N]; // Relative pos. & vel.
    double delta_nx[N][N], delta_ny[N][N];         // Normal overlap vector (projected).
    double v_nx[N][N], v_ny[N][N];                 // Normal relative velocity.
    double nx[N][N], ny[N][N];                     // Normal unit vector.
    double fnx[N][N], fny[N][N];                   // Normal contact force (components).
    double delta_tx[N][N], delta_ty[N][N];         // Tangential displacement.
    double delta_txnew[N][N], delta_tynew[N][N];   //(Updated) tangential displacement.
    double deltat[N][N];                           //|tangential displacement|
    double vtx[N][N], vty[N][N];                   // Tangential relative velocity.
    double vt[N][N];                               //|tangential relative velocity|
    double tx[N][N], ty[N][N];                     // Tangential unit vector.
    double ftx[N][N], fty[N][N];                   // Tangential frictional force (components).
    double T[N];                                   // Torque about cylinder center.
    double ft[N][N], fn[N][N];                     // tangential and normal contact force.

    /*===============Cylinder-wall interactions==============*/
    double delta_nx_wall[N][wall], delta_ny_wall[N][wall];       // Normal overlap vs walls.
    double vnx_wall[N][wall], vny_wall[N][wall];                 // Normal relative velocity vs walls.
    double fnx_wall[N][wall], fny_wall[N][wall];                 // Normal force vs walls.
    double l_wall[N][wall];                                      // Distance to wall plane.
    double delta_tx_wall[N][wall], delta_ty_wall[N][wall];       // Tangential displacement vs walls.
    double delta_txnew_wall[N][wall], delta_tynew_wall[N][wall]; //(Updated) tangential displacement.
    double deltat_wall[N][wall];                                 //|tangential displacement| vs walls.
    double vtx_wall[N][wall], vty_wall[N][wall];                 // Tangential relative velocity vs walls.
    double vt_wall[N][wall];                                     //|tangential velocity| vs walls.
    double tx_wall[N][wall], ty_wall[N][wall];                   // Tangential unit vector vs walls.
    double ftx_wall[N][wall], fty_wall[N][wall];                 // Tangential frictional force vs walls.
    double T_wall[N][wall];                                      // Torque from walls.
    double ft_wall[N][wall], fn_wall[N][wall];                   //|tangential|, |normal| forces vs walls.

    /*===============Wall geometry, velocity===================*/
    double nx_wall[wall], ny_wall[wall], d[wall]; // Wall normal vector (nx,ny) and offset d.
    double u_wall[wall], v_wall[wall];            // Wall velocities (translational).

    /*=============Material Parameters=================*/
    double m = 1.0;                                                         // Cylinder mass.
    double I = 0.5 * m * r * r;                                             // Moment of inertia of a cylinder.
    double etan = -2 * log(e) * sqrt(m * kn / (PI * PI + log(e) * log(e))); // Normal damping.
    double etat = -2 * log(e) * sqrt(m * kt / (PI * PI + log(e) * log(e))); // Tangential damping.

    /*==========Time stepping and output cadence=======*/
    int itend = tend / dt;                    // Total number of steps
    int countgraph = (int)(tend / dt / data); // Output every "countgraph" steps

    /*=============File destinations==================*/
    char filename[100];

    /*cylinder positions (x,y,theta)*/
    sprintf(filename, "data/p_A%.1e.txt", A);
    fp = fopen(filename, "w");

    /*contact forces (cylinder-cylinder, cylinder-wall)*/
    sprintf(filename, "data/f_A%.1e.txt", A);
    fp1 = fopen(filename, "w");

    /*wall geometry*/
    sprintf(filename, "data/w_A%.1e.txt", A);
    fp2 = fopen(filename, "w");

    /*final cylinder state*/
    sprintf(filename, "data/init_p_A%.1e.txt", A);
    fp3 = fopen(filename, "w");

    /*final tangential displacements*/
    sprintf(filename, "data/init_deltat_A%.1e.txt", A);
    fp4 = fopen(filename, "w");

    /*================Initial conditions===================*/
    t = 0.0;
    it = 0;

    /*cylinder initial states*/
    init(0, x, y, theta, -r, r, 0);
    init(0, u, v, omega, 0, 0, 0);
    init(1, x, y, theta, r, r, 0);
    init(1, u, v, omega, 0, 0, 0);
    init(2, x, y, theta, 0, r * (1 + sqrt(3)) + hight, 0);
    init(2, u, v, omega, 0, 0, 0);

    /*reset tangential state arrays*/
    matrix_init0(N, N, delta_tx);
    matrix_init0(N, N, delta_ty);
    matrix_init0(N, N, delta_txnew);
    matrix_init0(N, N, delta_tynew);

    /*==========Wall setup=========
        k=0: floor (normal vector (0,1)), offset d=0
        k=1: left wall (normal vector (1,0)), offset d=2r
        k=2: right wall (normal vector (-1,0)), offset d=2r*/
    init(0, nx_wall, ny_wall, d, 0.0, 1.0, 0.0);
    init(1, nx_wall, ny_wall, d, 1.0, 0.0, 2 * r);
    init(2, nx_wall, ny_wall, d, -1.0, 0.0, 2 * r);

    /*wall velocities*/
    u_wall[0] = 0.0;
    v_wall[0] = 0.0;
    u_wall[1] = 0.0;
    v_wall[1] = 0.0;
    u_wall[2] = 0.0;
    v_wall[2] = 0.0;

    /*reset wall tangential displacement*/
    matrix_init0(N, wall, delta_tx_wall);
    matrix_init0(N, wall, delta_ty_wall);
    matrix_init0(N, wall, delta_txnew_wall);
    matrix_init0(N, wall, delta_tynew_wall);

    /*==================Time evolution==================*/
    while (t < tend)
    {
        /*==========Compression of walls k = 1, 2===========*/
        if (flag == 0)
        {
            if (t < A / vwall)
            {
                d[1] = 2 * r - vwall * t;
                d[2] = 2 * r - vwall * t;
                u_wall[1] = vwall;
                u_wall[2] = -vwall;
            }
            else if (t >= A / vwall)
            {
                d[1] = 2 * r - A;
                d[2] = 2 * r - A;
                u_wall[1] = 0;
                u_wall[2] = 0;
            }
        }
        else if (flag == 1)
        {
            if (t >= tau + 1)
            {
                d[1] = 2 * r - A + vwall * (t - tau - 1);
                d[2] = 2 * r - A + vwall * (t - tau - 1);
                u_wall[1] = -vwall;
                u_wall[2] = vwall;
            }
        }

        for (i = 0; i < N; i++)
        {
            /*============Compute cylinder-cylinder forces=============*/
            for (j = 0; j < N; j++) // Interaction between cylinder i and j
            {
                if (j == i)
                {
                    continue;
                }

                /*==================Relative distance and relative velocity===================*/
                lx[i][j] = x[i] - x[j];
                ly[i][j] = y[i] - y[j];
                du[i][j] = u[i] - u[j];
                dv[i][j] = v[i] - v[j];
                l[i][j] = sqrt(pow(lx[i][j], 2) + pow(ly[i][j], 2));

                /*=====================Normal unit vector==================*/
                nx[i][j] = lx[i][j] / l[i][j];
                ny[i][j] = ly[i][j] / l[i][j];

                /*===============Normal displacement vector=================*/
                delta_nx[i][j] = (l[i][j] - 2 * r) * nx[i][j];
                delta_ny[i][j] = (l[i][j] - 2 * r) * ny[i][j];

                /*=============Normal relative velocity vector================*/
                v_nx[i][j] = (du[i][j] * nx[i][j] + dv[i][j] * ny[i][j]) * nx[i][j];
                v_ny[i][j] = (du[i][j] * nx[i][j] + dv[i][j] * ny[i][j]) * ny[i][j];

                /*================Tangential relative velocity vector=============*/
                vtx[i][j] = du[i][j] - (du[i][j] * nx[i][j] + dv[i][j] * ny[i][j]) * nx[i][j] + r * (omega[i] + omega[j]) * ny[i][j];
                vty[i][j] = dv[i][j] - (du[i][j] * nx[i][j] + dv[i][j] * ny[i][j]) * ny[i][j] - r * (omega[i] + omega[j]) * nx[i][j];

                /*==============Magnitudes of tangential displacement and velocity=================*/
                deltat[i][j] = sqrt(pow(delta_tx[i][j], 2) + pow(delta_ty[i][j], 2));
                vt[i][j] = sqrt(pow(vtx[i][j], 2) + pow(vty[i][j], 2));

                /*================Tangential unit vector=======================*/
                if (vt[i][j] == 0 && deltat[i][j] > 0)
                {
                    tx[i][j] = delta_tx[i][j] / deltat[i][j];
                    ty[i][j] = delta_ty[i][j] / deltat[i][j];
                }
                else if (vt[i][j] > 0)
                {
                    tx[i][j] = vtx[i][j] / vt[i][j];
                    ty[i][j] = vty[i][j] / vt[i][j];
                }
                else
                {
                    tx[i][j] = 0;
                    ty[i][j] = 0;
                }

                /*==================Normal contact force================*/
                if (l[i][j] - 2 * r < 0) // Contact
                {
                    fnx[i][j] = -kn * delta_nx[i][j] - etan * v_nx[i][j];
                    fny[i][j] = -kn * delta_ny[i][j] - etan * v_ny[i][j];
                }
                else // No contact
                {
                    fnx[i][j] = 0;
                    fny[i][j] = 0;
                }

                /*===============Trial tangential force====================*/
                if (l[i][j] - 2 * r < 0)
                {
                    ftx[i][j] = -kt * delta_tx[i][j] - etat * vtx[i][j];
                    fty[i][j] = -kt * delta_ty[i][j] - etat * vty[i][j];
                }
                else
                {
                    ftx[i][j] = 0;
                    fty[i][j] = 0;
                }

                ft[i][j] = sqrt(pow(ftx[i][j], 2) + pow(fty[i][j], 2));
                fn[i][j] = sqrt(pow(fnx[i][j], 2) + pow(fny[i][j], 2));

                /*============Determine tangential force and update tangential displacement========*/
                if (ft[i][j] >= mu * fn[i][j] && l[i][j] - 2 * r < 0) // Slip
                {
                    ftx[i][j] = -mu * fn[i][j] * tx[i][j];
                    fty[i][j] = -mu * fn[i][j] * ty[i][j];
                    delta_txnew[i][j] = mu * kn * (2 * r - l[i][j]) / kt * tx[i][j];
                    delta_tynew[i][j] = mu * kn * (2 * r - l[i][j]) / kt * ty[i][j];
                }
                else if (ft[i][j] < mu * fn[i][j] && l[i][j] - 2 * r < 0) // No-slip
                {
                    delta_txnew[i][j] = delta_tx[i][j] + vtx[i][j] * dt;
                    delta_tynew[i][j] = delta_ty[i][j] + vty[i][j] * dt;
                }
                else // No contact
                {
                    delta_tx[i][j] = 0;
                    delta_ty[i][j] = 0;
                }

                delta_tx[i][j] = delta_txnew[i][j];
                delta_ty[i][j] = delta_tynew[i][j];
            }

            /*==================Interactions with walls====================*/
            for (k = 0; k < wall; k++)
            {
                /*Distance from cylinder i to wall k*/
                l_wall[i][k] = nx_wall[k] * x[i] + ny_wall[k] * y[i] + d[k];

                /*Normal overlap vector*/
                delta_nx_wall[i][k] = (l_wall[i][k] - r) * nx_wall[k];
                delta_ny_wall[i][k] = (l_wall[i][k] - r) * ny_wall[k];

                /*Normal relative velocity*/
                vnx_wall[i][k] = ((u[i] - u_wall[k]) * nx_wall[k] + (v[i] - v_wall[k]) * ny_wall[k]) * nx_wall[k];
                vny_wall[i][k] = ((u[i] - u_wall[k]) * nx_wall[k] + (v[i] - v_wall[k]) * ny_wall[k]) * ny_wall[k];

                /*Tangential relative velocity*/
                vtx_wall[i][k] = (u[i] - u_wall[k]) - vnx_wall[i][k] + r * omega[i] * ny_wall[k];
                vty_wall[i][k] = (v[i] - v_wall[k]) - vny_wall[i][k] - r * omega[i] * nx_wall[k];

                /*Tangential displacement magnitude*/
                deltat_wall[i][k] = sqrt(pow(delta_tx_wall[i][k], 2) + pow(delta_ty_wall[i][k], 2));
                vt_wall[i][k] = sqrt(pow(vtx_wall[i][k], 2) + pow(vty_wall[i][k], 2));

                /*Tangential unit vector*/
                if (vt_wall[i][k] == 0 && deltat_wall[i][k] > 0)
                {
                    tx_wall[i][k] = delta_tx_wall[i][k] / deltat_wall[i][k];
                    ty_wall[i][k] = delta_ty_wall[i][k] / deltat_wall[i][k];
                }
                else if (vt_wall[i][k] > 0)
                {
                    tx_wall[i][k] = vtx_wall[i][k] / vt_wall[i][k];
                    ty_wall[i][k] = vty_wall[i][k] / vt_wall[i][k];
                }
                else
                {
                    tx_wall[i][k] = 0;
                    ty_wall[i][k] = 0;
                }

                /*Normal contact force*/
                if (l_wall[i][k] - r < 0)
                {
                    fnx_wall[i][k] = -kn * delta_nx_wall[i][k] - etan * vnx_wall[i][k];
                    fny_wall[i][k] = -kn * delta_ny_wall[i][k] - etan * vny_wall[i][k];
                }
                else
                {
                    fnx_wall[i][k] = 0;
                    fny_wall[i][k] = 0;
                }

                /*Trial tangentail frictional force*/
                if (l_wall[i][k] - r < 0)
                {
                    ftx_wall[i][k] = -kt * delta_tx_wall[i][k] - etat * vtx_wall[i][k];
                    fty_wall[i][k] = -kt * delta_ty_wall[i][k] - etat * vty_wall[i][k];
                }
                else
                {
                    ftx_wall[i][k] = 0;
                    fty_wall[i][k] = 0;
                }

                fn_wall[i][k] = sqrt(pow(fnx_wall[i][k], 2) + pow(fny_wall[i][k], 2));
                ft_wall[i][k] = sqrt(pow(ftx_wall[i][k], 2) + pow(fty_wall[i][k], 2));

                /*Determine tangential force and update tangential displacement*/
                if (ft_wall[i][k] >= mu * fn_wall[i][k] && l_wall[i][k] - r < 0) // Slip
                {
                    ftx_wall[i][k] = -mu * fn_wall[i][k] * tx_wall[i][k];
                    fty_wall[i][k] = -mu * fn_wall[i][k] * ty_wall[i][k];
                    delta_txnew_wall[i][k] = mu * kn * (r - l_wall[i][k]) / kt * tx_wall[i][k];
                    delta_tynew_wall[i][k] = mu * kn * (r - l_wall[i][k]) / kt * ty_wall[i][k];
                }
                else if (ft_wall[i][k] < mu * fn_wall[i][k] && l_wall[i][k] - r < 0) // No slip
                {
                    delta_txnew_wall[i][k] = delta_tx_wall[i][k] + vtx_wall[i][k] * dt;
                    delta_tynew_wall[i][k] = delta_ty_wall[i][k] + vty_wall[i][k] * dt;
                }
                else
                {
                    delta_txnew_wall[i][k] = 0;
                    delta_tynew_wall[i][k] = 0;
                }

                delta_tx_wall[i][k] = delta_txnew_wall[i][k];
                delta_ty_wall[i][k] = delta_tynew_wall[i][k];
            }
        }

        /*===========Update positions and velocities, record to file==============*/
        for (i = 0; i < N; i++)
        {
            fx[i] = 0.0;
            fy[i] = 0.0;
            T[i] = 0.0;
        }

        for (i = 0; i < N; i++)
        {

            /*Add cylinder-cylinder interactions to the resultant force/torque*/
            for (j = 0; j < N; j++)
            {
                if (i == j)
                {
                    continue;
                }
                fx[i] += fnx[i][j] + ftx[i][j];
                fy[i] += fny[i][j] + fty[i][j];
                T[i] += -r * (nx[i][j] * fty[i][j] - ny[i][j] * ftx[i][j]);
            }

            /*Add cylinder-wall interactions to the resultant*/
            for (k = 0; k < wall; k++)
            {
                if (k == 0)
                {
                    fx[i] += fnx_wall[i][k] + ftx_wall[i][k];
                    fy[i] += fny_wall[i][k] + fty_wall[i][k];
                    T[i] += -r * (nx_wall[k] * fty_wall[i][k] - ny_wall[k] * ftx_wall[i][k]);
                }
                else /*Compressing walls k = 1, 2 are frictionless*/
                {
                    fx[i] += fnx_wall[i][k];
                    fy[i] += fny_wall[i][k];
                }
            }

            /*Add gravity*/
            fy[i] += -g;

            /*=========================Update positions and velocities======================*/
            /*=========Gently place cylinder 2=========*/
            if (flag == 0)
            {
                if (t < A / vwall)
                {
                    fy[2] = 0;
                }
                else if (t > A / vwall && l[0][2] > 2 * r)
                {
                    fy[2] = 0;
                    v[2] = -vfall;
                }
                else
                {
                    flag = 1;
                    tau = t;
                }
            }

            /*Velocity update (leapfrog)*/
            unew[i] = u[i] + fx[i] * dt;
            vnew[i] = v[i] + fy[i] * dt;
            omeganew[i] = omega[i] + T[i] * dt / I;

            /*Position update (leapfrog)*/
            xnew[i] = x[i] + unew[i] * dt;
            ynew[i] = y[i] + vnew[i] * dt;
            thetanew[i] = theta[i] + omeganew[i] * dt;

            x[i] = xnew[i];
            y[i] = ynew[i];
            u[i] = unew[i];
            v[i] = vnew[i];
            theta[i] = thetanew[i];
            omega[i] = omeganew[i];
        }

        /*====================Output to files=========================*/
        if (count / countgraph * countgraph == count)
        {
            /*Write cylinder positions*/
            for (igraph = 0; igraph < N; igraph++)
            {
                fprintf(fp, "%f %f %f\n", x[igraph], y[igraph], theta[igraph]);
            }
            fprintf(fp, "\n\n");

            /*Write forces (pair and wall)*/
            for (igraph = 0; igraph < N; igraph++)
            {
                for (jgraph = 0; jgraph < N; jgraph++)
                {
                    if (igraph == jgraph)
                    {
                        continue;
                    }
                    fprintf(fp1, "%.20e %.20e %.20e %.20e\n", x[igraph] - nx[igraph][jgraph], y[igraph] - ny[igraph][jgraph], (fnx[igraph][jgraph] + ftx[igraph][jgraph]), (fny[igraph][jgraph] + fty[igraph][jgraph]));
                }
                for (kgraph = 0; kgraph < wall; kgraph++)
                {
                    if (kgraph == 0)
                    {
                        fprintf(fp1, "%.20e %.20e %.20e %.20e\n", x[igraph] - nx_wall[kgraph], y[igraph] - ny_wall[kgraph], (fnx_wall[igraph][kgraph] + ftx_wall[igraph][kgraph]), (fny_wall[igraph][kgraph] + fty_wall[igraph][kgraph]));
                    }
                    else
                    {
                        fprintf(fp1, "%.20e %.20e %.20e %.20e\n", x[igraph] - nx_wall[kgraph], y[igraph] - ny_wall[kgraph], fnx_wall[igraph][kgraph], fny_wall[igraph][kgraph]);
                    }
                }
            }
            fprintf(fp1, "\n\n");

            /*Write wall geometry*/
            fprintf(fp2, "%d %f %d %d\n", -10, d[0], 20, 0);
            fprintf(fp2, "%f %d %d %d\n", -d[1], 0, 0, 10);
            fprintf(fp2, "%f %d %d %d\n\n\n", d[2], 0, 0, 10);
        }

        /*===========Advance time==========*/
        t += dt;
        count++;
    }
    /*Dump initial conditions*/
    for (i = 0; i < N; i++)
    {
        fprintf(fp3, "%.10e %.10e %.10e %.10e %.10e %.10e\n", x[i], y[i], theta[i], u[i], v[i], omega[i]);
    }
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            fprintf(fp4, "%.10e ", delta_tx[i][j]);
        }
        fprintf(fp4, "\n");
    }
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            fprintf(fp4, "%.10e ", delta_ty[i][j]);
        }
        fprintf(fp4, "\n");
    }
    for (i = 0; i < N; i++)
    {
        fprintf(fp4, "%.10e ", delta_tx_wall[i][0]);
    }
    fprintf(fp4, "\n");
    for (i = 0; i < N; i++)
    {
        fprintf(fp4, "%.10e ", delta_ty_wall[i][0]);
    }

    return 0;
}