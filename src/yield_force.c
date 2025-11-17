/*
Overview:
    Simulation of three vertically stacked cylinders compressed slowly from above by a wall.
    The normal force acting on the wall is computed over time, and its maximum value
    is recorded as the yield force.

Conditions:
    Perform simulations with spring constants kn = 1e+2, 1e+3, and 1e+4,
    corresponding to ratio = 1, 10, and 100, respectively.

Output:
    Wall friction coefficient (mu_wall), yield force (fmax), and
    the time (t) at which the maximum occurs.
*/

#include <stdio.h>
#include <math.h>

/*=========User-defined parameters (to be manually changed for each simulation)==========*/
#define ratio 10.0                // Stiffness scale; kn and kt are multiplied by this factor.
#define mu_wall 0.2               // Floor-cylinder maximum static friction coefficient.
#define PI 3.14159265358979323846 // π

/*=========Material/Model parameters (dimensionless)==========*/
#define N 3               // Number of cylinders.
#define wall 2            // Number of walls (0 = fixed floor, 1 = moving top wall).
#define r 1.0             // Cylinder radius.
#define e 0.3             // Coefficient of restitution
#define mu 0.7            // Cylinder-cylinder maximum static friction coefficient.
#define kn 1e+2 * ratio   // Normal spring constant.
#define kt 4e+1 * ratio   // Tangential spring constant.
#define g 1.0             // Gravitational acceleration.
#define vwall 1e-6        // Compressing velocity of the top wall.
#define dt 1e-3           // Time step (for ratio = 1 or 10, dt = 1e-2 is also acceptable).
const double TMAX = 1e10; // Limit of time

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
    FILE *file;
    int i, j, k;
    double t; // Simulation time.

    /*==============Particle state==============*/
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

    /*=============物性パラメータ=================*/
    double m = 1.0;                                                         // Cylinder mass.
    double I = 0.5 * m * r * r;                                             // Moment of inertia of a cylinder.
    double etan = -2 * log(e) * sqrt(m * kn / (PI * PI + log(e) * log(e))); // Normal damping.
    double etat = -2 * log(e) * sqrt(m * kt / (PI * PI + log(e) * log(e))); // Tangential damping.

    /* Determine indentation depth A based on ratio (1 → 1e-3, 10 → 1.2e-4, 100 → 2e-5) */
    double A;
    if (ratio == 1.0)
    {
        A = 1.0e-3;
    }
    else if (ratio == 10.0)
    {
        A = 1.2e-4;
    }
    else if (ratio == 100.0)
    {
        A = 2.0e-5;
    }
    else
    {
        printf("Error: invalid ratio value (%.2f). Please choose ratio = 1, 10, or 100.\n", ratio);
        return 1; // または exit(EXIT_FAILURE);
    }

    /*Yield-force tracking*/
    double fmax = -1; // Running maximum of monitored normal force
    double cur = 0.0; // Current monitored value

    /*================Initialization===================*/
    t = 0.0;

    /*--------Load initial positions & velocities----------*/
    char filename[100];
    sprintf(filename, "init/init_p_A%.1e.txt", A);
    file = fopen(filename, "r");
    if (!file)
    {
        fprintf(stderr, "Failed to open %s\n", filename);
        return 1;
    }
    for (i = 0; i < N; i++)
    {
        fscanf(file, "%lf %lf %lf %lf %lf %lf", &x[i], &y[i], &theta[i], &u[i], &v[i], &omega[i]);
    }
    fclose(file);

    /*--------Load initial tangential displacements----------*/
    sprintf(filename, "init/init_deltat_A%.1e.txt", A);
    file = fopen(filename, "r");
    if (!file)
    {
        fprintf(stderr, "Failed to open %s\n", filename);
        return 1;
    }

    /*cylinder-cylinder*/
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            fscanf(file, "%lf", &delta_tx[i][j]);
        }
    }

    /*cylinder-wall (only for wall index 0)*/
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            fscanf(file, "%lf", &delta_ty[i][j]);
        }
    }
    for (i = 0; i < N; i++)
    {

        fscanf(file, "%lf", &delta_tx_wall[i][0]);
    }
    for (i = 0; i < N; i++)
    {

        fscanf(file, "%lf", &delta_ty_wall[i][0]);
    }
    fclose(file);

    /*--------initial updated tangential displacement---------*/
    matrix_init0(N, N, delta_txnew);
    matrix_init0(N, N, delta_tynew);

    /*--------Wall setup: floor (k=0), moving top wall (k=1)--------*/
    init(0, nx_wall, ny_wall, d, 0.0, 1.0, 0.0); // Floor: offset y=0, normal vector (0,1)
    u_wall[0] = 0.0;
    v_wall[0] = 0.0;

    init(1, nx_wall, ny_wall, d, 0.0, -1.0, 2 + sqrt(3)); // Top: offset y=2+sqrt(3), normal vector (0,-1)
    u_wall[1] = 0.0;
    v_wall[1] = 0.0;

    matrix_init0(N, wall, delta_txnew_wall);
    matrix_init0(N, wall, delta_tynew_wall);

    /*==================Time evolution==================*/
    while (1)
    {
        /*---Top-wall motion (hold for t<1, then move downward at vwall)---*/
        if (t < 1)
        {
            d[1] = 2 + sqrt(3);
            v_wall[1] = 0.0;
        }
        else
        {
            d[1] = 2 + sqrt(3) - vwall * (t - 1);
            v_wall[1] = -vwall;
        }

        /*============Contact calculations=============*/
        for (i = 0; i < N; i++)
        {
            /*---------cylinder-cylinder interactions----------*/
            for (j = 0; j < N; j++) // interaction with cylinder j
            {
                if (j == i)
                {
                    continue;
                }
                // printf("%d %d\n", i, j);

                /*Relative position/velocity and distance*/
                lx[i][j] = x[i] - x[j];
                ly[i][j] = y[i] - y[j];
                du[i][j] = u[i] - u[j];
                dv[i][j] = v[i] - v[j];
                l[i][j] = sqrt(pow(lx[i][j], 2) + pow(ly[i][j], 2));

                /*Normal unit vector (i ← j)*/
                nx[i][j] = lx[i][j] / l[i][j];
                ny[i][j] = ly[i][j] / l[i][j];

                /*Normal overlap vector (negative when overlapping)*/
                delta_nx[i][j] = (l[i][j] - 2 * r) * nx[i][j];
                delta_ny[i][j] = (l[i][j] - 2 * r) * ny[i][j];

                /*Normal relative velocity*/
                v_nx[i][j] = (du[i][j] * nx[i][j] + dv[i][j] * ny[i][j]) * nx[i][j];
                v_ny[i][j] = (du[i][j] * nx[i][j] + dv[i][j] * ny[i][j]) * ny[i][j];

                /*Tangential relative velocity (includes rotation)*/
                vtx[i][j] = du[i][j] - (du[i][j] * nx[i][j] + dv[i][j] * ny[i][j]) * nx[i][j] + r * (omega[i] + omega[j]) * ny[i][j];
                vty[i][j] = dv[i][j] - (du[i][j] * nx[i][j] + dv[i][j] * ny[i][j]) * ny[i][j] - r * (omega[i] + omega[j]) * nx[i][j];

                /*Magnitudes used for unit tangents*/
                deltat[i][j] = sqrt(pow(delta_tx[i][j], 2) + pow(delta_ty[i][j], 2));
                vt[i][j] = sqrt(pow(vtx[i][j], 2) + pow(vty[i][j], 2));

                /*Tangential unit vector*/
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

                /*Normal contact force*/
                if (l[i][j] - 2 * r < 0) // contact
                {
                    fnx[i][j] = -kn * delta_nx[i][j] - etan * v_nx[i][j];
                    fny[i][j] = -kn * delta_ny[i][j] - etan * v_ny[i][j];
                }
                else // no contact
                {
                    fnx[i][j] = 0;
                    fny[i][j] = 0;
                }

                /*Trial tangential force*/
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

                /*Update of tangential displacement*/
                if (ft[i][j] >= mu * fn[i][j] && l[i][j] - 2 * r < 0) // slip
                {
                    ftx[i][j] = -mu * fn[i][j] * tx[i][j];
                    fty[i][j] = -mu * fn[i][j] * ty[i][j];
                    delta_txnew[i][j] = mu * kn * (2 * r - l[i][j]) / kt * tx[i][j];
                    delta_tynew[i][j] = mu * kn * (2 * r - l[i][j]) / kt * ty[i][j];
                }
                else if (ft[i][j] < mu * fn[i][j] && l[i][j] - 2 * r < 0) // no-slip
                {
                    delta_txnew[i][j] = delta_tx[i][j] + vtx[i][j] * dt;
                    delta_tynew[i][j] = delta_ty[i][j] + vty[i][j] * dt;
                }
                else // no contact
                {
                    delta_tx[i][j] = 0;
                    delta_ty[i][j] = 0;
                }

                /*Commit tangential displacements*/
                delta_tx[i][j] = delta_txnew[i][j];
                delta_ty[i][j] = delta_tynew[i][j];
            }

            /*-------------Cylinder-wall interactions--------------*/
            for (k = 0; k < wall; k++)
            {
                /*Signed distance to wall plane*/
                l_wall[i][k] = nx_wall[k] * x[i] + ny_wall[k] * y[i] + d[k];

                /*Normal overlap (negative in compression)*/
                delta_nx_wall[i][k] = (l_wall[i][k] - r) * nx_wall[k];
                delta_ny_wall[i][k] = (l_wall[i][k] - r) * ny_wall[k];

                /*Normal relative velocity (cylinder vs wall)*/
                vnx_wall[i][k] = ((u[i] - u_wall[k]) * nx_wall[k] + (v[i] - v_wall[k]) * ny_wall[k]) * nx_wall[k];
                vny_wall[i][k] = ((u[i] - u_wall[k]) * nx_wall[k] + (v[i] - v_wall[k]) * ny_wall[k]) * ny_wall[k];

                /*Tangential relative velocity (includes rotation)*/
                vtx_wall[i][k] = (u[i] - u_wall[k]) - vnx_wall[i][k] + r * omega[i] * ny_wall[k];
                vty_wall[i][k] = (v[i] - v_wall[k]) - vny_wall[i][k] - r * omega[i] * nx_wall[k];

                /*Magnitudes for unit tangents*/
                deltat_wall[i][k] = sqrt(pow(delta_tx_wall[i][k], 2) + pow(delta_ty_wall[i][k], 2));
                vt_wall[i][k] = sqrt(pow(vtx_wall[i][k], 2) + pow(vty_wall[i][k], 2));

                /*Tangential unit vector at wall*/
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

                /*Normal contact force with wall*/
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

                /*Trial tangential force with wall*/
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

                /*Update tangential displacement*/
                if (ft_wall[i][k] >= mu_wall * fn_wall[i][k] && l_wall[i][k] - r < 0) // slip
                {
                    ftx_wall[i][k] = -mu_wall * fn_wall[i][k] * tx_wall[i][k];
                    fty_wall[i][k] = -mu_wall * fn_wall[i][k] * ty_wall[i][k];
                    delta_txnew_wall[i][k] = mu_wall * kn * (r - l_wall[i][k]) / kt * tx_wall[i][k];
                    delta_tynew_wall[i][k] = mu_wall * kn * (r - l_wall[i][k]) / kt * ty_wall[i][k];
                }
                else if (ft_wall[i][k] < mu_wall * fn_wall[i][k] && l_wall[i][k] - r < 0) // no-slip
                {
                    delta_txnew_wall[i][k] = delta_tx_wall[i][k] + vtx_wall[i][k] * dt;
                    delta_tynew_wall[i][k] = delta_ty_wall[i][k] + vty_wall[i][k] * dt;
                }
                else
                {
                    delta_txnew_wall[i][k] = 0;
                    delta_tynew_wall[i][k] = 0;
                }

                /*Commit wall tangential displacement*/
                delta_tx_wall[i][k] = delta_txnew_wall[i][k];
                delta_ty_wall[i][k] = delta_tynew_wall[i][k];
            }
        }

        /*===========Accumulate forces and integrate==============*/
        for (i = 0; i < N; i++)
        {
            fx[i] = 0.0;
            fy[i] = 0.0;
            T[i] = 0.0;
        }

        for (i = 0; i < N; i++)
        {

            /*Cylinder-cylinder resultant*/
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

            /*Cylinder-wall resultant*/
            for (k = 0; k < wall; k++)
            {
                fx[i] += fnx_wall[i][k] + ftx_wall[i][k];
                fy[i] += fny_wall[i][k] + fty_wall[i][k];
                T[i] += -r * (nx_wall[k] * fty_wall[i][k] - ny_wall[k] * ftx_wall[i][k]);
            }

            /*Gravity*/
            fy[i] += -g;

            /*Leapfrog*/
            unew[i] = u[i] + fx[i] * dt;
            vnew[i] = v[i] + fy[i] * dt;
            omeganew[i] = omega[i] + T[i] * dt / I;

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

        /*Monitor the normal force at top wall on cylinder 2: stop after peak*/
        cur = -fny_wall[2][1]; // compressive normal component from top wall

        if (fmax < 0.0)
        {
            fmax = cur;
        }
        else
        {
            if (cur < fmax)
            {
                printf("mu_wall: %f, fmax: %f, tend: %f\n", mu_wall, fmax, t);
                break;
            }
            else
            {
                fmax = cur;
            }
        }

        /*Advance time*/
        t += dt;

        if (t > TMAX)
        {
            fprintf(stderr, "Warning: peak not found before limit. mu_wall=%g, fmax=%g, t=%g\n", mu_wall, fmax, t);
            break;
        }
    }

    return 0;
}