/*
    Animation of stacking three cylinders (2D).
    This program uses gnuplot to generate an animated GIF from simulation snapshots.
*/
#include <stdio.h>
#include <assert.h>

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

int main(void)
{
    int itend = tend / dt;                    // Total steps
    int countgraph = (int)(tend / dt / data); // write to files every "countgraph" steps

    /*=============Launch gnuplot==============*/
    FILE *gp;

    gp = popen("gnuplot", "w");
    assert(gp);

    fprintf(gp, "set terminal gif animate optimize delay 2 size 1000, 1000 \n"); // delay n size <size> : pose n*0.001 sec
    fprintf(gp, "set output 'fig/A%.1e.gif'\n", A);
    fprintf(gp, "set xrange[-3:3]\n");
    fprintf(gp, "set yrange[0:6]\n");
    fprintf(gp, "set size square\n");

    for (int j = 0; j < data; j++)
    {
        fprintf(gp, "set title 't = %f'\n", j * dt * countgraph); // title

        fprintf(gp, "plot 'data/p_A%.1e.txt' index %d using 1:2:(1) with circles lc -1 lw 3 notitle,'data/p_A%.1e.txt' index %d using 1:2:(cos($3)):(sin($3)) with vectors nohead lc -1 lw 3 notitle, 'data/f_A%.1e.txt' index %d using 1:2:(1*$3):(1*$4) with vectors lc 'red' lw 3 notitle, 'data/w_A%.1e.txt' index %d using 1:2:3:4 with vectors lc -1 lw 3 notitle\n", A, j, A, j, A, j, A, j);
    }
    pclose(gp);

    return 0;
}