/* This program will compute the velocities and positions of a 2-body problem
 *
 * Written by Kyle Kanos, 2015
 */
 
#include <stdio.h>
#include <math.h>

/* Kepler is the position & velocity update. 
   Note that u[0] & u[2] are the positions 
             u[1] & u[3]  are the velocities
*/
void kepler(double t, double *u, double *up){
    double r3 = sqrt(pow(u[0]*u[0] + u[2]*u[2], 3.0));
    up[0] = u[1]; up[1] = -u[0]/r3;
    up[2] = u[3]; up[3] = -u[2]/r3;
}

/* Runge Kutta integrator
   This function does the core work in solveing the problem. It returns
   an int in case you want to do some error checking here
*/
int rkStep(int neqn, double *y, double *yp, double t, double dt){
    double k1[neqn], k2[neqn], k3[neqn], k4[neqn], yt[neqn];
    int i;
  // first step
    kepler(t, y, k1);
    
    for(i=0;i<neqn;i++){
        yt[i] = y[i] + 0.5*dt*k1[i];
    }
  // second step
    kepler(t, yt, k2);
    
    for(i=0;i<neqn;i++){
        yt[i] = y[i] + 0.5*dt*k2[i];
    }
  // third step
    kepler(t, yt, k3);
    
    for(i=0;i<neqn;i++){
        yt[i] = y[i] + dt*k3[i];
    }
  // fourth step
    kepler(t, yt, k4);
    
    for(i=0;i<neqn;i++){
        yp[i] = y[i] + dt*(k1[i] + 2.0*(k2[i] + k3[i]) + k4[i])/6.0;
    }
    return 0;
}

/* Runge Kutta caller
   tend is an arbitrary end point, choose what you want really
*/
void rk4Solve(int neqn, int nsteps, double *y){
    double beg=0.0, tend=7.79;
    
  // initial values of positions and velocities.
  // we chose it starting at x and moving purely in the y dir
     y[0] = 1.0; y[1] = 0.0;  // x, vx
     y[2] = 0.0; y[3] = 1.0;  // y, vy
     
     double dt, t=0.0, yp[4];
   // set time stepper
     dt = (tend - tbeg)/(double) nsteps;
     
     kepler(t, y, yp); // initializes yp is y'
     
     int istep, i;
     FILE *fp;
     fp = fopen("two_body.dat","w");
     
     for(istep=1; istep <= nsteps; istep++){
         rkStep(neqn, y, yp, t, dt, relerr, abserr);
       
       // print results to file
         fprintf(fp,"%e    %e\n",yp[0],yp[2]);
       // update y from yp found in rkStep
         for(i=0; i<neqn; i++) y[i] = yp[i];
     }
     fclose(fp);
 }
     
// obvious ;)
int main(){
    int neqn=4,numstep=1000;
    double x[neqn];
    
    printf("beginning solver...\n");
    rk4Solve(neqn, numstep, x);
    printf("ending program.\n");
    
    return 0;
    
}
    
