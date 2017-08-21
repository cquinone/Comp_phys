#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

/*****************************************
 * Notes-
 *
 *****************************************/

#define RDIM 100
#define INNERDIM 10
#define PHIDIM 64

#define SOLVE_CENTER 0

/*
 * Scaling constants
 *
 * You'll have to pick a value for dt which produces stable evolution
 * for your stencil!
 */
double tMax= 50.0;
double dr= 0.1;
double dphi= 2.0*M_PI/PHIDIM;
double dt= 0.01; /* change this! */

/* 
 * Storage for the arrays.  Ut is U-transformed, the Fourier space version!
 */
double U[RDIM][PHIDIM];
fftw_complex Ut[RDIM][(PHIDIM/2)+1];
fftw_complex UtNew[RDIM][(PHIDIM/2)+1];
fftw_complex UtOld[RDIM][(PHIDIM/2)+1];

/*
 * FFTW 'plans'.  If you need additional plans you can put the globals
 * for them here.
 */
//fftw_plan uForwardPlan= NULL;
//fftw_plan uBackwardPlan= NULL;

double initval( double r, double phi ) {
  /*
   * This routine gives the initial values as a function of r and phi.
   * You'll need to set up additional initial data which is constant 
   * in phi and in r as you are debugging your program!
   */
  double x= r*cos(phi);
  double y= r*sin(phi);
  double ctrX= -2.0;
  double ctrY= 0.0;
  double sigma= 0.25;
  double maxU= 5.0;
  double distSqr= (x-ctrX)*(x-ctrX) + (y-ctrY)*(y-ctrY);
  double U= maxU*exp(-distSqr/(sigma*sigma));

  /*
  if (U >= 0.00001)
    fprintf(stderr,"%f %f -> %f %f -> %f\n",r,phi,x,y,U);
  */
  return U;
}

void copyComplex( fftw_complex* out, fftw_complex* in, long n ) {
  /* Just for clarity */
  long i;
  for (i=0; i<n; i++) out[i]= in[i];
}

void initializeU() {
  /* Apply the initial condition everywhere on the grid */
  long i,j;
  for (i=0; i<RDIM; i++) {
    double r= i*dr;
    for (j=0; j<PHIDIM; j++) {
      double phi= j*dphi;
      // YOOOOO set r, phi as const for now
      //r=1.0;
      //phi=1.0;
      //phi =0;
      U[i][j]= initval(r,phi);
    }
  }
}

void calculateUtFromU() {
  /* Use the forward real-to-complex FFT to calculate U-transformed from U */
  fftw_plan p;
  double Ur[PHIDIM];
  fftw_complex Utr[(PHIDIM/2)+1];
  //p = fftw_plan_dft_r2c_2d(RDIM, PHIDIM, U, Ut,FFTW_ESTIMATE );
  //fftw_plan p = makeplan();
  //fftw_execute(p);
  //fftw_destroy_plan(p);
  int i,j;
  for (i=10; i<RDIM; i++) {
    //grab  row
    for (j=0; j<PHIDIM; j++){
      Ur[j] = U[i][j];
    }

    p = fftw_plan_dft_r2c_1d(PHIDIM, Ur, Utr,FFTW_ESTIMATE );
    fftw_execute(p);
    fftw_destroy_plan(p);
    // add utr to ut
    for (j=0; j<(PHIDIM/2)+1; j++) {
      Ut[i][j] = Utr[j];
    }
  }

}

void calculateUFromUt() {
  /* Use the inverse transform to get back from Fourier space to real space */

  fftw_plan g;
  double Ur[PHIDIM];
  fftw_complex Utr[(PHIDIM)/2 +1];

  int i,j;
  for (i=0; i<RDIM; i++) {
    for (j=0; j<(PHIDIM)/2 +1; j++){
      Utr[j] = Ut[i][j];
    }
    //printf("%s\n", "size of Utr: ");
    //printf("%zu\n", sizeof(Utr)/sizeof(fftw_complex));
    
    g = fftw_plan_dft_c2r_1d(PHIDIM, Utr, Ur,FFTW_ESTIMATE );
    fftw_execute(g);
    fftw_destroy_plan(g);
    for (j=0; j<PHIDIM; j++) {
      U[i][j] = Ur[j];
    }
 
  }
  // shift by 1/n to correct values
  for (i=0; i<PHIDIM; i++) {
    for (j=0; j<RDIM; j++) {
      if (j < 10 || j == 10) {
        U[j][i] = .025;
      }
      U[j][i]= U[j][i] / (double)(PHIDIM); 
    }
  }

}


void doTimeStep() {
long i,j;
for (j=10; j<RDIM-1; j++) {
  double r= j*dr; 
  for (i=0; i<(PHIDIM/2)+1; i++) {
    // I have 1*1 just co i remember that c^2 is included
    UtNew[j][i]= (dt*dt)*(1*1)*((Ut[j+1][i]+Ut[j-1][i]-2*Ut[j][i])/(dr*dr) + (1/(2*dr*r))*(Ut[j+1][i] - Ut[j-1][i]) - ((i*i)/(r*r))*Ut[j][i]) - UtOld[j][i] + 2*Ut[j][i];
    //if (j == RDIM -2 || j > RDIM - 2){
      //printf("%s", "hey");
    //  UtNew[j][i] = 1;
    //}
  }
}

}

int timeToOutput(double t, int count) {
  /* A little test to tell how often to dump output */
  return (count%5 == 0);
}

void outputU(double t, int count) {
  /*
   * You may need to adapt this to the needs of the tool you use
   * to make your animation.
   */
  char buf[256];
  char buf2[256];
  double Udiag[PHIDIM][RDIM];
  FILE* ofile= NULL;
  FILE* ofile2= NULL;
  double max;
  double min;
  int i,j;
  static int n= 0;

  sprintf(buf,"out_%03d.raw",n);

  /*
   * FFTW expects data in row-major order, but VisIt wants it in column-
   * major order.  We just make a flipped copy for output.
   */
  for (i=0; i<RDIM; i++)
    for (j=0; j<PHIDIM; j++) 
      Udiag[j][i]= U[i][j];

  ofile= fopen(buf,"w");
  (void)fwrite(Udiag,sizeof(double),RDIM*PHIDIM,ofile);
  (void)fclose(ofile);

  snprintf(buf2,sizeof(buf2),"out_%03d.bov",n);
  ofile2= fopen(buf2,"w");
  if (!ofile2) {
    perror(buf2);
    exit(-1);
  }
  //fprintf(stderr,"Opened %s\n",buf2);
  fprintf(ofile2,"TIME: %g\n",t);
  fprintf(ofile2,"DATA_FILE: %s\n",buf);
  fprintf(ofile2,"DATA_SIZE: %d %d 1\n",RDIM,PHIDIM);
  fprintf(ofile2,"DATA_FORMAT: DOUBLE\n");
  fprintf(ofile2,"VARIABLE: U\n");
  fprintf(ofile2,"DATA_ENDIAN: LITTLE\n");
  fprintf(ofile2,"CENTERING: ZONAL\n");
  fprintf(ofile2,"BRICK_ORIGIN: 0. 0. 0.\n");
  fprintf(ofile2,"BRICK_SIZE: 1.0 1.0 1.0\n");
  fclose(ofile2);
  n++;

  for (i=0; i<RDIM; i++)
    for (j=0; j<PHIDIM; j++) {
      if (i==0 && j==0) {
         min= max= U[i][j];
      }
      else {
  if (U[i][j]<min) min= U[i][j];
  if (U[i][j]>max) max= U[i][j];
      }
    }

  //printf("Output %s; min= %f, max= %f\n",buf,min,max);
}

int main( int argc, char* argv[] ) {
  double t= 0.0;
  int count= 0;

  //makePlans();

  initializeU();
  calculateUtFromU();

  /* 
   * We start out with two copies of the initial condition, 
   * so you'll have something to leap from if you use a leapfrog
   * stencil.
   */
  copyComplex(&(UtOld[0][0]),&(Ut[0][0]),RDIM*((PHIDIM/2)+1));

  for (t=0.0; t<tMax; t+=dt) {
    if (timeToOutput(t,count)) {
      calculateUFromUt();
      outputU(t,count);
    }
    doTimeStep();
    copyComplex(&(UtOld[0][0]),&(Ut[0][0]),RDIM*((PHIDIM/2)+1));
    copyComplex(&(Ut[0][0]),&(UtNew[0][0]),RDIM*((PHIDIM/2)+1));
    count += 1;
  }
}