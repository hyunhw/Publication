///////////////////////////////////////////////////////////////////////
//          < APPENDED TO SHEAR.CPP >                                //
// Singletrajotf.cpp-----------------> energy.dat / post-processing  //
//  |<- coarseningtraj.xyz          Coarsen step max = 0;            //
//  |<- rnl.dat                                                      //
//                                                                   //
// This file reads in the desired time frame(n) of my simulation,    //
//  Time frames are based off of MATLAB post processing values (i)'s //
//  This means the same index in C will be (i-1)'s                   //
//                                                                   //
//   -get the xyz coordinates at (i-2)-initial, (n-1)-final          //
//   -ripens the particles via original rules,                       //
//   -and relaxes particles via steepest descent                     //
//                                                                   //
//                                                                   //
//    SD -> delt = 0.1                                               //
//    FIRE -> delt = 0.002                                           //
//                                                                   //
//    Manually switch between calling steepest_descent() or fire()   //
//                                                                   //
///////////////////////////////////////////////////////////////////////


#include <cstdio>
#include <time.h>
#include <cstdlib>
#include <cmath>
#define PI 3.141592653589793238462643383

//  xx -> initial reading from coarseningtraj.xyz
//  |---> store in *init for final dr comparison (after minimization)
//  x -> variable I work with to minimize (from original code)
//   |---> FIRE
//   |---> Steepest descent
//  xxx -> variable I work with to post-process the trajectory


/////////////////////////////
// Define global variables //
int frame, nf, nf2, fr_max, fr_max2, input, Zstart, Zend;
double **oldx, ***vv, ***xx, ***xxx, **dd, **init, *cl, *Rend, *cl_frame, *Usys, *Udiff;
double timesum;
int origin; //variable used in MSD2 routine



////////////////////////////
// Define global function //
void read_input( void );
void read_traj( void );
void read_rnl( void );
void store_coordinate( void );
void remove_PBC( void );
void steepest_descent( void );
void MSD( void );
void MSD2( void );
double calc_MSD2( void );
void write_energy( void );
void get_contour( void );




/*------------------------*/
/* Define Global Variable */
/*------------------------*/
//----> original SHEAR.CPP variables
int N, N_tracer;                                                                        // Main MD
int print_freq, neighbor_freq, neighbor_freq2;
double T, PE, KE, vol, delt, box[3], boxh[3];
double **x, **v, **f, *d, **s, **elastic, **viscous;
double rcut, rcut2, b;
/*------------------------------------------------------------------------------------------------*/
long int idum; //Random number seed
int *c, avgc, mn, mx; 
int do_hessian; 
/*------------------------------------------------------------------------------------------------*/
int **n, *n_neighbor, *n_cell, ***n_in_cell, ****cell, **a, n_max;         // Neighbor * Cell list 
double *l_cell, rnl, rnl2, max_cl;                                           
/*------------------------------------------------------------------------------------------------*/
int start, end, id;                                                                       // Shear 
int *l; 
double **w; 
/*------------------------------------------------------------------------------------------------*/
int Nmin, step, firestep_max, fire_flag, rcount;                                           // FIRE 
double finc, fdec, deltmax, falpha, alphastart, *rmsd;    
/*------------------------------------------------------------------------------------------------*/
int coarsen_step;
int coarsen_print_freq;
int n_flag, *flag;
double ZZ, ZC, Ziso;
double C, Cf, r_res, *Q;                                                             // Coarsening
/*------------------------------------------------------------------------------------------------*/
double **fin, **pre_fin, meandistance, *msd, Pr;                                            // MSD 
/*------------------------------------------------------------------------------------------------*/
double **H, **ff, **fb, **minx, **minf, delx;                                           // Hessian
/*------------------------------------------------------------------------------------------------*/
double ***stracer;                                                                       // Stress

double inivol, inivol2, cdum;
int mass_freq;



/*-------------------------*/
/* Define Global Functions */
/*-------------------------*/

double U_system (void);
double pbc_vdr (double*, double*, double*);
double MSD (double*, double*, double*);
double calc_MSD (void);

void neighbor (void);
void neighbor_N2 (void);
void cell_memory (void);
void assign_cell (void);
void rnl_update (void);

void fire (void);
void coarsening (void);
void mass_consv (void);
void hessian_ana (void);
void hessian_num (void);
void hessian_memory (void);
void forward_force (int, int, double);
void backward_force (int, int, double);

void print_hessian (void);
void print_coordi (void);
void print_msd(void);
void print_Pr(void);
void print_diameter (void);
void print_pov_info (void);
void print_initial_diameter (void);
void print_net_diameter (void);
void print_max_diam (void);
void print_coordinates (void);
void print_tracer_coordination (void);
void print_tracer_stress (void);
void print_iso(void);

void write_xyz (void);
void write_coarsen_xyz(void);
void write_tracer_xyz(void);
double ran2 (void);
double dgauss (double, double);

/*------------------------------------------------*/
/*                 Start of Main                  */ 
/*------------------------------------------------*/

int main ( int argc, char** argv) {

  FILE *coarsen, *Zinput, *inp;
  double sum3;
  int i, j, k;
  time_t timer; //type long int

  //read in Z window index so you can pick random coarsening step inside the window
  Zinput = fopen("zwindow.dat","r");
  fscanf(Zinput, "%d %d\n", &Zstart, &Zend);
  fclose(Zinput);
  if (Zstart == 0 && Zend == 0){
    printf("window is not defined, terminating simulation\n");
    exit(1);
  }
    
  time(&timer);
  idum = -timer;
  printf("seed: %ld\n", idum);

//  input = int((Zend-Zstart)*ran2())+Zstart;
//  printf("input coarsening step: %d\n", input);

/*    if (argc < 3){
        printf("Desired input frame / origin is missing!\n");
        exit(1);
    }
    input = atoi (argv[1]);
    origin = atoi (argv[2]);
*/


  // Define simulation parameters //
  /*------------------------------*/
  //idum = -long(time(0));
    fr_max = 5000; //some number to allocate memory for initial trajectory (ORIGINAL SHEAR RUN)
    fr_max2 = 0; //frame number of FIRE or STEEPEST descnt run
    // |---> this number could be very large depending on how long the FIRE or STEEPEST DESCENT runs for
    // |---> this number is determined during FIRE or STEEPEST DESCENT algorithm. Thus xxx memory is allocated afterwards
    timesum = 0.0;
    
  //idum = -3140;
  do_hessian = 0;

  firestep_max = 1000000000;
  end = 100000;
  start = 1000000000;

  //N = 2000;
  n_flag = 0;
  //vol = 2503.5544426715;
  b = 0.0;
  C = 0.05;
  Cf = 0.002;
  ZC = 6.0;

  rcut = 2.5;
  rcut2 = rcut * rcut ;
  rnl = 3.0;
  rnl2 = rnl * rnl ;

  delt = 0.001;
  deltmax = 0.1;

  print_freq = 1;
  neighbor_freq = 1;
  neighbor_freq2 = 10;
  coarsen_print_freq = 1;
  mass_freq = 49;

  alphastart = 0.1;
  finc = 1.1;
  fdec = 0.5;
  falpha = 0.99;
  Nmin = 5;
  /*------------------------------*/
  meandistance = 0.0;
  coarsen_step = -1;

  inp = fopen("coarseningtraj.xyz","r");
  fscanf(inp, "%d %lf\n", &N, &vol);
  fclose(inp);
  printf("N: %d vol: %lf\n", N, vol);
  N_tracer = int(N*0.005);

  // Set up simulation box
  for (i=0; i<3; i++){
    box[i] = pow ( vol, 1.0/3.0 );
    boxh[i] = box[i] * 0.5;
  }

  //define max cell length (min 3 cells)
  max_cl = box[0]/3.0;

  // Initialize center of mass velocity to 0
  double cm[3];
  cm[0] = cm[1] = cm[2] = 0.0;


  //cl_frame = (double*) calloc (nf, sizeof (double));
  // Allocate memory for various parameters - singletraj
  cl = (double*) calloc (N, sizeof (double));
  Rend = (double*) calloc (N, sizeof (double));
  xx = (double***) calloc(fr_max, sizeof(double**));
  vv = (double***) calloc(fr_max, sizeof(double**));
  dd = (double**) calloc(fr_max, sizeof(double*));
    for (i=0; i<fr_max; i++){
        xx[i] = (double**) calloc ( N , sizeof ( double* ));
        vv[i] = (double**) calloc ( N , sizeof ( double* ));
        dd[i] = (double*) calloc ( N , sizeof ( double ));
        for (j=0; j<N; j++){
            xx[i][j] = (double*) calloc ( 3 , sizeof ( double ));
            vv[i][j] = (double*) calloc ( 3 , sizeof ( double ));
        }
    }
    
     init = (double**) calloc( N , sizeof(double*));
    for (i=0; i<N; i++){
        init[i] = (double*) calloc ( 3, sizeof ( double ));
    }

  // Allocate memory for various parameters
  //
  //  **s - total stress
  //  **elastic - elastic stress
  //  **viscous - viscous stress
  //  ***stracer - stress on tracer particles
  //  **H - Hessian matrix
  //
  //  **x, **v and **f, *d (diameter for particle type), **w (initial position), *l (label)
  //  *c - coordination number
  //  *area, *delP - overlapping area and pressure difference for coarsening
  //
  s = (double**) calloc (3, sizeof (double*));
  elastic = (double**) calloc (3, sizeof (double*));
  viscous = (double**) calloc (3, sizeof (double*));
  stracer = (double***) calloc (N_tracer, sizeof (double**));

  for (i=0; i<3; i++){
    s[i] = (double*) calloc (3, sizeof (double));
    elastic[i] = (double*) calloc (3, sizeof (double));
    viscous[i] = (double*) calloc (3, sizeof (double));
  }

  for (i=0; i<N_tracer; i++){
    stracer[i] = (double**) calloc (3, sizeof (double*));
    for (j=0; j<3; j++){
      stracer[i][j] = (double*) calloc (3, sizeof (double));
    }
  }
  
  //Allocate memory for coarsening
  Q = (double*) calloc (N, sizeof (double));
  rmsd = (double*) calloc (N, sizeof (double));
  flag = (int*) calloc (N, sizeof (int));

  x = (double**) calloc (N, sizeof (double*));
  w = (double**) calloc (N, sizeof (double*));
  v = (double**) calloc (N, sizeof (double*));
  f = (double**) calloc (N, sizeof (double*));
  fin = (double**) calloc (N, sizeof (double*));
  pre_fin = (double**) calloc (N, sizeof (double*));
  msd = (double*) calloc (N, sizeof (double*));
  d = (double*) calloc (N, sizeof (double));
  l = (int*) calloc (N, sizeof (int));
  c = (int*) calloc (N, sizeof (int));
  
  // Initialize velocity and particles
  for (i=0; i<N; i++){
    x[i] = (double*) calloc (3, sizeof (double));
    w[i] = (double*) calloc (3, sizeof (double));
    v[i] = (double*) calloc (3, sizeof (double));
    f[i] = (double*) calloc (3, sizeof (double));
    fin[i] = (double*) calloc (3, sizeof (double));
    pre_fin[i] = (double*) calloc (3, sizeof (double));
    rmsd[i] = 0.0;
    flag[i] = 1;

    for (j=0; j<3; j++){
      // Randomly place particles in the box with the origin set to the middle of the box
      x[i][j] = (ran2() * box[j]) - boxh[j];
      // Assign random velocities to the particles
      v[i][j] = dgauss (0.0, 1.0);
      // Accumulate center of mass velocity
      cm[j] += v[i][j];
    }
  }


  ///////////////////
  //  SINGLE TRAJ  //
  ///////////////////

  cell_memory();

  int maxf, d1, d2, d3, d4;
  double d5, maxpr;
  //Read in ncstats.dat
  inp = fopen("ncstats.dat","r");

  maxpr = -1.0;
  while(!feof(inp)){
    fscanf(inp, "%d %d %d %d %lf\n", &d1, &d2, &d3, &d4, &d5); 
    if ( double(d3)/double(d4) > maxpr ){
      maxpr = double(d3)/double(d4);
      maxf = d1;
    }
  }//end of while loop
  fclose(inp);

  printf("frame: %d, pr: %lf\n", maxf, maxpr);
  input = maxf-1;


  read_input();
  store_coordinate();
  read_rnl();

  //store initial coordinates for later MSD calculation inside fire routine
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      pre_fin[i][j] = init[i][j];
    }
  }
 
  printf("rnl: %1.16lf\n", rnl);

  int a1,b1,c1;
  a1 = b1 = c1 = 0; //a1=flag 0, b1=flag 1, c1=flag 2;
  //check to see if evreything is "read" correctly
  for (i=0; i<N; i++){
    if (flag[i] == 2) c1+=1;
    if (flag[i] == 0) a1+=1;
    b1 = N - c1 - a1;
  }
  printf("flag[i]:  0: %d  1: %d 2: %d\n", a1, b1, c1);

  printf("check1\n"); fflush(stdout);

  //memory allocated for all neighbor routines in cell_memory
  //assign_cell();
  neighbor_N2();
  //neighbor();

  printf("check1\n"); fflush(stdout);
  PE = U_system();
  printf("Initial PE (before coarsening): %1.16lf\n", PE);

  printf("check1\n"); fflush(stdout);
  coarsening();
  rnl_update();






  //b/c my coarsening iteration is for input + 1 step, (input+1) is what's considered for mass update
  if ((input+1) % mass_freq == 0) {
    mass_consv();
  }

  PE = U_system();
  printf("Initial PE (after coarsening): %1.16lf\n", PE);


  //store initial coordinates for total MSD calculation 
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      pre_fin[i][j] = init[i][j];
    }
  }


  fire();
//  steepest_descent();



/*  do_hessian = 1;
  hessian_memory();
  hessian_ana();
  print_hessian();  */

  //MEMORY DECLARED IN READ_TRAJ FUNCTION
//  cl_frame = (double*) calloc (nf, sizeof (double));
//  Usys = (double*) calloc (nf, sizeof (double));
//  Udiff = (double*) calloc (nf, sizeof (double));



  return 0;
}

//Calculate dr_ij with PBC
//Don't read input file in pbc-vdr for efficiency because this is in the U_system loop
double pbc_vdr ( double x1[3], double x2[3] , double dr[3] ){
  double mdr2 = 0.0;
  int i;

  //PBC in all x,y,z directions
  for (i=0; i<3; i++){
    dr[i] = x1[i] - x2[i];

    if (dr[i] > boxh[i]){
      dr[i] -= box[i];
    }
    else if ( dr[i] < -boxh[i] ){
      dr[i] += box[i];
    }
    mdr2 += dr[i] * dr[i];
  }
  return mdr2;
}

void hessian_ana (){
  int i, j, k, z, zz;
  double dr[3], dr2, rij, dij, dij2, h1, h2;
  if (do_hessian){

    // Initialize Hessian to 0
    for (i=0; i<3*N; i++){
      for (j=0; j<3*N; j++){
        H[i][j]= 0.0;
      }
    }

    // Loop over all particles (using neighbor list) to calculate the hessian matrix
    for (i=0; i<N; i++){
      if (flag[i] ==0)
        continue;
      for (k=0; k<n_neighbor[i]; k++){

        j = n[i][k];

        if (flag[j] ==0)
          continue;
        if (j <= i)
          continue;

        dr2 = pbc_vdr(x[i], x[j], dr);
        dij = (d[i] + d[j]) / 2.0;
        dij2 = dij * dij;

        if (dr2 >= dij2)
          continue;

        rij = sqrt (dr2);

        // Coefficients calculated based on i != j hessian derivation
        h1 = -(2.0/dij2) + (2.0/dij*rij); //added only on diagonal components (z==zz)
        h2 = -(2.0/dij*rij*dr2);

        // Loop over x,y,z of particle i and j
        for (z=0; z<3; z++){
          for (zz=0; zz<3; zz++){

            h1 = 2.0 / dij2 / dr2 * dr[z] * dr[zz]
               + 2.0 * ( 1.0 - rij / dij ) / ( dij * dr2 * rij )  * dr[z] * dr[zz] ;
            
            if ( z == zz )
              h2 = -2.0 * ( 1.0 - rij / dij ) / ( dij * rij ) ;
            else
              h2 = 0.0 ;

            H[3*i+z][3*i+zz] += h1 + h2 ;
            H[3*j+z][3*j+zz] += h1 + h2 ;
       
            H[3*i+z][3*j+zz] -= h1 + h2 ;
            H[3*j+zz][3*i+z] -= h1 + h2 ;
          }
        }
      }
    }
  }
}

void print_hessian(){
  int i, j;
  FILE *hessian;

  hessian = fopen ("hessian.dat","w");

  fprintf (hessian, "%d\n", N);

  for (i=0; i<3*(N); i++){
    for (j=0; j<3*(N); j++){
      fprintf (hessian, "%1.16lf ", H[i][j]);
    }
    fprintf (hessian, "\n");
  }
  fclose (hessian);
}


//Calculate the total potential energy of system and forces on all particles
double U_system () {
  int i,j,k,z,zz;
  double Utotal, dr[3], dr2, rij, dij, dij2, dv[3], p;
  double F_pre;
  Utotal = 0.0;
  //Initialize force on all particles to 0
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      f[i][j] = 0.0;
    }
  }

  //Initialize stress on tracers to 0
  for (i=0; i<N_tracer; i++){
    for (j=0; j<3; j++){
      for (k=0; k<3; k++){
        stracer[i][j][k] = 0.0;
      }
    }
  }

  //Initialize stress on all particles to 0
  for (i=0; i<3; i++){
    for (j=0; j<3; j++){
      s[i][j] = 0.0;
      elastic[i][j] = 0.0;
      viscous[i][j] = 0.0;
    }
  }

  //Initialize coordination number to 0
  for (i=0; i<N; i++){
    c[i] = 0;
  }

  mn = 50000;
  mx = -5;
  avgc = 0;

  //for every particle i....
  for (i=0; i<N; i++){
    if (flag[i] == 0)
      continue;

    //and its neighboring particles
    for (k=0; k<n_neighbor[i]; k++){
      j = n[i][k];

      if (j <= i)
        continue;

      if (flag[j] == 0)
        continue;

      dr2 = pbc_vdr ( x[i], x[j], dr);
      dij = ( d[i] + d[j] ) / 2.0;
      dij2 = dij * dij;

      if ( dr2 > dij2 )
        continue;

      rij = sqrt(dr2);
      
      //count coordination number when they overlap
      c[i] += 1;
      c[j] += 1;
      avgc += 2; 

      Utotal += ( 1.0 - 2.0*(rij/dij) + dr2/dij2 );
      F_pre = (2.0 * (1.0/(rij*dij) - 1.0/dij2));
      
      //Calculate prefactor for dissipation force (vi-vj) dot rij
      p=0.0;
      for (z=0; z<3; z++){
        dv[z] = v[i][z] - v[j][z];
        p += dv[z]*dr[z]/rij;
      }

      //loop over each cartesian coordinate, f_ij = -f_ji
      for (z=0; z<3; z++) {
        f[i][z] += (F_pre * dr[z]) - b*p*dr[z]/rij ;
        f[j][z] -= (F_pre * dr[z]) - b*p*dr[z]/rij ;


        //Calculate Total/Elastic/Viscous Stress, using the microscopic stress equation
        for(zz=0; zz<3; zz++){
          if (z==zz){
            s[z][zz] += -(N*T/vol) - (1.0/vol)*((F_pre*dr[z])-b*p*dr[z]/rij)*dr[zz];
            elastic[z][zz] += -(N*T/vol) - (1.0/vol)*(F_pre*dr[z])*dr[zz];
            viscous[z][zz] += -(N*T/vol) + (1.0/vol)*(b*p*dr[z]/rij)*dr[zz];

            // if we are looking at a neighbor of a tracer particle
            if (flag[i] == 2){
              stracer[i-(N-N_tracer)][z][zz] += -(N*T/vol) - (1.0/vol)*((F_pre*dr[z])-b*p*dr[z]/rij)*dr[zz];
            }
          }
          else {
            s[z][zz] += -(1.0/vol)*((F_pre*dr[z])-b*p*dr[z]/rij)*dr[zz];
            elastic[z][zz] += -(1.0/vol)*(F_pre*dr[z])*dr[zz];
            viscous[z][zz] += (1.0/vol)*(b*p*dr[z]/rij)*dr[zz];

            //if we are looking at a neighbor of a tracer particle
            if (flag[i] == 2){
              stracer[i-(N-N_tracer)][z][zz] += -(1.0/vol)*((F_pre*dr[z])-b*p*dr[z]/rij)*dr[zz];
            }
          }
        }
      }
    }//end of neighbor particle loop
    if (c[i] < mn)
      mn = c[i];
    if (c[i] > mx)
      mx = c[i];
  }//end of particle loop
  avgc = avgc / (N-n_flag);

  cdum = 0,0;
  ZZ = 0.0;
  // to calculate Z-Zc
  for (i=0; i<N; i++){
    if (flag[i] == 0) continue; //to exclude ghost particles
    if (c[i] == 0) continue; //to exclude rattlers
    ZZ += c[i];
    cdum += 1.0;
  }
  ZZ = ZZ/cdum; //averaged coordination number per particle
  Ziso = ZZ - ZC;  

  return Utotal;
}

double dgauss (double ave, double var){
  double x= ran2();
  double y= ran2();

  double dgauss = ave + var * sqrt (-2.0*log(x)) * cos(2.0*PI*y);

  return (dgauss);
}

/*==========*/
/* random.c */
/*==========*/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

// generates random number from [0,1]
double ran2 (void) {
  int j;
  long int k;
  static long idum2 = 123456789;
  static long iy =0;
  static long iv[NTAB];
  double temp;
  extern long idum;

  if (idum <= 0){
    if (-(idum) <1) idum=1;
    else idum = -(idum);
    idum2 = (idum);
    for (j=NTAB + 7; j>=0; j--){
      k=(idum)/IQ1;
      idum = IA1 * (idum - k * IQ1) - k * IR1;
      if (idum <0) idum += IM1;
      if (j<NTAB) iv[j] = idum;
    }

    iy= iv[0];
  }

  k = (idum) / IQ1;
  idum = IA1 * (idum - k * IQ1) - k * IR1;
  if (idum<0) idum += IM1;
  k = idum2 / IQ2;
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
  if (idum2< 0) idum2 += IM2;
  j = iy/NDIV;
  iy = iv[j] - idum2;
  iv[j] = idum;
  if (iy<1) iy += IMM1;
  if ((temp = AM*iy) > RNMX ) return RNMX;
  else return temp;

}

// Neighbor calculation routine that goes through N x N loop
void neighbor_N2 () {
  int i, j;
  double dr[3], dr2;

  for (i=0; i<N; i++){
    if (flag[i] == 0) continue;

    //initialize number of neighbor
    n_neighbor[i] = 0;

    for (j=0; j<N; j++){
      if (j == i) continue;
      if (flag[j] == 0) continue;
  
      dr2 = pbc_vdr ( x[i], x[j], dr);

      //if distance i-j > rnl, terminate iteration
      if (dr2 > rnl2)
        continue;

      //Assign j to identify neighbor
      n[i][n_neighbor[i]] = j;
            
      //add 1 to count the number of neighbor of particle i
      n_neighbor[i] += 1;
 
    }
  }
}


void neighbor () {
  int i, j, k;
  double dr[3], dr2;
  int cix, ciy, ciz, idx, idy, idz, ncix, nciy, nciz;

  //Count neighbors around particle i in the nearest +/-1 cells
  for (i=0; i<N; i++){
    if (flag[i] == 0) 
      continue;

    //initialize number of neighbor
    n_neighbor[i] = 0;
    
    //cell that i is in
    cix = a[i][0];
    ciy = a[i][1];
    ciz = a[i][2];

    //X
    for (idx=-1; idx<2; idx++){
      //# of cell we are currently in (which is +/- cell of cix)
      ncix = cix + idx ;
      
      //if we are in the 1st cell, we loop back to the last cell for -1 cell for PBC
      if (ncix <0) {
        ncix = n_cell[0] -1;
      }

      //if we are in the last cell, we loop back to the first cell for +1 cell for PBC
      if (ncix == n_cell[0]) {
        ncix = 0;
      }

      //Y
      for (idy=-1; idy<2; idy++){
        //# of cell we are currently in (which is +/- cell of cix)
        nciy = ciy + idy ;
        
        //if we are in the 1st cell, we loop back to the last cell for -1 cell for PBC
        if (nciy <0) {
          nciy = n_cell[1] -1;
        }

        //if we are in the last cell, we loop back to the first cell for +1 cell for PBC
        if (nciy == n_cell[1]) {
          nciy = 0;
        }
    
        //Z
        for (idz=-1; idz<2; idz++){
          //# of cell we are currently in (which is +/- cell of cix)
          nciz = ciz + idz ;
      
          //if we are in the 1st cell, we loop back to the last cell for -1 cell for PBC
          if (nciz <0) {
            nciz = n_cell[2] -1;
          }

          //if we are in the last cell, we loop back to the first cell for +1 cell for PBC
          if (nciz == n_cell[2]) {
            nciz = 0;
          }
  
          //Labeling j as the 'k'th particle in the cell we are in
          for (k=0; k<n_in_cell[ncix][nciy][nciz]; k++){
            j = cell[ncix][nciy][nciz][k];

            if (j == i) continue;
            if (flag[j] == 0) continue;

            // if j is less than i, terminate for loop and do not calculate interaction because we are looking at wrong particle
            //if (j <= i)
              //continue;
  
            dr2 = pbc_vdr ( x[i], x[j], dr);

            //if distance i-j > rnl, terminate iteration
            if (dr2 > rnl2)
              continue;

            //Assign j to identify neighbor
            n[i][n_neighbor[i]] = j;
            
            //add 1 to count the number of neighbor of particle i
            n_neighbor[i] += 1;
          }
        }
      }
    }
  }
}

// Allocate memory to arrays that will be used for the cell list
void cell_memory () {
  int i, j, k;
  int cell_max = N/2; 

  //Memory allocation for # of cells in each direction
  n_cell = (int*) calloc (3, sizeof (int));
  l_cell = (double*) calloc (3, sizeof (double));

  //Calculate number of cells in x,y,z
  //and using that calculate length of cell in x,y,z
  for (i=0; i<3; i++) {
    n_cell[i] = int (box[i] / rnl);
    l_cell[i] = box[i]/double (n_cell[i]);
  }

  //Define max number of neighbors you can have
  n_max = 20.0 * (N/vol) * (4.0/3.0) * PI * rnl * rnl * rnl;

  //Allocate memory for **n, identification of the neighbor
  n = (int**) calloc (N, sizeof (int*));
  for (i=0; i<N; i++){
    n[i] = (int*) calloc (n_max, sizeof (int));
  }

  //Allocate memory for *n, # of neighbors a particle can have
  n_neighbor = (int*) calloc (N, sizeof (int));

  //Memory allocation for **a, to store cell # particle i belongs to
  a = (int**) calloc ( N, sizeof (int*));
    for (i=0; i<N; i++){
      a[i] = (int*) calloc (3, sizeof (int));
    }

  //Allocate memory for # of particles in one cell
  n_in_cell = (int***) calloc (n_cell[0], sizeof (int**));
  for (i=0; i<n_cell[0]; i++){
    n_in_cell[i] = (int**) calloc (n_cell[1], sizeof(int*));
    for (j=0; j<n_cell[1]; j++){
      n_in_cell[i][j] = (int*) calloc (n_cell[2], sizeof (int));
    }
  }
  
  //Allocate memory for ****cell, identification of particle in cell 
  cell = (int****) calloc (n_cell[0], sizeof (int***));
    for (i=0; i<n_cell[0]; i++){
      cell[i] = (int***) calloc (n_cell[1], sizeof(int**));
      for (j=0; j<n_cell[1]; j++){
        cell[i][j] = (int**) calloc (n_cell[2], sizeof (int*));
        for (k=0; k<n_cell[2]; k++){
          cell[i][j][k] = (int*) calloc (cell_max, sizeof (int));
        }
      }
    }
}

//Assign particles to cell
void assign_cell() {
  int i, j, k;

  //Initialize # of particles in cell to 0
  for (i=0; i<n_cell[0]; i++){
    for (j=0; j<n_cell[1]; j++){
      for (k=0; k<n_cell[2]; k++){
        n_in_cell[i][j][k] = 0;
      }
    }
  }

  //Calculate which box particle i is in, and label particle i 
  for (i=0; i<N; i++){
    a[i][0] = int (( x[i][0] + boxh[0]) / l_cell[0] ); 
    a[i][1] = int (( x[i][1] + boxh[1]) / l_cell[1] ); 
    a[i][2] = int (( x[i][2] + boxh[2]) / l_cell[2] ); 

    for ( j=0 ; j<3 ; j++ )
      if ( a[i][j] < 0 || a[i][j] >= n_cell[j] ) {
        printf("Cell index error on particle %d in direction %d\n" , i , j ) ;
        printf("Found cell index %d, should be between 0 and %d\n" , a[i][j] , n_cell[j] );
        printf("position of particle: %lf %lf %lf\n" , x[i][0], x[i][1], x[i][2] ) ;
        print_net_diameter();
        exit(1);
      }

    cell[a[i][0]][a[i][1]][a[i][2]][n_in_cell[a[i][0]][a[i][1]][a[i][2]]] = i;

    n_in_cell[a[i][0]][a[i][1]][a[i][2]] += 1;
  }
}

void write_coarsen_xyz(){
  FILE *coarsen;
  int i;
  
  if (coarsen_step ==0)
    coarsen = fopen ("coarseningtraj.xyz", "w");
  else
    coarsen = fopen ("coarseningtraj.xyz", "a");

  // Writes number of particles and an extra blank line
  fprintf (coarsen, "%d\n\n", N);

  for (i=0; i<N; i++){
    fprintf (coarsen, "H %1.16lf %1.16lf %1.16lf %1.16lf\n", x[i][0], x[i][1], x[i][2], d[i]);
  }

  fclose (coarsen);
}

void write_tracer_xyz(){
  FILE *coarsen;
  int i;
  
  if (coarsen_step ==0)
    coarsen = fopen ("tracertraj.xyz", "w");
  else
    coarsen = fopen ("tracertraj.xyz", "a");

  // Writes number of tracer particles and an extra blank line
  fprintf (coarsen, "%d\n\n", N_tracer);

  for (i=0; i<N; i++){
    if (flag[i] == 0)
      continue;
    if (flag[i] == 1)
      continue;
    fprintf (coarsen, "H %1.16lf %1.16lf %1.16lf\n", x[i][0], x[i][1], x[i][2]);
  }

  fclose (coarsen);
}

void print_tracer_stress(){
  FILE *otp;
  int i;

  if (coarsen_step == 0)
    otp = fopen ("tracerstress.dat", "w");
  else 
    otp = fopen ("tracerstress.dat", "a");

  // Writes number of tracer particles and an extra blank line
  fprintf (otp, "%d\n\n", N_tracer);

  for (i=0; i<N_tracer; i++){
    fprintf (otp, "%d %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf\n", i+(N-N_tracer), stracer[i][0][0], stracer[i][1][1], stracer[i][2][2], stracer[i][0][1], stracer[i][0][2], stracer[i][1][2]);
  }

  fclose(otp);
}

void write_xyz(){
  FILE *otp;
  int i;

  //If this is the first step, erases old file by opening with w
  //for other steps, file is appended

  if (step ==0)
    otp = fopen ("traj.xyz", "w");

  else 
    otp = fopen ("traj.xyz", "a");

  //writes number of particles and an extra blank line
  fprintf (otp, "%d\n\n", N);

  //writes the particles as H atoms and He atoms
  for (i=0; i<N; i++){
    fprintf (otp, "H %1.16lf %1.16lf %1.16lf\n", x[i][0], x[i][1], x[i][2] );
   
  }
  fclose (otp);
}

void fire(){
  int i, j, z;
  int Ncount, doprint;
  double conv, alpha;
  double vel, force, P, Uo;
  double min, max, maxx, force2, s2, dr2, dr[3], Stotal;
  FILE *otp, *stp, *coordi, *otp2;
double check;
check = 0.0;
Stotal = 0.0;

  otp2 = fopen("EvsS_fire.dat","w");

  for (i=0;i<N;i++){
    rmsd[i] = 0.0;
  }

  if (coarsen_step == 0){
    otp = fopen ("fire.dat", "w");
  }
  if (coarsen_step > 0){
    otp = fopen ("fire.dat", "a");
  } 
  //stp = fopen ("coordination_fire.dat", "w");
  //coordi = fopen ("coordination_fire_perparti.dat", "w");

  vel = 0.0;
  conv = 0.0;
  force = 0.0;
  Uo = 0.0;
  b = 0.0;
  Ncount = 0;
  alpha = alphastart;
  min = max = 0.0;
    
  timesum = 0.0;
  fr_max2 = 0;
  nf = 0;
/*
  //before calculating, define initial coordinates 
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      pre_fin[i][j] = x[i][j];
    }
  }
*/ 
  for (step=0; step<firestep_max; step++){
    
    doprint = 1;

    
    //update neighbor list at given frequency for initial minimization
    if (step % neighbor_freq == 0 && fire_flag == 1){
      assign_cell();
      neighbor();
    }

    //update neighbor at regular frequency 
    else if (max > (0.5*rnl) || step % neighbor_freq2 == 0){
      if (rnl > max_cl){
        neighbor_N2();
        check = U_system();
        if (fabs(check-Uo)/double(N-n_flag) > 0.00001 && step != 0){
          printf ("Change in PE detected, neighbor list is faulty\n");
          printf ("check: %1.16lf, Uo: %1.16lf\n", check, Uo);
          exit(1);
        }
      }
      else {
        assign_cell();
        neighbor();
        check = U_system();
        if (fabs(check-Uo)/double(N-n_flag) >0.00001 && step != 0){
          printf ("Change in PE detected2, neighbor list is faulty\n");
          printf ("check: %1.16lf, Uo: %1.16lf\n", check, Uo);
          exit(1);
        }
      }
    }


    /*------------------------------------------------------------------------MD STEPS------*/
    // For every particle i ...
    // Update v(t + delt/2) and x(t + delt)

    max = -5.0;
    min = 50000;

    for (i=0; i<N; i++){
      for (j=0; j<3; j++){
        init[i][j] = x[i][j];
      }
    }

    for (i=0; i<N; i++){
      // If particle is flagged as 0, turn off the velocity
      if (flag[i] == 0){
        v[i][0]=v[i][1]=v[i][2] = 0.0;
      }
      // If particle is flagged as 0, continue onto next iteration
      if (flag[i] == 0) continue;

      rmsd[i] = 0.0;
      for (j=0; j<3; j++){

        v[i][j] = v[i][j] + delt * f[i][j] / 2.0;
        
        x[i][j] = x[i][j] + delt * v[i][j] ;

        rmsd[i] += (delt*v[i][j])*(delt*v[i][j]);

        //Ensure particles remain in the box
        if (x[i][j] > boxh[j])
          x[i][j] -= box[j];
        else if (x[i][j] < -boxh[j])
          x[i][j] += box[j];
      }

      if (rmsd[i] < min)
        min = rmsd[i];
      else if (rmsd[i] > max)
        max = rmsd[i];
    }//end of particle loop

    //Calculate new energy and forces from updated **x
    PE = U_system();

    //Using updated force, Update velocity to v(t + delt)
    for (i=0; i<N; i++){
      if (flag[i] ==0)
        continue;
      for (j=0; j<3; j++){
        v[i][j] = v[i][j] + delt * f[i][j] / 2.0;
      }
    }

    //Calculate the temperature and KE based on the new velocities
    KE = 0.0;
    for (i=0; i<N; i++){
      if (flag[i] ==0) continue;
      for (j=0; j<3; j++){
        KE += v[i][j] * v[i][j];
      }
    }
    KE *= 0.5;

    T = 2.0 * KE / ( 3.0 * double(N-n_flag) );

/*
    //accumulate max force for convergence criterion
    maxx = -5.0;
    for (i=0; i<N; i++){
      force2 = 0.0;
      for (j=0; j<3; j++){
        force2 = f[i][j] * f[i][j];
      }
      if (force2 > maxx) maxx = force2;
    }
    
    if (maxx < 0.0) {
      printf("max force calc is incorrect! negative max value\n");
      exit(1);
    }
    maxx = sqrt(maxx);
    if (maxx < 0.0000001) break;
*/
    conv = fabs(PE-Uo);

    if (conv < 0.00000000000001){
      break;
    }

    Uo = PE;
    /*-----------------------------------------------------------------------END OF MD------*/

    /*--------------------------------------------------------------------FIRE ALGORITHM----*/
    // adjusting v, delt, and alpha
    // doesn't change PE
    //  still in step loop ...

    P = 0.0;

    // for every particle
    for (i=0; i<N; i++){
      if (flag[i] == 0)
        continue;
      //calculate power and magnitude of velocity and force on particle i
      for (z=0; z<3; z++){
        P += v[i][z] * f[i][z];
        vel += v[i][z] * v[i][z];
        force += f[i][z] * f[i][z];
      }
    }//end of particle loop

    //sqrt to get magnitude
    vel = sqrt(vel);
    force = sqrt(force);

    //for every particle, update velocity according to EOM of FIRE algorithm
    for (i=0; i<N; i++){
      if (flag[i] == 0) continue;
      for (z=0; z<3; z++){
        v[i][z] = (1.0-alpha) * v[i][z] + alpha * vel * f[i][z] / force;
      }
    }

    //if we are accelerating (going downhill)..
    if (P>0.0){

      //count how many steps we went last since power was negative (went uphill)
      Ncount += 1;

      //if number of steps since last uphill is greater than some min value,
      if (Ncount > Nmin){
        //increase the time step
        delt = delt * finc;

        //put a cap on how high the time step could be
        if (delt > deltmax){
          delt = deltmax;
        }

        //decrease alpha ("fraction of current velocity info?")
        alpha = alpha * falpha;

        //initialize Ncount back to 0
        Ncount = 0;

      }//end of Ncount > Nmin
    }//end of P>0

    //if we are going uphill, or at a minimum!
    else if (P<=0.0){
      //decrease timestep and initialize alpha
      delt = delt * fdec;
      alpha = alphastart;

      //kill momentum by setting v=0.0
      for (i=0; i<N; i++){
        if (flag[i] == 0) continue;
        for (z=0; z<3; z++){
          v[i][z] = 0.0;
        }
      }
      Ncount = 0;
      doprint = 0;
    }//end of P<=0

//    printf ("%lf\n", delt);

    /*--------------------------------------------------------------------------------------*/

  //Store final coordinates for MSD calculation
  if (doprint == 1){
    for (i=0; i<N; i++){
      for (j=0; j<3; j++){
        fin[i][j] = x[i][j];
      }
    }
  }
  FILE *rstep;
  if (step ==0 ) rstep = fopen ("del_r.dat","w");
  if (step !=0 ) rstep = fopen ("del_r.dat","a");
  meandistance=calc_MSD();
  if (doprint == 1){
    fprintf(rstep, "%d %1.16lf\n", step, meandistance);
  }
  fclose(rstep);


    if (step % print_freq == 0 && doprint == 1){
      nf += 1;
      fr_max2 += 1;
      //nf +=1;
      //fprintf (otp, "%d %lf %e %e %e\n" , step, KE+PE, s[2][0], elastic[2][0], viscous[2][0]);
      //fprintf (stp, "%d %d %d\n", avgc, mx, mn);
      printf ("Step: %d, Total E: %1.12lf, PE: %1.12lf, KE: %lf, T: %lf, szx: %e, sxz: %e\n", step, KE+PE, PE, KE, T, s[2][0], s[0][2] );

      s2 = 0.0;
      for (i=0; i<N; i++){
	      if (flag[i] ==0) continue;
	      dr2 = pbc_vdr(x[i], init[i], dr);
	      s2 += dr2;
      }
      Stotal += sqrt(s2);

      //calculate new PE
      PE = U_system();

      fprintf(otp2, "%1.16lf %1.16lf\n", Stotal, PE);




    }

    //save step number in which initilization ended
//    start = step;
 
double dummy;
dummy=U_system();
if (Uo != dummy){
printf ("PE calling is numerical error \n");
exit(1);
}
  
  }//End of step loop
/*
  //Store final coordinates for MSD calculation
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      fin[i][j] = x[i][j];
    }
  }
*/
  printf ("FIRE: step:%d PE:%1.16lf\n", step, PE);

  if (coarsen_step >= 0){
    fprintf (otp, "Coarsen step: %d FIRE step: %d PE:%1.16lf N: %d\n", coarsen_step, step, PE, (N-n_flag));
    fclose(otp);

    if (step <= 2)
      rcount +=1;
    else if(step>10) 
      rcount = 0;
  }
printf("rcount: %d\n", rcount);

  if (rcount > 50){
    printf("50 steps of no rearrangement! Terminating simulation\n");
print_net_diameter();
print_pov_info();
    exit(1);
  }
  /*
  for (i=0; i<N; i++){
    fprintf (coordi, "%d\n", c[i]);
  }
  */
//  start += 1;
  //fclose (stp);
  //fclose (coordi);

}//End of FIRE function

void print_iso(){
  FILE *otp;
  int i;

  if (coarsen_step == 0) 
    otp = fopen("Z.dat", "w");
  else 
    otp = fopen("Z.dat", "a");

  fprintf (otp, "%d %lf %lf\n", coarsen_step, Ziso, cdum);
  
  fclose(otp);
}

void print_tracer_coordination(){
  FILE *otp;
  int i;

  if (coarsen_step == 0 ){
    otp = fopen("tracercoordi.dat", "w");
  }
  else {
    otp = fopen("tracercoordi.dat", "a");
  }

  fprintf (otp, "%d ", coarsen_step);
  for (i=(N-N_tracer); i<N; i++){
    fprintf (otp, "%d ", c[i]);
  }
  fprintf(otp, "\n");

  fclose(otp);
}

//prints coordination number
void print_coordi(){
  FILE *otp;
  int i;
  otp = fopen("coordination.dat","w");

  fprintf (otp, "step: %d\n", step);
  for (i=0; i<N; i++){
    fprintf (otp, "%d\n", c[i]);
  }

  fclose(otp);
}

//prints diameter at a particular coarsening step
// "diameter_#.dat" where # is the particular coarsening step
void print_diameter(){
  int i;
  FILE *otp;
  char out[40];

  //store the file name as a string
  sprintf (out, "diameter_%d.dat", coarsen_step);

  //writes a file named as the string
  otp = fopen ( out ,"w");

  for (i=0; i<N; i++){
    fprintf (otp, "%1.16lf\n", d[i]);
  }
  fclose(otp);
}

void print_initial_diameter(){
  int i;
  FILE *otp;

  otp = fopen ("initial_diameter.dat","w");

  for (i=0; i<N; i++){
    fprintf (otp, "%1.16lf\n" , d[i]);
  }
  fclose(otp);
}

void print_Pr(){
  int i;
  FILE *otp;
  char out[40];

  if (coarsen_step == 0)
    otp = fopen ("Pr.dat", "w");
  else 
    otp = fopen ("Pr.dat", "a");

  fprintf (otp, "%d %lf\n", coarsen_step, Pr);
  fclose (otp);

}

void print_msd () {
  int i;
  FILE *otp;
  char out[40];
                  
  //store the file name as a string
  sprintf (out, "msd_%d.dat", coarsen_step);

  //writes a file named as the string
  otp = fopen (out ,"w");
 
  for (i=0; i<N; i++){
    if (flag[i] ==0) continue;
    fprintf (otp, "%lf\n", msd[i]);
  }
  fclose(otp);
}

//prints diameter at a particular coarsening step
// "net_diameter_#.dat" where # is the particular coarsening step
void print_pov_info(){
  int i;
  FILE *otp;
  char out[40];

  //store the file name as a string
  sprintf (out, "pov_info_%d.dat", coarsen_step);

  //writes a file named as the string
  otp = fopen ( out ,"w");

  for (i=0; i<N; i++){
    if (flag[i]== 0) continue;
fprintf(otp, "%d\n\n", N-n_flag);
    fprintf (otp, "%lf %lf %lf %1.16lf\n", x[i][0], x[i][1], x[i][2], d[i]);
  }
  fclose(otp);
}


//prints diameter at a particular coarsening step
// "net_diameter_#.dat" where # is the particular coarsening step
void print_net_diameter(){
  int i;
  FILE *otp;
  char out[40];

  //store the file name as a string
  sprintf (out, "net_diameter_%d.dat", coarsen_step);

  //writes a file named as the string
  otp = fopen ( out ,"w");

  for (i=0; i<N; i++){
    if (flag[i]){
      fprintf (otp, "%1.16lf\n", d[i]);
    }
  }
  fclose(otp);
}

void print_coordinates(){
  FILE *otp;
  int i;
  char out[40];

  //store the file name as a string
  sprintf (out, "coordinates_%d.dat", coarsen_step);

  //writes a file named as the string
  otp = fopen ( out ,"w");

  for (i=0; i<N; i++){
    fprintf (otp, "%lf %lf %lf\n", x[i][0], x[i][1], x[i][2]);
  }
  fclose(otp);
}

//MSD Caculation
// Accumulate MSD for all particles
// before calling calc_MSD for first time, have to initially store pre_fin
double calc_MSD (){
  int i, j, count;
  double dr2 = 0.0;
  double MSD_total, dr[3], top, bot;

  Pr = MSD_total = top = bot = 0.0;
  for(i=0;i<N;i++){
      msd[i] = 0.0;
  }
/*
  //before calculating, define final coordinates 
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      fin[i][j] = xxx[nf-1][i][j];
    }
  }

  //before calculating, define initial coordinates 
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      pre_fin[i][j] = init[i][j];
    }
  }
    
*/

  count = 0;
  // accumulate r2 for every particle
  for (i=0; i<N; i++){
    if (flag[i] == 0) continue;

    dr2 = MSD (pre_fin[i], fin[i], dr);
    //dr2 = sqrt(dr2);
    msd[i] = dr2;
    MSD_total += dr2;
    count+=1;
  }

  //normalize MSD by number of particles
  MSD_total = MSD_total;

  // calculate participation ratio
  for (i=0; i<N; i++){
    if (flag[i] == 0) continue;
    top += msd[i];
    bot += msd[i]*msd[i];
  }
  top = top*top;
  Pr = top / bot / (N-n_flag);

  // store current coordinate to use for initial coordinate next run
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      pre_fin[i][j] = fin[i][j];
    }
  }

  return MSD_total;
}

//routine to calculate original and post-processed end trajectory
// to see if FIRE went down similar landscape
double calc_MSD2 (){
  int i, j, count;
  double dr2 = 0.0;
  double MSD_total, dr[3], top, bot;

  Pr = MSD_total = top = bot = 0.0;
  for(i=0;i<N;i++){
      msd[i] = 0.0;
  }

  //before calculating, define final coordinates 
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      fin[i][j] = xx[input+1][i][j];
    }
  }

  //before calculating, define initial coordinates 
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      pre_fin[i][j] = x[i][j];
    }
  }
    

  count = 0;
  // accumulate r2 for every particle
  for (i=0; i<N; i++){
    if (flag[i] == 0) continue;

    dr2 = MSD (pre_fin[i], fin[i], dr);
    //dr2 = sqrt(dr2);
    msd[i] = dr2;
    MSD_total += dr2;
    count+=1;
  }

  //normalize MSD by number of particles
  MSD_total = MSD_total;
  
  return MSD_total;

}


// Return distance between fin and pre_fin
double MSD (double pre_fin[3], double fin[3], double dr[3]){
  int i;
  double mdr2 = 0.0;

  for (i=0; i<3; i++){
    dr[i] = fin[i] - pre_fin[i];

    if (dr[i] > boxh[i]){
      dr[i] -= box[i];
    }
    else if (dr[i] < -boxh[i]) {
      dr[i] += box[i];
    }

    mdr2 += dr[i] * dr[i];
  }
  return mdr2;
}

void hessian_memory () {
  int i;

  H = (double**) calloc (3*N, sizeof (double*));
  for (i=0; i<3*N; i++){
    H[i] = (double*) calloc (3*N, sizeof (double));
  }

  minx = (double**) calloc (N, sizeof (double*));
  minf = (double**) calloc (N, sizeof (double*));
  ff = (double**) calloc (N, sizeof (double*));
  fb = (double**) calloc (N, sizeof (double*));
  for (i=0; i<N; i++){
    minx[i] = (double*) calloc (3, sizeof (double));
    minf[i] = (double*) calloc (3, sizeof (double));
    ff[i] = (double*) calloc (3, sizeof (double));
    fb[i] = (double*) calloc (3, sizeof (double));
  }
}

// Calculates forward force for hessian calculation
void forward_force(int i, int z, double delx) {
  int ii, jj;

  // initializes forward force 
  for (ii=0; ii<N; ii++){
    for (jj=0; jj<3; jj++){
      ff[ii][jj] = 0.0;
    }
  }
  
  // displaces coordinates by (+) delta x
  x[i][z] = x[i][z] + delx;
  
  //apply PBC
  if (x[i][z] < -boxh[z])
    x[i][z] += box[z];
  else if (x[i][z] > boxh[z])
    x[i][z] -= box[z];
  
  // Calculate the force of the new configuration
  PE = U_system();
  
  // store the force to **ff and restore coordinates to original minimized coordinates
  for (ii=0; ii<N; ii++){
    for (jj=0; jj<3; jj++){
      ff[ii][jj] = f[ii][jj];
      x[ii][jj] = minx[ii][jj];
    }
  }
}

// Calculates backward force for hessian calculation
void backward_force(int i, int z, double delx) {
  int ii, jj;

  // initializes backward force 
  for (ii=0; ii<N; ii++){
    for (jj=0; jj<3; jj++){
      fb[ii][jj] = 0.0;
    }
  }
  
  // displaces coordinates by (-) delta x
  x[i][z] = x[i][z] - delx;

  //apply PBC
  if (x[i][z] < -boxh[z])
    x[i][z] += box[z];
  else if (x[i][z] > boxh[z])
    x[i][z] -= box[z];
  
  // Calculate the force of the new configuration
  PE = U_system();
  
  // store the force to **fb and restore coordinates to original minimized coordinates
  for (ii=0; ii<N; ii++){
    for (jj=0; jj<3; jj++){
      fb[ii][jj] = f[ii][jj];
      x[ii][jj] = minx[ii][jj];
    }
  }
}

void hessian_num () {
  int i, j, z, zz;
  double delx = 0.00001;

  // Stores coordinates and forces at minimum to separate arrays: **minx and **minf
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      minx[i][j] = x[i][j];
      minf[i][j] = f[i][j];
    }
  }

  // start of upper half hessian calculation
  for (i=0; i<N; i++){ //particle indices for the rows
    for (z=0; z<3; z++){
      // Calculates forward and backward force for hessian calculation
      forward_force (i,z,delx);
      backward_force (i,z,delx);
      for (j=i; j<N; j++){ //particle indices for the columns
        for (zz=0; zz<3; zz++){
          H[3*i+z][3*j+zz] = (fb[j][zz] - ff[j][zz]) / (2.0*delx);
        }
      }
    } 
  }

  // using symmetry of Hessian, calculate lower half
  for (i=0; i<N-1; i++){
    for (j=i+1; j<N; j++){
      for (z=0; z<3; z++){
        for (zz=0; zz<3; zz++){
          H[3*j+zz][3*i+z] = H[3*i+z][3*j+zz];
        }
      }
    }
  }
}

// Original coarsening + constant flux (reservoir) method
void coarsening () {
  int i, j, k, closest_n;
  double aa, bb, cc2, dd;
  double dr[3], dr2, rij, dij, dij2;
  double volume, volume2, delv, diam, max;
  double area, delP;

  //Initialize area of overlap to 0
  for (i=0; i<N; i++){
    Q[i] = 0.0;
  }

  //For all i-j interactions...
  for (i=0; i<N; i++){
    // If the particle is flagged (0 volume), skip to next loop iteration
    if (flag[i] == 0 ) 
      continue;
    if (flag[i] == 2 ) 
      continue;

    for (k=0; k<n_neighbor[i]; k++){

      j = n[i][k];

      if (j<i)
        continue;

      // If the particle is flagged (0 volume), skip to next loop iteration
      if (flag[j] == 0)
        continue;
      if (flag[j] == 2)
        continue;

      dr2 = pbc_vdr(x[i], x[j], dr);
      dij = (d[i] + d[j])/2.0;
      dij2 = dij * dij;

      if (dr2 >= dij2)
        continue;

      rij = sqrt(dr2);

      bb = -(d[i]*d[i] - d[j]*d[j])/(rij*8.0) + rij/2.0;
      aa = rij - bb;
      cc2 = (d[i]*d[i]/4.0) - (aa*aa);

      //Accumulate total area of overlap for all i-j interactions
      area = cc2*PI;
      delP = (2.0/d[i] - 2.0/d[j]);
      Q[i] -= (C * area * delP);
      Q[j] += (C * area * delP);
    }
  }//End of particle loop

  /*----- ATTACH RESERVOIR METHOD ------*/
  // Caclulate <R> where R is particle radius
  // Sum (R) = dd   <R> = dd/(N-n_flag)
  dd = 0.0;
  for (i=0; i<N; i++){    // Loops all non-zero diameter particles (n_flag initially 0) 
    if (flag[i] == 0 ) continue;
    if (flag[i] == 2) continue;
    dd += d[i]/2.0;
  }
  // Average the value
  dd = dd/double(N-n_flag-N_tracer);

  // Calculate the radius of reservoir
  r_res = dd; 

  for (i=0; i<N; i++){
    max = 50000.0;
    closest_n = -1;
    if (flag[i] == 0)
      continue;
    if (flag[i] == 2)
      continue;

    // calculate volume change using EOM: delv = -b (1/ri - 1/r_res) and sum delv = 0
    delv = -Cf*(2.0/d[i] - 1.0/r_res)*(d[i]/2.0);

    // calculate new volume of particle i
    volume = (d[i]*d[i]*d[i]*PI/6.0) + delv + Q[i];

    // IF the volume goes down to 0 or below,
    if (volume <= 0.0){
      // set d[i] to 0
      d[i] = 0.0;
      // flag it as "ghost" and add 1 to number of ghost particles
      flag[i]=0;
      n_flag += 1;

      // most likely the volume would not be 0.0, but a negative number, so calculate how much (-) it is
      // and the corresponding diameter
      delv = -volume;
        
      // identify i's closeset neighbor 
      for (k=0; k<n_neighbor[i]; k++){
        j=n[i][k];
        if (flag[j] == 0)
          continue;
        if (flag[j] == 2) continue;
        dr2 = pbc_vdr(x[i], x[j], dr);
        if (dr2 < max && delv < d[j]*d[j]*d[j]*PI/6.0){
          closest_n = j;
          max = dr2;
        }
      }

      // checking to see if neighbors were found
      if (closest_n == -1){
        printf ("No neighbors were found for particle %d! Terminating simulation\n", i);
        exit(1);
      }

      // calculate the volume of the closest neighbor and subtract how much (-) we went for the volume of i
      volume2 = (d[closest_n]*d[closest_n]*d[closest_n]*PI/6.0) - delv;
      // update the diameter for the closest  neighbor
      d[closest_n] = pow ((volume2*6.0/PI), 1.0/3.0);

      if (d[closest_n] < 0.0){
        printf("Negative diameter!\n");
        printf("particle %d(j) became negative after interacting with particle %d(i)\n", closest_n, i);
        printf("d[j]: %lf, d[i]: %lf\n",d[closest_n], d[i]);
        exit(1);
      }
    }// end of negative volume loop

    else{ 
      d[i] = pow ( (volume*6.0/PI), 1.0/3.0);
    }

  }//end of particle loop

}//end of coarsening

void mass_consv() {
  int i;
  double drift;
  double add, sum;

  drift = add = sum = 0.0;

  // Accumulate the total volume
  for (i=0; i<N; i++){
    if (flag[i] == 0) continue;
    if (flag[i] == 2) continue;
    sum += d[i]*d[i]*d[i]*PI/6.0;
  }

  // Calculate "drifted" volume
  drift = inivol - sum;

  // divide by number of particles present
  drift = drift / double(N-n_flag-N_tracer);

  // update diameter to account for the drift
  for (i=0; i<N; i++){
    if (flag[i] == 0) continue;
    if (flag[i] == 2) continue;
    
    add = d[i]*d[i]*d[i]*PI/6.0 + drift;
    d[i] = pow ( add*6.0/PI, 1.0/3.0);
  }
}

// ADAPTIVE neighbor list
void rnl_update() {
  int i;
  double maxx = -5.0; //variable to store max diameter

  // Determine MAX diameter of the system
  for (i=0; i<N; i++){
    if (flag[i] == 0)
      continue;
 
    if (d[i] > maxx){
      maxx = d[i];
    }
  }//end of determing max diam

  if (maxx > 0.95*max_cl*3.0){
    printf ("max diameter is 95 percent of box length!\n");
    printf ("max diameter: %lf\n", maxx);
    print_net_diameter();
    exit(1);
  }

  if (rnl < 1.2 * maxx){
    // Re-define rnl if rnl is smaller than 1.2 * max diameter 
    rnl = 1.2 * maxx;
    rnl2 = rnl * rnl;

  }// end of updating
}//end of routine

void print_max_diam () {
  int i, j;
  double maxx = -5.0;
  FILE *otp;

  // Determine MAX diameter of the system
  for (i=0; i<N; i++){
    if (flag[i] == 0)
      continue;
 
    if (d[i] > maxx){
      maxx = d[i];
    }
  }//end of determing max diam

  if (coarsen_step == 0){
    otp = fopen ("max_diameter.dat", "w");
  }

  else {
    otp = fopen ("max_diameter.dat", "a");
  }
  
  fprintf (otp, "%lf\n", maxx);
  fclose(otp);
}


//reads entire trajectory, and stores the coordinates of (input-2) into separate array that i can work to minimize
void read_input () {
    int i;
    FILE *inp, *inp2, *otp;
    
    inp = fopen("coarseningtraj.xyz","r");
    inp2 = fopen("velocity.dat","r");
    
    if (inp == NULL || inp2 == NULL){
        printf("No input file detected!\n");
        exit(1);
    }
    //otp = fopen("test.xyz","w");
    
    nf2 = 0;    
    frame = 0;
    
    //while it is not the end of the file
    while(!feof(inp)){
        //scan for the number of particles
        fscanf (inp, "%d %lf\n\n", &N, &vol);
        //fprintf (otp, "%d\n\n", N);
        
        //then scsan for the coordinates of N atoms in that frame
        for (i=0;i<N;i++){
            fscanf(inp, "H %lf %lf %lf %lf\n", &xx[frame][i][0], &xx[frame][i][1], &xx[frame][i][2], &dd[frame][i]);
            //fprintf(otp, "H %lf %lf %lf %lf\n", x[frame][i][0], x[frame][i][1], x[frame][i][2], d[i]);
        }
        
        //add 1 to frame to store next set of coordinates to next frame
        frame += 1;
        nf2 += 1;
    }
    
    fclose(inp);
     

    //while it is not the end of the file
    while(!feof(inp2)){
        //scan for the number of particles
        fscanf (inp2, "%d\n\n", &N);
        //fprintf (otp, "%d\n\n", N);
        
        //then scsan for the coordinates of N atoms in that frame
        for (i=0;i<N;i++){
            fscanf(inp2, "H %lf %lf %lf %lf\n", &vv[frame][i][0], &vv[frame][i][1], &vv[frame][i][2], &dd[frame][i]);
            //fprintf(otp, "H %lf %lf %lf %lf\n", x[frame][i][0], x[frame][i][1], x[frame][i][2], d[i]);
        }
        
    }
    
    fclose(inp2);  
}

//reads entire trajectory, and stores the coordinates of (input-2) into separate array that i can work to minimize
void read_rnl() {
    int i;
    FILE *inp, *otp;
    double *nrl;

    inp = fopen("rnl.dat","r");
    
    if (inp == NULL){
        printf("No input file detected!\n");
        exit(1);
    }
   
    nrl = (double*) calloc (nf2, sizeof (double));

    for (i=0; i<nf2; i++){
      fscanf(inp, "%lf\n", &nrl[i]);
    }
    
    fclose(inp);
    //fclose(otp);

    rnl = nrl[input];
    rnl2 = rnl * rnl;
}


void store_coordinate() {
    int i, j, k;
    int dummy=0;
    n_flag = 0;
    
    
    for (i=0; i<N; i++){
        x[i][0] = xx[input][i][0];
        x[i][1] = xx[input][i][1];
        x[i][2] = xx[input][i][2];

        v[i][0] = vv[input][i][0];
        v[i][1] = vv[input][i][1];
        v[i][2] = vv[input][i][2];
        
        init[i][0] = xx[input][i][0];
        init[i][1] = xx[input][i][1];
        init[i][2] = xx[input][i][2];
        
        d[i] = dd[input][i];
        
        //count "vanished" particles
        if (d[i] == 0.0){
            flag[i] = 0;
            n_flag += 1;
        }
        
        //count tracers
        if (d[i] == 2.5){
            flag[i] = 2;
            dummy += 1;
        }
    }
    
    if (dummy != N_tracer) {
        printf("Incorrectly calculated number of tracer particles!\n");
        printf("tracer: %d sim value: %d\n", N_tracer, dummy);
        exit(1);
    }
    
    FILE *otp, *otp2;
    otp = fopen("test.xyz","w");
    otp2 = fopen("test2.xyz","w");
    for (i=0; i<N; i++){
        fprintf(otp, "%1.16lf %1.16lf %1.16lf %1.16lf %d\n", x[i][0], x[i][1], x[i][2], d[i], flag[i]);
        fprintf(otp2, "H %1.16lf %1.16lf %1.16lf %1.16lf\n", x[i][0], x[i][1], x[i][2], d[i]);
    }
    fclose(otp);
    fclose(otp2);

    FILE *flagg;
    flagg = fopen("flag.dat","w");
    for (i=0; i<N; i++){
      fprintf(flagg, "%d\n", flag[i]);
    }
    fclose(flagg);
}


//this reads the trajectory file generated from the minimization from this step
void read_traj() {
    int i, j;
    FILE *inp, *otp, *inp2;
    double dummmy;
    

    xxx = (double***) calloc(fr_max2, sizeof (double**));
    for (i=0; i<fr_max2; i++){
        xxx[i] = (double**) calloc(N, sizeof(double*));
        for (j=0; j<N; j++){
            xxx[i][j] = (double*) calloc(3, sizeof(double));
        }
    }
    cl_frame = (double*) calloc (fr_max2, sizeof (double));
    Usys = (double*) calloc (fr_max2, sizeof (double));
    Udiff = (double*) calloc (fr_max2, sizeof (double));

    
    
    inp = fopen("traj.xyz","r");
    inp2 = fopen("energy.dat","r");
    
    printf("check3\n"); fflush(stdout);
    
    
    if (inp == NULL || inp2 == NULL){
        printf("No input file detected!\n");
        exit(1);
    }
    //otp = fopen("test.xyz","w");
    
    frame = 0;
    
    //SCAN TRAJECTORY
    //while it is not the end of the file
    while(!feof(inp)){
        //scan for the number of particles
        fscanf (inp, "%d\n\n", &N);
        
        //then scsan for the coordinates of N atoms in that frame
        for (i=0;i<N;i++){
            fscanf(inp, "H %lf %lf %lf\n", &xxx[frame][i][0], &xxx[frame][i][1], &xxx[frame][i][2]);
        }
        
        //add 1 to frame to store next set of coordinates to next frame
        frame += 1;
    }

    fclose(inp);
//    printf("scanned trajectory\n"); fflush(stdout);

    frame = 0;
    //SCAN ENERGY 
    //while it is not the end of the file
    while(!feof(inp2)){
        //scan for the number of particles
        fscanf (inp2, "%lf %lf", &dummmy, &Usys[frame]);
        
//    printf("read %d energy\n",frame); fflush(stdout);
        //add 1 to frame to store next set of coordinates to next frame
        frame += 1;
    }
    
    fclose(inp2);
    
//    printf("read all energy\n"); fflush(stdout);
/*
    FILE *test;
    test = fopen("testread.dat","w");
    for (frame = 0; frame<nf; frame++){
      fprintf(test, "%d\n\n",N);
      for (i=0; i<N; i++){
        fprintf(test, "%1.16lf %1.16lf %1.16lf\n", xxx[frame][i][0], xxx[frame][i][1], xxx[frame][i][2]);
      }
    }
    fclose(test);
*/
}

void remove_PBC(){
    int i, j, frame;
    double dr[3];
    
    //printf("nf: %d\n", nf); fflush(stdout);

    for (frame=1; frame< fr_max2; frame++){
        for (i=0; i<N; i++){
            for (j=0; j<3; j++){
                
                dr[j] = xxx[frame][i][j] - xxx[frame-1][i][j];
                
                if (dr[j] > boxh[j]){
                    xxx[frame][i][j] -= box[j];
                }
                else if (dr[j] < -boxh[j]){
                    xxx[frame][i][j] += box[j];
                }
            }
        }
    }
}

void steepest_descent(){
  int i, j, k;
  double force, maxx, dr2, dr[3], s2, Stotal;
  FILE *otp;

  otp = fopen("EvsS.dat","w");

  //lambda = 0.002;
  force = 0.0;
  fr_max2 = 0;
  nf = 0;
  Stotal = 0.0;

  //get current force info
  PE = U_system();
  for (step=0; step<firestep_max; step++){
    if (step % neighbor_freq == 0 ){
      neighbor_N2();
    }

    //before SD step, store coordinates for dr2 calc
    for (i=0; i<N; i++){
      for (j=0; j<3; j++){
        init[i][j] = x[i][j];
      }
    }

    for (i=0; i<N; i++){
      if (flag[i] == 0) continue;
      for (j=0; j<3; j++){
        x[i][j] += delt * f[i][j];

        //apply PBC
        if (x[i][j] > boxh[j])
          x[i][j] -= box[j];
        else if (x[i][j] < -boxh[j])
          x[i][j] += box[j];
      }
    }

    s2 = 0.0;
    for (i=0; i<N; i++){
      if (flag[i] ==0) continue;
      dr2 = pbc_vdr(x[i], init[i], dr);
      s2 += dr2;
    }
    Stotal += sqrt(s2);

    //calculate new PE
    PE = U_system();

    fprintf(otp, "%1.16lf %1.16lf\n", Stotal, PE);

    //accumulate max force for convergence criterion
    maxx = -5.0;
    for (i=0; i<N; i++){
      force = 0.0;
      for (j=0; j<3; j++){
        force += f[i][j]*f[i][j];
      }
      if (force > maxx) maxx = force;
    }

    if(maxx < 0.0){
      printf("max force calculation is incorrect! negative max value\n");
      exit(1);
    }

/*
    if (step % print_freq == 0 ){
      printf("step: %d PE: %1.16lf\n", step, PE);
      nf += 1;
      fr_max2 += 1;
      write_xyz();
      write_energy();
    }
*/
    maxx = sqrt(maxx);
    if (maxx < 0.0000001) break; 

  }//end of step loop
  fclose(otp);
    
}

//regular MSD routine that calculates MSD per lag time
void MSD() {
  int i, j, k;
  char tt[80];


  ///////////////////
  //  MSD routine  //
  double MSD_total, dr2, dr[3], delx, dely, delz, delx2, dely2, delz2, delr2;
  int count, tau;
  //contour length - cl, end to end R
  double cl_total, Rend_total;

    
  FILE *msd, *msd2, *msd3;
  msd = fopen ( "MSD.dat" , "w" );
  msd2 = fopen ( "MSD_xyz.dat" , "w" );
  msd3 = fopen ( "Rend_Cl.dat" , "w" );
  
  /*
  for (i=0; i<N; i++){
    cl[i] = 0.0;
  }

  for (i=0; i<nf; i++){
    cl_frame[i] = 0.0;
  }
  */
    
  FILE *otp2;
  FILE *lagtest;
  for (tau = 1; tau < fr_max2; tau++){
    MSD_total = 0.0;
    count = 0;
    cl_total = 0.0;
    
    //to print separate lag time values
    sprintf(tt, "tau_%d.dat", tau);
        
    //MSD loop for each frame of tau intervals
    for (step=tau; step < fr_max2; step++){
      Rend_total = 0.0;
      cl_total = 0.0;
     
/*
      for (i=0; i<N; i++){
        cl[i] = 0.0;
      }
*/ 

      //for all particles
      for (i=0; i<N; i++){

        if (flag[i] ==0) continue;

        dr2 = 0.0;
        for (j=0; j<3; j++){
                    
          dr[j] = xxx[step][i][j] - xxx[step-tau][i][j];
                    
    
          //this is like counting the Rend of particles in adjacent tau frames
          dr2 += dr[j] * dr[j];
        }//end of j (dimension)
/*
        //get contour length of particle i in separate array
        if (tau == 1){
          cl[i] = sqrt(dr2); // this is contour length of N bonds. (i (step+1) - i(step) -> N # of i's)
          if (cl[i] < 0.0){
            printf("cl[i] is negative for particle %d!\n", i);
            exit(1);
          }
        }
*/
        //if (step == nf-1 ){
          Rend[i] = dr2;
          Rend[i] = sqrt(Rend[i]);

          if (Rend[i] < 0.0){
            printf("negative end to end distance for particle %d at step %d for tau of %d\n", i, step, tau);
            exit(1);
          }
        //}

                
        MSD_total += dr2;
                
        //accumulate delx, dely, .. etc
        delx += dr[0];
        dely += dr[1];
        delz += dr[2];
        delx2 += dr[0] * dr[0];
        dely2 += dr[1] * dr[1];
        delz2 += dr[2] * dr[2];
        count += 1;
                
        }//end of i
/*
      //get contour length of frame
      if (tau == 1){
        //accumulate contour length for frame = (sum of contour length of i's between two [step] and [step-tau] frames)
        for (i=0; i<N; i++){
          cl_frame[step-tau] += cl[i];
          //printf("cl_frame[%d] = %1.16lf\n", step-tau, cl_frame[step-tau]);
          if (cl_frame[step-tau] < 0.0){
            printf("cl_frame went negative for %d frame\n", step-tau);
            exit(1);
          }
        }
      }
*/
      //accumulate total Rend
      for (i=0; i<N; i++){

        if (flag[i] == 0) continue;

        Rend_total += Rend[i];
        if (Rend_total < 0.0){
          printf("Rend_total is negative at step %d and tae %d\n", step, tau);
          exit(1);
        }
      }

      // accumulate total contour length by summing over individual frame Cl (which is again the sum of inidividual cl of i's)
      for (k=step-tau; k<step; k++){
        cl_total += cl_frame[step-tau];
        if (cl_total< 0.0){
          printf("Cl_total is negative at step %d and tae %d\n", step, tau);
          exit(1);
        }
      }
 

      if (tau == 1 & step == 1){
        printf("Rend: %1.16lf Cl: %1.16lf\n", Rend_total, cl_total);
        if (Rend_total != cl_total){
          printf("Rend_Cl calculation is incorrect!\n");
          exit(1);
        }
      }

      //print Rend (just simply dr2) and contour length 
      //
      //fprintf(msd3, "%1.16lf %1.16lf\n", cl_total, Rend_total+(100.0*(tau-1))); 
      fprintf(msd3, "%1.16lf %1.16lf\n", cl_total, Rend_total); 
      //fprintf(msd3, "%1.16lf %1.16lf\n", cl_total+(1000.0*(tau-1)), Rend_total); 
      //

      if (tau == 2 || tau == 300){
      
        if (step-tau == 0 )
        lagtest = fopen(tt, "w");
        else 
        lagtest = fopen(tt, "a");
      
        fprintf(lagtest, "%1.16lf %1.16lf\n", cl_total, Rend_total); 
        fclose(lagtest);
      
      }



    }//end of step (end of looping over one tau value)

   
    delr2 = MSD_total / double(count);
        
    delx = delx / double(count);
    dely = dely / double(count);
    delz = delz / double(count);
    delx2 = delx2 / double(count);
    dely2 = dely2 / double(count);
    delz2 = delz2 / double(count);
        
    fprintf (msd2, "%d %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf\n", tau, delx, dely, delz, delx2, dely2, delz2, delr2);

    //prints delr^2 per tau value
    fprintf (msd, "%d %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf\n", tau, ((delx2)-(delx*delx)), ((dely2)-(dely)*(dely)), ((delz2)-(delz)*(delz)), ((delx2)-(delx*delx))+((dely2)-(dely)*(dely))+((delz2)-(delz)*(delz)),delr2 );
        
  }//end of tau loop
    
  fclose(msd);
  fclose(msd2);
  fclose(msd3);
    
}



void MSD2() {
  int i, j, k;

  ////////////////////////////////////////
  //  MSD routine with no time average  //
  double MSD_total, dr2, dr[3], delx, dely, delz, delx2, dely2, delz2, delr2;
  int count, tau;
  //contour length - cl, end to end R
  double cl_total, Rend_total;

    
    
  FILE *msd, *msd2, *msd3;
  msd = fopen ( "MSD.dat" , "w" );
  msd2 = fopen ( "MSD_xyz.dat" , "w" );
  msd3 = fopen ( "Cl_Rend.dat" , "w" );
  
  for (i=0; i<N; i++){
    cl[i] = 0.0;
  }
    
  printf("check2\n"); fflush(stdout);
  FILE *otp2;
  for (tau = origin+1; tau < fr_max2; tau++){
    MSD_total = 0.0;
    count = 0;
    cl_total = 0.0;
        
    //MSD loop for each frame of tau intervals
 //   for (step=tau; step < nf ; step++){
      Rend_total = 0.0;
      cl_total = 0.0;

 
      //for all particles
      for (i=0; i<N; i++){

        //if (flag[i] == 0) continue;

        dr2 = 0.0;
        for (j=0; j<3; j++){
                    
          dr[j] = xxx[tau][i][j] - xxx[origin][i][j];
                    
    
          //this is like counting the Rend of particles in adjacent tau frames
          dr2 += dr[j] * dr[j];
        }//end of j (dimension)

        //get contour length of particle i in separate array
//        if (tau == 1){
//          cl[i] += sqrt(dr2); // this is contour length of N bonds. (i (step+1) - i(step) -> N # of i's)
//        }

        //if (step == nf-1 ){
//          Rend[i] = dr2;
//          Rend[i] = sqrt(Rend[i]);
        //}

                
        MSD_total += dr2;
        Rend_total += dr2;
                
        //accumulate delx, dely, .. etc
        delx += dr[0];
        dely += dr[1];
        delz += dr[2];
        delx2 += dr[0] * dr[0];
        dely2 += dr[1] * dr[1];
        delz2 += dr[2] * dr[2];
        count += 1;
                
        }//end of i

      //get contour length of frame
//     if (tau == 1){
        //accumulate contour length for frame = (sum of contour length of i's between two [step] and [step-tau] frames)
//        for (i=0; i<N; i++){
//          cl_frame[step-tau] += cl[i];
//        }
//      }

/*      
      //accumulate total Rend
      for (i=0; i<N; i++){
        if (flag[i] == 0 ) continue;
        Rend_total += Rend[i];
      }
*/
      // accumulate total contour length by summing over individual frame Cl (which is again the sum of inidividual cl of i's)
      for (k=origin; k<tau; k++){
        cl_total += cl_frame[k];
      }

      //Rend_total = sqrt(Rend_total);

      if (tau == 1 & step == 1){
        printf("Rend: %1.16lf Cl: %1.16lf\n", Rend_total, cl_total);
        if (Rend_total != cl_total){
          printf("Rend_Cl calculation is incorrect!\n");
          exit(1);
        }
      }


//    }//end of step (end of looping over one tau value)

   
    delr2 = MSD_total / double(count);

    //print Rend (just simply dr2) and contour length 
    fprintf(msd3, "%1.16lf %1.16lf\n", cl_total, Rend_total); 

   
    delx = delx / double(count);
    dely = dely / double(count);
    delz = delz / double(count);
    delx2 = delx2 / double(count);
    dely2 = dely2 / double(count);
    delz2 = delz2 / double(count);
        
    fprintf (msd2, "%d %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf\n", tau, delx, dely, delz, delx2, dely2, delz2, delr2);

    //prints delr^2 per tau value
    fprintf (msd, "%d %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf\n", tau, ((delx2)-(delx*delx)), ((dely2)-(dely)*(dely)), ((delz2)-(delz)*(delz)), ((delx2)-(delx*delx))+((dely2)-(dely)*(dely))+((delz2)-(delz)*(delz)),delr2 );
        
  }//end of tau loop
    
  fclose(msd);
  fclose(msd2);
  fclose(msd3);
    
}

void write_energy(){
  FILE *otp;
  int i;

  if (step ==0)
    otp = fopen("energy.dat","w");
  else
    otp = fopen("energy.dat","a");

  timesum += delt;

  fprintf(otp, "%lf %1.16lf\n", timesum, PE);
  fclose(otp);

}


void get_contour(){
  int i, j, k;
  int count;
  double dr2, dr[3], r2, rigidity;
  double **olddr, **newdr;
  FILE *otp, *otp2, *otp3;

  olddr = (double**) calloc(N, sizeof (double*));
  newdr = (double**) calloc(N, sizeof (double*));
  for (i=0; i<N; i++){
    olddr[i] = (double*) calloc(3, sizeof (double));
    newdr[i] = (double*) calloc(3, sizeof (double));
  }

  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
        newdr[i][j] = 0.0;
        olddr[i][j] = 0.0;
    }
  }

  otp = fopen("Ediff.dat","w");
  otp2 = fopen("R2_Cl.dat","w");
  otp3 = fopen("rigidity.dat","w");


  for (i=0; i<N; i++){
    cl[i]=0.0;
  }

  for (i=0; i<fr_max2; i++){
    cl_frame[i]=0.0;
  }

  for (step=1; step<fr_max2; step++){

    for (i=0; i<N; i++){
      for (j=0; j<3; j++){
        newdr[i][j] = 0.0;
      }
    }


    for (i=0; i<N; i++){
      cl[i] = 0.0;
    }

    //initialize r2 for msd calc for each frame
    r2 = 0.0;

    //for all aprticles
    for(i=0; i<N; i++){

      if (flag[i]==0) continue;

      dr2 = 0.0;
      for (j=0; j<3; j++){
        dr[j] = xxx[step][i][j] - xxx[step-1][i][j];
        dr2 += dr[j] * dr[j];

        //accumulate 3N displacement vector for frame
        newdr[i][j] = dr[j];

      }//end of j

      
      //get contour length of particle i in separate array
      cl[i] = dr2;


//      printf("cl[i]: %1.16lf\n", cl[i]);
      //get r2 for all particles
      r2 += dr2;
    
    }//end of i

//    printf("clx/y/z: %1.16lf %1.16lf %1.16lf\n", clx,cly,clz);
    //get contour length of frame
    for (i=0; i<N; i++){
      if (flag[i]==0) continue;
      cl_frame[step-1] += cl[i];
    }
    cl_frame[step-1] = sqrt(cl_frame[step-1]);

    printf("cl_frame: %1.16lf\n", cl_frame[step-1]);

    double dumsum;
    dumsum = 0.0;
    //normalize to get unit vector
    for (i=0; i<N; i++){
      for (j=0; j<3; j++){
        newdr[i][j] = newdr[i][j] / cl_frame[step-1];
        dumsum += newdr[i][j] * newdr[i][j];
      }
    }
    printf("sum of unit vector: %1.16lf\n", dumsum);

/*    
    //unit vector of new cl, current cl
    newcl[0] = clx / cl_frame[step-1]; 
    newcl[1] = cly / cl_frame[step-1]; 
    newcl[2] = clz / cl_frame[step-1]; 
*/

    rigidity = 0.0;
    //take dot product of newcl and oldcl
    for (i=0; i<N; i++){
      for (j=0; j<3; j++){
        rigidity += olddr[i][j] * newdr[i][j];
      }
    }

    //the rigidity should always be 1 at very first, b/c old cl is 0 at first frame
    rigidity = 1.0 - rigidity;


    //store into olddr for next frame calculation
    for (i=0; i<N; i++){
      for (j=0; j<3; j++){
        olddr[i][j] = newdr[i][j];
      }
    }


    //get energy difference of frame
    Udiff[step-1] = Usys[step] - Usys[step-1];

    //print contour length vs energy
    fprintf(otp, "%1.16lf %1.16lf\n", cl_frame[step-1], Udiff[step-1]);
    fprintf(otp2, "%1.16lf %1.16lf\n", cl_frame[step-1], r2);
    fprintf(otp3, "%1.16lf\n", rigidity);

  }//end of step

  fclose(otp);
  fclose(otp2);
}
