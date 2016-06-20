////////////////////////////////////
// coarsening simulation          //
// Repulsive Interactions - foam  //
//                                //
// Jennifer Hwang                 //
////////////////////////////////////

#include <cstdio>
#include <time.h>
#include <cstdlib>
#include <cmath>
#define PI 3.141592653589793238462643383

/*------------------------*/
/* Define Global Variable */
/*------------------------*/

int N, N_tracer;                                                                        // Main MD
int print_freq, neighbor_freq, neighbor_freq2;
double T, PE, KE, vol, delt, box[3], boxh[3];
double **x, **v, **f, *d, **s, **elastic, **viscous, **sparticle;
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

double inivol, inivol2, cdum, **temp, **temp2, timesum, **prestep, **nowstep, oPE, clength, **newdr, **olddr, *vals, **coarseninginfo, volfrac;
int mass_freq, temp_frame, nf, *ins, cnf;

/*-------------------------*/
/* Define Global Functions */
/*-------------------------*/

double U_system (void);
double pbc_vdr (double*, double*, double*);
double MSD (double*, double*, double*);
double calc_MSD (void);

void qs(double*, int*, int, int);
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
void steepest_descent( void );

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

void write_energy( void);
void write_xyz (void);
void write_coarsen_xyz(void);
void write_tracer_xyz(void);
void read_input (void);
double ran2 (void);
double dgauss (double, double);

void differentiation(void);

/*------------------------------------------------*/
/*                 Start of Main                  */ 
/*------------------------------------------------*/

int main ( int argc, char** argv) {
  int i, j, k;
  double sum3;
  time_t timer; //type long int

  // Define simulation parameters //
  /*------------------------------*/
  //idum = -long(time(0));
  //store current time into timer
  time(&timer);
  
  // this is called later in the simulation
  // if no additional input was passed as an argument...
  // take timer time as seed
  if (argc < 2 ){
    idum = -timer;
    printf("seed: %ld\n", idum);
  }
  // if there is an additional input is passed
  // take input as seed
  if (argc > 1 ){
    idum = atoi(argv[1]);
    printf("seed: %ld\n", idum);
  }

  do_hessian = 0;

  firestep_max = 1000000;
  end = 100000;
  start = 1000000000;

  N = 4000;
  N_tracer = int(N*0.005);
  n_flag = 0;
//  vol = 2503.5544426715;
  volfrac = 0.75;
  b = 0.0;
  C = 0.05;
  Cf = 0.002;
  ZC = 6.0;

  rcut = 2.5;
  rcut2 = rcut * rcut ;
  rnl = 3.0;
  rnl2 = rnl * rnl ;

  delt = 0.002;
  deltmax = 0.1;

  temp_frame = 1500000;

  //used inside fire
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

  // Initialize center of mass velocity to 0
  double cm[3];
  cm[0] = cm[1] = cm[2] = 0.0;

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

  coarseninginfo = (double**) calloc ( 5000, sizeof (double*));
  for (i=0; i<5000; i++){
    coarseninginfo[i] = (double*) calloc (4, sizeof (double));
  }
  ins = (int *) calloc ( temp_frame, sizeof (int)); 
  vals = (double *) calloc ( temp_frame, sizeof (double)); 
  temp = (double**) calloc( temp_frame, sizeof (double*) );
  temp2 = (double**) calloc( temp_frame, sizeof (double*) );
  for (i=0; i<temp_frame; i++){
    temp[i] = (double*) calloc( 3, sizeof(double) );
    temp2[i] = (double*) calloc( 5, sizeof(double) );
  }
  prestep = (double**) calloc( N, sizeof(double*) );
  nowstep = (double**) calloc( N, sizeof(double*) );
  for (i=0; i<N; i++){
    prestep[i] = (double*) calloc( 3, sizeof(double) );
    nowstep[i] = (double*) calloc( 3, sizeof(double) );
  }
  
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
  sparticle = (double**) calloc (N, sizeof (double*));
  newdr = (double**) calloc (N, sizeof (double*));
  olddr = (double**) calloc (N, sizeof (double*));
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
    sparticle[i] = (double*) calloc (6, sizeof (double));
    newdr[i] = (double*) calloc (3, sizeof (double));
    olddr[i] = (double*) calloc (3, sizeof (double));
    w[i] = (double*) calloc (3, sizeof (double));
    v[i] = (double*) calloc (3, sizeof (double));
    f[i] = (double*) calloc (3, sizeof (double));
    fin[i] = (double*) calloc (3, sizeof (double));
    pre_fin[i] = (double*) calloc (3, sizeof (double));
    rmsd[i] = 0.0;
    flag[i] = 1;

//    for (j=0; j<3; j++){
//      // Randomly place particles in the box with the origin set to the middle of the box
//      x[i][j] = (ran2() * box[j]) - boxh[j];
//      // Assign random velocities to the particles
//      v[i][j] = dgauss (0.0, 1.0);
//      // Accumulate center of mass velocity
//      cm[j] += v[i][j];
//    }
  }


  // Make a gaussian distributed random mixture of particles. Mean diameter = 1 , std dev = 0.5 
  for (i=0; i<(N-N_tracer); i++){
    d[i] = dgauss(1,0.5);
    if (d[i] < 0.0)
      d[i] = fabs(d[i]); 
  }

  // Distribute N_tracer numbers of tracer particles with radius 1
  for (i=(N-N_tracer); i<N; i++){
    d[i] = 2.5;
    flag[i] = 2;
  }

//  print_initial_diameter();

  // Check for total volume in order to use it for mass consv
  inivol = inivol2 = 0.0;
  for (i=0;i<N; i++){
    if (flag[i] == 0) continue;
    if (flag[i] == 2) continue;
    inivol += d[i]*d[i]*d[i]*PI/6.0;
  }
  for (i=0;i<N; i++){
    inivol2 += d[i]*d[i]*d[i]*PI/6.0;
  }
  vol = inivol2/volfrac;
  printf ("Total Particle Volume: %1.16lf, Box volume:%1.16lf\n", inivol2, vol);


  // Set up simulation box
  for (i=0; i<3; i++){
    box[i] = pow ( vol, 1.0/3.0 );
    boxh[i] = box[i] * 0.5;
  }

  //define max cell length (min 3 cells)
  max_cl = box[0]/3.0;

  //assign velocity and position randomly
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      x[i][j] = (ran2() * box[j]) - boxh[j];
      v[i][j] = dgauss(0.0, 1.0);
      cm[j] += v[i][j];
    }
  }

  inivol = inivol2 = 0.0;
  for (i=0;i<N; i++){
    if (flag[i] == 0) continue;
    if (flag[i] == 2) continue;
    inivol += d[i]*d[i]*d[i]*PI/6.0;
  }

  // Normalize CM velocity
  for (j=0; j<3; j++){
    // Divide by N to make CMV an extensive variable
    cm[j] *= (1.0/ double(N));
  }

  // Subtract CM velocity from velocity of each particle so that total force/velocity is 0 on the system
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      v[i][j] -= cm[j];
    }
  }

  // Calculate the KE
  KE = 0.0;
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      KE += v[i][j] * v[i][j];
    }
  }
  KE *= 0.5;

  // Calculate the temperature (3 d.o.f., <KE> = 3/2 NkT)
  T = (2.0 / 3.0) * ( KE / double (N));
  //-----END of initializing system--------------------------------------//
  
  //To calculate PE of the initial system...
  // Allocate memory
  // Assign particles to cell
  // and gather neighborlist
  cell_memory();
  assign_cell();
  neighbor();

  printf("initial PE calculation\n"); fflush(stdout);

  //Calculate the initial PE
  PE = U_system();
  printf ("Initial PE: %1.16lf KE: %1.16lf T: %1.16lf szx: %e, sxz: %e\n", PE, KE, T, s[2][0], s[0][2]);


  //////////////////////
  // Start of MD loop //
  //////////////////////

  //------------------------------------------------COARSENING-----------------------------//


  // ff = 1 so that initial relaxation updates neighbor list DURING relaxation
  fire_flag = 1;
  // Minimize initial random configuration
  fire();
//steepest_descent();
  fire_flag = 0;
  coarsen_step = 0;

 
  // Obtain PE of system and print initial minimized PE 
  PE = U_system();

  // intialize pre_fin for intial calculation
  // pre_fin is the current xyz positions of particles
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      pre_fin[i][j] = x[i][j];
    }
  }


  cnf = 0;

  for (coarsen_step=0; coarsen_step<9999999; coarsen_step++){

    coarsening();
    cnf += 1;

    rnl_update(); //don't need neighbor update because coarsening doesn't change position

   //initial neighbor update is done in first step of fire
    print_max_diam();

    PE = U_system();
    //store coarsening PE
    coarseninginfo[coarsen_step][0] = PE;

    sum3=0.0;
    for (i=0;i<N;i++){
      if (flag[i] == 0) continue;
      if (flag[i] == 2) continue;
      sum3 += d[i]*d[i]*d[i]*PI/6.0;
    }
    printf ("At end of %d'th coarsening.. Total Particle Volume: %1.16lf\n", coarsen_step, sum3);


    // because double float precision is not enough to account for small changes, update mass at large intervals
    if (coarsen_step % mass_freq == 0){
      mass_consv(); 


     sum3=0.0;
      for (i=0;i<N;i++){
        if (flag[i] == 0) continue;
        if (flag[i] == 2) continue;
        sum3 += d[i]*d[i]*d[i]*PI/6.0;
      }
      printf ("After mass consv routine.. Total Particle Volume: %1.16lf\n",  sum3);
    }

    // Relax after Coarsening
    fire();
//steepest_descent();
neighbor_N2();
//store neighbor list cut off 
coarseninginfo[coarsen_step][3] = rnl;

    print_tracer_stress();
    meandistance = calc_MSD();

    differentiation();

    if (coarsen_step == 81){
      FILE *test1, *test2;
      test1 = fopen("temp.dat","w");   
      test2 = fopen("temp2.dat","w");   
      for (i=0; i<nf; i++){
        fprintf(test1, "%1.16e %1.16e %1.16e\n", temp[i][0], temp[i][1], temp[i][2]);
        fprintf(test2, "%1.16e %1.16e %1.16e %1.16e %1.16e\n", temp2[i][0], temp2[i][1], temp2[i][2], temp2[i][3], temp2[i][4]);
      }
      fclose(test1); 
      fclose(test2);
    }


    //print_Pr();
  
    print_iso();
 
    //print_tracer_coordination();

    //store FIRE PE
    coarseninginfo[coarsen_step][1] = PE;
    coarseninginfo[coarsen_step][2] = meandistance;

    // print xyz at frequency
    if (coarsen_step % coarsen_print_freq == 0){
      write_coarsen_xyz();
      //write_tracer_xyz();
    }
    
  }//END of coarsening loop


  //do_hessian = 1;
  //hessian_memory();
  //hessian_ana();
  //print_hessian();

  //print volume to check for mass conservation
  sum3 = 0.0;
  for (i=0;i<N;i++){
    if (flag[i] == 0 || flag[i] == 2 )
      continue;
    sum3 += d[i]*d[i]*d[i]*PI/6.0;
  }
  printf ("Total Particle Volume: %lf\n", sum3);
  printf ("Total Particle (%d - %d - %d): %d\n", N, n_flag, N_tracer, N-n_flag-N_tracer);

  return 0;
}

//Calculate dr_ij with PBC
//Don't read input file in pbc-vdr for efficiency because this is in the U_system loop
double pbc_vdr ( double x1[3], double x2[3] , double dr[3] ){
  double mdr2 = 0.0;
  int i;

  //PBC in all x,y,z directions
  if (step<start){
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
  }
      
  //PBC only in x and y directions (so that walls don't feel the force)
  else {

      for (i=0; i<2; i++){
        dr[i] = x1[i] - x2[i];

        if (dr[i] > boxh[i]){
          dr[i] -= box[i];
        }
        else if ( dr[i] < -boxh[i] ){
          dr[i] += box[i];
        }
        mdr2 += dr[i] * dr[i];
      }
      dr[2] = x1[2] - x2[2];
      mdr2 += dr[2] * dr[2];
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

  //calculate stress for the individual bubbles
  for (i=0; i<N; i++){
    for (j=0; j<6; j++){
      sparticle[i][j] = 0.0;
    }
  }

  for (i=0; i<N; i++){
    if (flag[i] == 0) continue;
    for (k=0; k<n_neighbor[i]; k++){
      j = n[i][k];
      if (j == i) continue;
      if (flag[j] == 0) continue;

      dr2 = pbc_vdr (x[i], x[j], dr);
      dij = (d[i] + d[j])/2.0;
      dij2 = dij * dij;
      if (dr2 > dij2) continue;
      rij = sqrt(dr2);
      F_pre = (2.0 * (1.0/(rij*dij) - 1.0/dij2));
      
      for (z=0; z<3; z++){
        for (zz=0; zz<3; zz++){
          if (z==zz){
            sparticle[i][z] += -(N*T/vol) - (8.0/(d[i]*d[i]*d[i]))*(F_pre*dr[z])*dr[zz];
          }
          else if (zz>z){
            sparticle[i][z+zz+2] += -(8.0/(d[i]*d[i]*d[i]))*(F_pre*dr[z])*dr[zz];
          }
        }//end of zz
      }//end of z
    }//end of k
  }//end of i

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

//Read coordinates of N particles from input.config file and store it into **x
void read_input () {
  int i;
  FILE *inp;
  inp = fopen ( "input.config", "r");

  if ( inp == NULL ){
    printf ("Failed to open input.config\n");
    exit(1);
  }

  fscanf ( inp, "%d\n", &N);

  //Allocate memory for **x
  x = (double**) calloc (N, sizeof (double*));
  //Allocate and scan input file and store in **x
  for (i=0; i<N; i++){
    x[i] = (double*) calloc (3, sizeof (double));

    fscanf (inp, "%lf %lf %lf\n", &x[i][0], &x[i][1], &x[i][2]);
  }

  fclose (inp);
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
  FILE *coarsen, *velocity;
  int i;
  
  if (coarsen_step ==0) {
    coarsen = fopen ("coarseningtraj.xyz", "w");
    velocity = fopen ("velocity.dat", "w");
  }
  else {
    coarsen = fopen ("coarseningtraj.xyz", "a");
    velocity = fopen ("velocity.dat", "a");
  }

  // Writes number of particles and an extra blank line
  fprintf (coarsen, "%d %1.16lf\n\n", N, vol);
  fprintf (velocity, "%d\n\n", N);

  for (i=0; i<N; i++){
    fprintf (coarsen, "H %1.16lf %1.16lf %1.16lf %1.16lf\n", x[i][0], x[i][1], x[i][2], d[i]);
    fprintf (velocity, "H %1.16lf %1.16lf %1.16lf %1.16lf\n", v[i][0], v[i][1], v[i][2], d[i]);
  }

  fclose (coarsen);
  fclose (velocity);
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
  FILE *otp, *otp2;
  int i;

  if (coarsen_step == 0){
    otp = fopen ("tracerstress.dat", "w");
    otp2 = fopen ("particlestress.dat", "w");
  }
  else { 
    otp = fopen ("tracerstress.dat", "a");
    otp2 = fopen ("particlestress.dat", "a");
  }

  // Writes number of tracer particles and an extra blank line
  fprintf (otp, "%d\n\n", N_tracer);
  fprintf (otp2, "%d\n\n", N);

  for (i=0; i<N_tracer; i++){
    fprintf (otp, "%d %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf\n", i+(N-N_tracer), stracer[i][0][0], stracer[i][1][1], stracer[i][2][2], stracer[i][0][1], stracer[i][0][2], stracer[i][1][2]);
  }

  for (i=0; i<N; i++){
    fprintf (otp2, "%1.16lf %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf\n", sparticle[i][0], sparticle[i][1], sparticle[i][2], sparticle[i][3], sparticle[i][4], sparticle[i][5]);
  }

  fclose(otp);
  fclose(otp2);
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
    if (i<(N/2))
      fprintf (otp, "He %lf %lf %lf\n", x[i][0], x[i][1], x[i][2] );
    else
      fprintf (otp, "H %lf %lf %lf\n", x[i][0], x[i][1], x[i][2] );
  }
  fclose (otp);
}

void fire(){
  int i, j, z;
  int Ncount, doprint;
  double vel, force, P, Uo, conv, alpha, min, max, dr2, dr[3];
  FILE *otp, *stp, *coordi;
  double check = 0.0;
  double clengthnow;
  
  
  for (i=0;i<N;i++){
    rmsd[i] = 0.0;
  }

  if (coarsen_step == 0 && fire_flag == 0){
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
  clength = 0.0;
  timesum = 0.0;
  nf = 0;
  doprint = 0;

  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      prestep[i][j] = x[i][j];
    }
  }

  for (i=0; i<temp_frame; i++){
    for(j=0; j<3; j++){
      temp[i][j] = 0.0;
    }
  }

  oPE = PE;
  
  for (step=0; step<firestep_max; step++){

    if (fire_flag==0) doprint = 1;

    if (nf > (temp_frame-1)){
      printf("nf exceeds number of frames allocated for temp array!\n");
      exit(1);
    }

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
    min = 50000.0;

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

    //accumulate fire time
    timesum += delt;

    //store updated coordinates into nowstep
    for (i=0; i<N; i++){
      for (j=0; j<3; j++){
        nowstep[i][j] = x[i][j];
      }
    }

    clengthnow = 0.0;

    //calculates contour length no matter what doprint value is
    //calculate contour length between previous and current step
    for (i=0; i<N; i++){
      if (flag[i] == 0) continue;
      dr2 = 0.0;
      for (j=0; j<3; j++){
        dr[j] = nowstep[i][j] - prestep[i][j];

        if (dr[j] > boxh[j]) dr[j] -= box[j];
        else if (dr[j] < -boxh[j]) dr[j] += box[j];

        dr2 += dr[j] * dr[j]; //r^2 of particle i
      }//end of j
      clengthnow += dr2;
    }//end of i

    clengthnow = sqrt(clengthnow);

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
    double maxx, force2;
    maxx = -5.0;
    for (i=0; i<N; i++){
      if (flag[i] == 0 ) continue;
      force2 = 0.0;
      for (j=0; j<3; j++){
        force2 += f[i][j] * f[i][j];
      }
      if (force2 > maxx) maxx = force2;
    }
    maxx = sqrt(maxx);
    if (maxx < 0.0000001) break;
*/    

    //Convergence check
    conv = fabs(PE - Uo);

    if ( conv < 0.00000000000001 ){
      break; //break out of step loop
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

    //only when velocity was not resetted, update clength and store these in temp array
    if (doprint == 1){
      for (i=0; i<N; i++){
        for (j=0; j<3; j++){
          prestep[i][j] = nowstep[i][j];
        }
      }
      clength += clengthnow;

      //store accumulated contour length into temp array
      temp[nf][2] = clength;

      //store energy
      temp[nf][1] = PE;

      //store accumulated fire time into temp array
      temp[nf][0] = timesum;

      nf += 1;
    }

    /*--------------------------------------------------------------------------------------*/

    /*
    if (step % print_freq == 0 && doprint == 1){
      //printf ("Step: %d, Total E: %1.12lf, PE: %1.12lf, KE: %lf, T: %lf, szx: %e, sxz: %e\n", step, KE+PE, PE, KE, T, s[2][0], s[0][2] );
      nf += 1;
    }*/

    //add to nf to count number of fire frames
  
  }//End of step loop
  printf ("FIRE: step:%d nf: %d  PE:%1.16lf\n", step, nf, PE);
  if (coarsen_step >= 0){
    fprintf (otp, "Coarsen step: %d FIRE step: %d PE:%1.16lf N: %d\n", coarsen_step, nf, PE, (N-n_flag));
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
    FILE *nrl, *coarsen;
    nrl = fopen("rnl.dat","w");
    coarsen = fopen("coarsening.dat","w");
    for (i=0; i<cnf-1; i++){
      fprintf(coarsen, "%d %1.16e %1.16e %1.16e\n", i, coarseninginfo[i][0], coarseninginfo[i][1], coarseninginfo[i][2]);
      fprintf(nrl, "%1.16lf\n", coarseninginfo[i][3]);
    }
    fclose(coarsen);
    fclose(nrl);

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
  double MSD_total, dr[3], top, bot, rigid, length;
  FILE *rigidity;

  if (coarsen_step == 0) rigidity = fopen("rigidity.dat","w");
  if (coarsen_step > 0) rigidity = fopen("rigidity.dat","a");

  Pr = MSD_total = top = bot = 0.0;
  for(i=0;i<N;i++){
      msd[i] = 0.0;
  }

  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      newdr[i][j] = 0.0;
    }
  }

  //before calculating, store current xyz as final coordinates
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      fin[i][j] = x[i][j];
    }
  }

  count = 0;
  // accumulate r2 for every particle
  for (i=0; i<N; i++){
    if (flag[i] == 0) continue;
    dr2 = 0.0;
    for (j=0; j<3; j++){
      dr[j] = fin[i][j] - pre_fin[i][j];

      if (dr[j] > boxh[j]){
        dr[j] -= box[j];
      }
      else if (dr[j] < -boxh[j]) {
        dr[j] += box[j];
      }
      //accumulate 3N displacement vector for this frame (calc_MSD is called after fire)
      newdr[i][j] = dr[j];
      dr2 += dr[j]*dr[j]; 
    }

    //dr2 = sqrt(dr2);
    msd[i] = dr2;
    MSD_total += dr2;
    count += 1;
  }//end of i

  //get length_r for unit vector cal
  length = sqrt(MSD_total);
  //normalize MSD by number of particles
//  MSD_total = MSD_total / double(count);

  double dumsum = 0.0;
  //get unit vector for this frame
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      newdr[i][j] = newdr[i][j] / length;
      dumsum += newdr[i][j] * newdr[i][j];
    }
  }
  
  //check to see if unit vector sum is 1
  if (dumsum > 1.01 || dumsum < .99 ){
    printf("unit vector calculation is wrong! exiting simulation\n");
    printf("sum: %1.16e\n", dumsum);
    exit(1);
  }

  rigid = 0.0;

  //take dot product of unit vectors (newcl and oldcl)
  if (coarsen_step > 0){
    for (i=0; i<N; i++){
      for (j=0; j<3; j++){
        rigid += olddr[i][j] * newdr[i][j];
      }
    }
    rigid = 1.0 - rigid;
  }

  fprintf(rigidity, "%d %1.16e\n", coarsen_step, rigid);
  fclose(rigidity);

  //store current unit vector to old unit vector
  for (i=0; i<N; i++){
    for (j=0; j<3; j++){
      olddr[i][j] = newdr[i][j];
    }
  }

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

    FILE *nrl, *coarsen;
    nrl = fopen("rnl.dat","w");
    coarsen = fopen("coarsening.dat","w");
    for (i=0; i<cnf-1; i++){
      fprintf(coarsen, "%d %1.16e %1.16e %1.16e\n", i, coarseninginfo[i][0], coarseninginfo[i][1], coarseninginfo[i][2]);
      fprintf(nrl, "%1.16lf\n", coarseninginfo[i][3]);
    }
    fclose(coarsen);
    fclose(nrl);

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

void steepest_descent(){
  int i, j, k;
  double force, lambda, oldPE, min, max, check;

  lambda = 0.002;
  force = 0.0;

  //get energy of system
  PE = U_system();
  oldPE = PE;

  //step is global variable
  for (step=0; step<start; step++){
/*    if (step % neighbor_freq == 0 ){
      assign_cell();
      neighbor();
    }
*/    

    //update neighbor at regular frequency 
    if (max > (0.5*rnl) || step % neighbor_freq2 == 0){
      if (rnl > max_cl){
        neighbor_N2();
        check = U_system();
        if (fabs(check-oldPE)/double(N-n_flag) > 0.00001 && step != 0){
          printf ("Change in PE detected, neighbor list is faulty\n");
          printf ("check: %1.16lf, Uo: %1.16lf\n", check, oldPE);
          exit(1);
        }
      }
      else {
        assign_cell();
        neighbor();
        check = U_system();
        if (fabs(check-oldPE)/double(N-n_flag) >0.00001 && step != 0){
          printf ("Change in PE detected2, neighbor list is faulty\n");
          printf ("check: %1.16lf, Uo: %1.16lf\n", check, oldPE);
          exit(1);
        }
      }
    }


    oldPE = PE;

    min = 5000000.0;
    max = -5.0;

    for (i=0; i<N; i++){

      if (flag[i] == 0) continue;
      rmsd[i] = 0.0;

      for (j=0; j<3; j++){
        if (flag[j] == 0 ) continue;
        x[i][j] += lambda * f[i][j];
        
        //apply PBC
        if (x[i][j] > boxh[j])
          x[i][j] -= box[j];
        else if (x[i][j] < -boxh[j])
          x[i][j] += box[j];
      }

      if (rmsd[i] < min)
        min = rmsd[i];
      else if (rmsd[i] > max)
        max = rmsd[i];

    }

//write_energy();
//write_xyz();

    PE = U_system();
    if (fabs(PE-oldPE) < 0.0000000001) break;

  }//end of step loop

}

void write_energy(){
  int i;
  FILE *otp;

  if (step ==0) otp = fopen("energy.dat","w");
  else otp = fopen("energy.dat","a");
  
  PE = U_system();

  fprintf(otp, "%d %1.16lf\n", step, PE);

  fclose(otp);
}

void differentiation( ){
  int i, j, k, count, lowest_third_index, highest_third_index;
  FILE *otp, *tpm;
  double max, min, maxt, mint, maxs, mins, avg, std, sqavg, sqstd, lowest_third_value, highest_third_value;

  max = -5000.0;
  min = 5000.0;
  lowest_third_index = highest_third_index = 0;
  lowest_third_value = highest_third_value = 0.0;

  if (coarsen_step == 0){ 
    otp = fopen("gradientstats.dat","w");
    tpm = fopen("tpm.dat","w");
  }
  else if (coarsen_step != 0){
    otp = fopen("gradientstats.dat","a");
    tpm = fopen("tpm.dat","a");
  }

  for (i=0; i<temp_frame; i++){
    ins[i] = 0;
    vals[i] = 0.0;
    for(j=0; j<5; j++){
      temp2[i][j] = 0.0;
    }
  }

  //////////////
  // gradient //
  //////////////


  //backwards difference in energy
  for (i=1; i<nf; i++){
    // dE/dS_backwards
    temp2[i][0] = (temp[i][1] - temp[i-1][1]) / (temp[i][2] - temp[i-1][2]);
  }//index: [1,nf]

  //forwards difference in energy
  for (i=0; i<nf-1; i++){
    // dE/dS_forwards
    temp2[i][1] = (temp[i+1][1] - temp[i][1]) / (temp[i+1][2] - temp[i][2]);
  }//index: [0,nf-1]

  //gradient obtained by averaging backwards and forwards difference
  for (i=1; i<nf-1; i++){
    temp2[i-1][2] = (temp2[i][0] + temp2[i][1]) / 2.0;
    temp2[i-1][3] = temp[i][2];  // corresponding contour length for gradient index
    temp2[i-1][4] = temp[i][0];  // corresponding time for gradient index
  }//index: [0,nf-2]
    

  int countt=0;
  int minindex, maxindex;

minindex = -1;
maxindex = -1;

  //get index for min and mx S (range that we want to calc max/min grad)
  for (i=0; i<nf-2; i++){
    if (temp2[i][3] > (0.25+temp2[0][3])){
      countt+=1;
      if (countt == 1){
        minindex = i;
      }
    }
    if ((0.25+temp2[i][3]) < temp2[nf-3][3]){
      maxindex = i;
    }
  }

if (maxindex > nf && maxindex > 0){
printf("Max index is greater than nf!\nTerminating simulation\n");
exit(1);
}
if (minindex > nf && minindex > 0){
printf("Min index is greater than nf!\nTerminating simulation\n");
exit(1);
}

/*
  FILE *otp2;
  FILE *otp3;

  if (coarsen_step == 0) otp3 = fopen("index.dat","w");
  if (coarsen_step != 0) otp3 = fopen("index.dat","a");
  fprintf(otp3, "%d %d %d\n", coarsen_step, minindex, maxindex);
  fclose(otp3);

  
  if (coarsen_step == 150){
    otp2 = fopen("firegradient.dat","w");
    for (i=0; i<nf-2; i++){
      fprintf(otp2, "%1.16e %1.16e\n", temp2[i][2], temp2[i][3]);
    }
    //fprintf(otp2, "%d %d %1.16e\n", minindex, maxindex);
    fclose(otp2);
  }
*/

  fprintf(tpm, "%d ", -coarsen_step); 
  double modd=0.0;
  int lastindex = minindex;
  int printcount=0;

  if ( (minindex != -1) && (maxindex != -1)){
    if( (temp2[maxindex][3] - temp2[minindex][3]) > 1.0  && maxindex > minindex ) {

      avg = std = 0.0;
      sqavg = sqstd = 0.0;
      count = 0;

      for (i=minindex; i<maxindex; i++){
        //accumulate sum to average
        avg += temp2[i][2];
        sqavg += sqrt(fabs(temp2[i][2]));
        
        count += 1;
    
        modd = (temp2[i][3] - temp2[lastindex][3]);//modd is 0.0 very first iteration

        if (modd == 0.0 || modd > 0.5){
          if (printcount == 1){
            fprintf(tpm, "%1.16lf %1.16lf ", temp2[minindex][2], temp2[i][2]);

            lastindex = i;
            printcount += 1;
          }
          else if (printcount > 1){
            fprintf(tpm, "%1.16lf %1.16lf ", temp2[lastindex][2], temp2[i][2]);
         
            lastindex = i;
            printcount += 1;
          }
          else {
            printcount += 1;
          }
        }//end of if modd
     

        if (temp2[i][2] > max){
          max = temp2[i][2];
          maxs = temp2[i][3];
          maxt = temp2[i][4];
        }
        else if (temp2[i][2] < min){
          min = temp2[i][2];
          mins = temp2[i][3];
          mint = temp2[i][4];
        }
      }//done with index loop
 
      //get the average gradient
      avg = avg / double(count);
      sqavg = sqavg / double(count);
      
      count = 0;
      //loop again to accumulate standard deviation
      for (i=minindex; i<maxindex; i++){
        std += (temp2[i][2] - avg)*(temp2[i][2] - avg);
        sqstd += (temp2[i][2] - sqavg)*(temp2[i][2] - sqavg);
        count +=1;
      }
      std = std / double(count);
      sqstd = sqstd / double(count);
      std = sqrt(std);
      sqstd = sqrt(sqstd);

      int dummmm = 0;
      //copy index and gradient values into ins and vals to sort
      for (i=minindex; i<maxindex; i++){
        ins[dummmm] = i;
        vals[dummmm] = temp2[i][2];
        dummmm += 1;
      }

      highest_third_index = int ( 0.97 * double(dummmm));
      lowest_third_index  = int ( 0.03 * double(dummmm));

      qs (vals, ins, 0, dummmm-1);

      highest_third_value = vals[highest_third_index];
      lowest_third_value  = vals[lowest_third_index];

      double dummp = 0.0;
      dummp = -temp2[0][2];
      if (dummp < 0.0){
        printf("start gradient is positive!\n");
      }

      fprintf(otp, "%d %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %d %1.16e %d %1.16e %1.16e %1.16e %lf %1.16e %1.16e %d\n", coarsen_step, max, maxs, maxt, min, mins, mint, temp2[0][2], temp2[nf-3][2], temp2[nf-3][4], fabs(PE-oPE), meandistance, clength , avg, std, sqrt(fabs(temp2[minindex][2])), temp2[minindex][3], sqrt(fabs(temp2[maxindex][2])), temp2[maxindex][3], sqrt(fabs(temp2[0][2])), sqavg, sqstd, lowest_third_index, lowest_third_value, highest_third_index, highest_third_value, temp2[lowest_third_index][3], temp2[highest_third_index][3], Ziso, (temp[maxindex][1] - temp[minindex][1]), (temp[maxindex][2]-temp[minindex][2]), (N-n_flag) );
    }//end of if maxindex > minindex

    else{ fprintf(otp, "%d %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %d %1.16e %d %1.16e %1.16e %1.16e %lf %1.16e %1.16e %d\n", coarsen_step, 0.0 , 0.0 , 0.0 , 0.0, 0.0, 0.0, temp2[0][2], 0.0, 0.0, fabs(PE-oPE), meandistance, clength , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0, 0.0, 0.0, 0.0, Ziso, 0.0, 0.0, (N-n_flag));
  }

  }//end of if maxindex !=-1 and minindex != -1 

  else{ fprintf(otp, "%d %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %d %1.16e %d %1.16e %1.16e %1.16e %lf %1.16e %1.16e %d\n", coarsen_step, 0.0 , 0.0 , 0.0 , 0.0, 0.0, 0.0, temp2[0][2], 0.0, 0.0, fabs(PE-oPE), meandistance, clength , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0, 0.0, 0.0, 0.0, Ziso, 0.0, 0.0, (N-n_flag));
  }

  fclose(otp);
  fclose(tpm);
}

void qs(double *dr, int *inds, int left, int right){

  int i, j, yi, xi;
  double x, y;
  i = left; j=right;
  x = dr[(left+right)/2];

  do {
    while((dr[i] < x) && (i < right)) i++;
    while((x < dr[j]) && (j > left )) j--;

    if (i <= j) {
      y = dr[i];
      dr[i] = dr[j];
      dr[j] = y;

      yi = inds[i];
      inds[i] = inds[j];
      inds[j] = yi;

      i++; j--;
    }
  } while (i <= j);

  if (left < j)  qs(dr, inds, left, j);
  if (i < right) qs(dr, inds, i, right);

}
