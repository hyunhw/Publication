/////////////////////////////////////////////////////////
// Once nf max is defined, reads in coarseningtraj.xyz //
// and prints out MSD.dat vs tau                       //
//                                                     //
// Arguments are Frame_min and Frame_Max               //
/////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cmath>

// Define global variables
int N, tau, taustar, step, nf, frmin, frmax;
double ***x, vol, box[3], boxh[3];
char tt[80];


// Define global functions
void read_input (void);
void remove_PBC (void);


///////////////////
// Start of Main //
///////////////////

int main (int argc, char** argv){
  int i, j;
  FILE *inp, *inp2; 
  double MSD_total, dr2; 
  double dr[3], delx, dely, delz, delx2, dely2, delz2, delr, delr2;
  double delxstar, delystar, delzstar;
  int count, countstar;
  char tt[80];

  /*-----Define Simulation Parameters-----*/
  //tau = 5;
  //vol = 2086.295368892912;
  /*--------------------------------------*/

  inp2 = fopen("zwindow.dat","r");
  fscanf(inp2, "%d %d\n", &frmin, &frmax);
  if (argc > 3){
    taustar = atoi (argv[3]);
  }
  nf = frmax - frmin;
  
  // Read in N
  inp = fopen ("coarseningtraj.xyz", "r");
  if (inp == NULL){
    printf("Failed to open input.config\n");
    exit(1);
  }
  fscanf (inp, "%d %lf\n", &N, &vol);
  N = int(N*0.005);
  printf("N: %d vol: %1.16lf\n", N, vol);
  fclose (inp);

  // Calculate box size
  for (i=0; i<3; i++){
    box[i] = pow (vol, 1.0/3.0);
    boxh[i] = box[i] / 2.0;
  }

  //Allocate memory for ***x
  x = (double***) calloc(nf, sizeof (double**));
  for (i=0; i<nf; i++){
    x[i] = (double**) calloc (N, sizeof (double*));
    for (j=0; j<N; j++){
      x[i][j] = (double*) calloc (3, sizeof (double));
    }
  }

  //////////////////
  // Start of MSD //
  //////////////////

  // Read coordinates of all frames
  read_input();
  printf ("number of frames: %d\n", nf); fflush(stdout);

  FILE *otp;
  otp = fopen ("readtest.dat","w");
  for (i=0;i<nf;i++){
    for (j=0; j<N;j++){
      fprintf(otp, "%lf %lf %lf\n", x[i][j][0], x[i][j][1], x[i][j][2]);
    }
    fprintf(otp, "\n");
  }
  fclose(otp);

  remove_PBC();
  
  printf ("read input and removed PBC\n"); fflush(stdout);

  FILE *msd, *msd2;
  msd = fopen ("MSD.dat" , "w");
  msd2 = fopen ("MSD_xyz.dat", "w");

  FILE *otp2;
  // van Hove function
  if (argc > 3){
    sprintf (tt, "vanhove_%d.dat", taustar);
    otp2 = fopen (tt,"w");
  }
  for (tau=1; tau < nf; tau++){
    MSD_total = 0.0;
    count = countstar = 0;
    delx = dely = delz = delx2 = dely2 = delz2 = delr2 = delxstar = delystar = delzstar = 0.0;
    // MSD loop for each frame of tau intervals
    for (step=tau; step < nf; step++){

      // For all particles
      for (i=0; i<N; i++){
        dr2 = 0.0;
        for (j=0; j<3; j++){

          dr[j] = x[step][i][j] - x[step-tau][i][j];

          dr2 += dr[j]*dr[j];
          
        }//end of j 

        MSD_total += dr2;
        //accumulate delx, dely, .. etc
        delx += dr[0]; 
        dely += dr[1]; 
        delz += dr[2]; 
        delx2 += dr[0] * dr[0] ;
        dely2 += dr[1] * dr[1] ;
        delz2 += dr[2] * dr[2] ;
        count += 1;

        // van Hove function
        if (argc > 3 && tau == taustar){ 
          delxstar = dr[0];
          delystar = dr[1];
          delzstar = dr[2];
          fprintf (otp2, "%lf %lf %lf\n", delxstar, delystar, delzstar);
        }

      }//end of i


    }//end of step. end of 1 tau value
    delr2 = MSD_total / double(count);
    
    delx = delx / double(count);
    dely = dely / double(count);
    delz = delz / double(count);
    delx2 = delx2 / double(count);
    dely2 = dely2 / double(count);
    delz2 = delz2 / double(count);

    fprintf (msd2, "%d %lf %lf %lf %lf %lf %lf %lf\n", tau, delx, dely, delz, delx2, dely2, delz2, delr2); 
    fprintf (msd, "%d %lf %lf %lf %lf %lf\n", tau, ((delx2)-(delx*delx)), ((dely2)-(dely)*(dely)), ((delz2)-(delz)*(delz)), ((delx2)-(delx*delx))+((dely2)-(dely)*(dely))+((delz2)-(delz)*(delz)),delr2 ); 

  }//end of tau loop
  
  fclose(msd);
  fclose(msd2);

  if (argc > 3){
    fclose(otp2);
  }


  return 0;

}//end of MAIN

void read_input() {

  int i, frame;
  FILE *inp;
  inp = fopen ("tracertraj.xyz", "r");

  for (frame=0; frame<nf+frmin; frame++){

    //Read the initial frames we don't want
    if (frame < frmin) {
      //read the first 2 lines
      fgets(tt,80,inp);
      fgets(tt,80,inp);

      for (i=0;i<N;i++){
        fgets (tt,80,inp);
      }
    }

    else {
      //red the first two lines
      fgets(tt,80,inp);
      fgets(tt,80,inp);

      // read one frame
      for(i=0;i<N;i++){
        fscanf (inp, "H %lf %lf %lf\n", &x[frame-frmin][i][0], &x[frame-frmin][i][1], &x[frame-frmin][i][2]);
      }
    }
  }

}

void remove_PBC(){
  int i, j, frame;
  double dr[3];

  for (frame=1; frame < nf; frame++){
    for (i=0; i<N; i++){
      for (j=0; j<3; j++){

        /*if (frame == 0){
          dr[j] = 0.0;
        }*/
        //else {
          dr[j] = x[frame][i][j] - x[frame-1][i][j];

          if (dr[j] > boxh[j]){
            x[frame][i][j] -= box[j];
            //dr[j] = x[frame][i][j] - x[frame-1][i][j];
          }
          else if (dr[j] < -boxh[j]){
            x[frame][i][j] += box[j];
            //dr[j] = x[frame][i][j] - x[frame-1][i][j];
          }
        //}
      }
    }
  }
}
