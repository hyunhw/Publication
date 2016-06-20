//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  displacement2.cpp  ----------> This prints (tau vs MSD) (tau vs S)  //
//    |                                                                 //
//    |<------ coarseningtraj.xyz                                       //
//    |<------ zwindow.dat                                              //
//                                                                      //
//    This code post-processes the simulation trajectory to measure     //
//    mean squared displacements as a function of lag time and          //
//    contour length                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <cstdio>
#include <cstdlib>
#include <cmath>

/////////////////////////////
// Define global variables //
int N, nf, start, end, frame_max, *particles;
double ***x, **d, *S, vol, box[3], boxh[3], *diameter;


/////////////////////////////
// Define global functions //
void read_input( void );
void read_file ( void );
void remove_pbc( void );

int main (int argc, char** argv){
  int i, j, k, frame, lastframe, tauc, step;
  double MSD_total, dr2, dr[3], S_total, Ssum;

  frame_max = 8000;

  
  //read and scan for N
  FILE *inp;
  inp = fopen("coarseningtraj.xyz","r");
  fscanf(inp, "%d", &N);
  printf("N: %d\n",N);
  fclose(inp);


  //allocate memory
  x = (double***) calloc(frame_max, sizeof(double**));
  d = (double**) calloc(frame_max, sizeof(double*));
  S = (double*) calloc(frame_max, sizeof(double));
  diameter = (double*) calloc(frame_max, sizeof(double));
  particles = (int*) calloc(frame_max, sizeof(int));
  for (i=0; i<frame_max; i++){
    x[i] = (double**) calloc(N, sizeof(double*));
    d[i] = (double*) calloc(N, sizeof(double));
    for (j=0; j<N; j++){
      x[i][j] = (double*) calloc(3, sizeof(double));
    }
  }

  read_input();
  read_file();

  if (start == 0 && end == 0){
    printf("start and end are both zero, null file\nTerminating simulation\n");
    exit(1);
  }

  for (i=0; i<3; i++){
    box[i] = pow(vol, 1.0/3.0);
    boxh[i] = box[i] / 2.0;
  }

  remove_pbc();
  printf("removed pbc\n"); fflush(stdout);

  //start of MSD calculation for lag coarsening steps

  int loop;
  FILE *otp, *otp2;
  otp = fopen("cMSD.dat","w");
  otp2 = fopen("cS.dat","w");
  int count, count2;
  for (tauc=1; tauc < (end-start+1); tauc++){
    MSD_total = 0.0;
    S_total = 0.0;
    count = count2 = 0;

    for (step = tauc+start; step < end+1; step++){

      //accumulate Ssum
      Ssum = 0.0;
      loop = step-tauc+1;
      while (loop < step +1){ 
        Ssum += S[loop]/diameter[loop];
        loop++;
      }//accumulated all the S micros
      S_total += Ssum;
      count2 += 1;

      //for all particles
      for (i=0; i<N; i++){
        if (d[step][i] == 0) continue;
        dr2 = 0.0;
        for (j=0; j<3; j++){

          dr[j] = x[step][i][j] - x[step-tauc][i][j];

          dr2 += dr[j]*dr[j];
        }//end of j

        MSD_total += dr2;
        count += 1;

      }//end of i

    }//end of step (looping thru frames for one tauc value)
    MSD_total = MSD_total / double(count);
    S_total = S_total / double(count2);

    fprintf(otp, "%d %lf\n", tauc, MSD_total);
    fprintf(otp2, "%d %lf\n", tauc, S_total);
  }//end of tauc

  fclose(otp);
  fclose(otp2);

  return 0;

}

void read_input(){
  int i, j, k;
  FILE *inp;

  inp = fopen("coarseningtraj.xyz","r");

  nf = 0;

  for (i=0; i<frame_max; i++){
    for (j=0; j<N; j++){
      for (k=0; k<3; k++){
        x[i][j][k] = 0.0;
      }
    }
  }

  //while it is not the end of the file
  while (!feof(inp)){
    fscanf(inp, "%d %lf\n\n", &N, &vol);

    for (i=0; i<N; i++){
      //scan for the coordinates
      fscanf(inp, "H %lf %lf %lf %lf\n", &x[nf][i][0], &x[nf][i][1], &x[nf][i][2], &d[nf][i]);
    }

    nf += 1;

  }

  printf("coarsening nf: %d\n", nf);

  fclose(inp);

  inp = fopen("zwindow.dat","r");
  fscanf(inp, "%d %d", &start, &end);
  fclose(inp);

}

void read_file() {
  FILE *inp, *inp2;
  int nf, dummy, dummy22, dummy24;
  double dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7, dummy8, dummy9, dummy10, dummy11, dummy13, dummy14, dummy15, dummy16, dummy17, dummy18, dummy19, dummy20, dummy21, dummy23, dummy25, dummy26, dummy27, dummy28, dummy29, dummy30;


  inp = fopen("gradientstats.dat","r");
  inp2 = fopen("avgdiam.dat","r");
  nf = 0;

  //while it is not the end of file
  while(!feof(inp)){

    //scan for various stats..
    fscanf(inp, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %d %lf %lf %lf %lf %lf %lf %d\n", &dummy, &dummy1, &dummy2, &dummy3, &dummy4, &dummy5, &dummy6, &dummy7, &dummy8, &dummy9, &dummy10, &dummy11, &S[nf], &dummy13, &dummy14, &dummy15, &dummy16, &dummy17, &dummy18, &dummy19, &dummy20, &dummy21, &dummy22, &dummy23, &dummy24, &dummy25, &dummy26, &dummy27, &dummy28, &dummy29, &dummy30, &particles[nf]);
    fscanf(inp2, "%d %lf\n", &dummy, &diameter[nf]);

    nf += 1;

  }
  printf("gradientstats nf : %d\n", nf);
  fclose(inp);
  fclose(inp2);

}

void remove_pbc(){
  int i, j, frame;
  double dr[3];

  for (frame=1; frame<nf; frame++){
    for (i=0; i<N; i++){
      for (j=0; j<3; j++){

        dr[j] = x[frame][i][j] - x[frame-1][i][j];

        if (dr[j] > boxh[j]){
          x[frame][i][j] -= box[j];
        }
        else if (dr[j] < -boxh[j]){
          x[frame][i][j] += box[j];
        }
      }
    }
  }
}
