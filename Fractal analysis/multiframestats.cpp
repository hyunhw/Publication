//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  multiframestats.cpp  ------------> S_MSD_%d.dat , avgdiam.dat           //
//     ^                                                                    //
//     |---- coarseningtraj.xyz :     input(  interval (n)  )               //
//     |---- zwindow.dat                                                    //
//     |---- gradientstats.dat                                              //
//                                                                          //
//         start / end = Z-Zc window                                        //
//         0 = write, 1 = append                                            //
//                                                                          //
//    *this file reads in coarseningtraj.xyz and calculates                 //
//     the S and MSD for every n intervals between [start,end]              // 
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


#include <cstdio>
#include <cstdlib>
#include <cmath>


/////////////////////////////
// define global variables //
int N, nf, start, end, nn, frame_max, *particles;
double ***x, **d, **info, *S, vol, box[3], boxh[3];



/////////////////////////////
// define global functions //
void read_input ( void );
void read_file ( void );
void remove_pbc ( void );

int main(int argc, char** argv){
  int i, j, k, frame, lastframe, infocount;
  double clengthnow, clength, SD, dr2, dr[3], Ssum, *avgd;


  frame_max = 8000;


  //read and store N
  FILE *inp;
  inp = fopen("coarseningtraj.xyz","r");
  fscanf(inp, "%d", &N);
  printf("N: %d\n",N);
  fclose(inp);

  if (argc < 2){
    nn=10;
    printf("frame number not specified! setting nn to %d\n", nn);
  }
  else {
    nn = atoi(argv[1]);
  }

 // start = atoi(argv[1]);
 // end = atoi(argv[2]);
 // vol = atof(argv[4]);


  //allocate memory for x
  x = ( double*** ) calloc ( frame_max, sizeof (double**) );
  d = ( double** ) calloc ( frame_max, sizeof (double*) );
  S = ( double* ) calloc ( frame_max, sizeof (double) );
  particles = ( int* ) calloc ( frame_max, sizeof (int) );
  avgd = ( double* ) calloc ( frame_max, sizeof (double) );
  info = ( double** ) calloc ( frame_max, sizeof (double*) );
  for (i=0; i<frame_max; i++){
    x[i] = ( double** ) calloc ( N, sizeof (double*)) ;
    d[i] = ( double* ) calloc ( N, sizeof (double)) ;
    info[i] = ( double* ) calloc ( 2, sizeof (double)) ;
    for (j=0; j<N; j++){
      x[i][j] = ( double* ) calloc ( 3, sizeof (double) );
    }
  }


  

  read_input();

  if (start == 0 && end == 0){
    printf("Start and end are both zero, null file\nTerminating simulation\n"); 
    exit(1);
  }

  if (nn >= (end-start)){ 
    printf("nn is greater than window!\nTerminating simulation\n");
    exit(1);
  }

  for (i=0; i<3; i++){
    box[i] = pow(vol, 1.0/3.0);
    boxh[i] = box[i] / 2.0;
  }
  printf("box: %lf %lf %lf\n", box[0], box[1], box[2]);
  printf("vol: %lf, start: %d, end: %d\n", vol, start, end); fflush(stdout);
  remove_pbc();
  printf("removed pbc\n"); fflush(stdout);
  read_file();
  printf("read gradient stats\n"); fflush(stdout);


  double dummysum;
  int dummyflag=0;

  FILE *otpd;
  otpd = fopen("diameterforfire.dat","w");
  //
  //
  //calculate average diameter per frame
  for (i=0; i<nf; i++){
    for (j=0; j<N; j++){
      if (i == 1249){
        fprintf(otpd, "%1.16lf\n", d[i][j]);
      } 
      //accumulate sum of diameter for N particles
      dummysum += d[i][j];
      if (d[i][j] != 0.0) dummyflag += 1;
    }//end of j (particles)
    avgd[i] = dummysum / double(dummyflag);
    if (avgd[i] != dummysum / double(particles[i])){
      printf("particle calculation, and thus average diameter is wrong!\n");
      exit(1);
    } 
    dummysum = 0.0;
    dummyflag = 0;
  }//end of calculating avg diameter
  fclose(otpd);
   
  FILE *otp2;
  otp2 = fopen("avgdiam.dat","w");
  for (i=0; i<nf; i++){
    fprintf(otp2, "%d %1.16lf\n", i, avgd[i]);
  }

  fclose(otp2);
 
  exit(1);  


  int flag=0;
  //calculate Macro S(S_coarsening) and MSD between nn frames
  infocount = 0;
  //last frame is start - nn because we are taking MSD.
  //lastframe = start-nn;
  lastframe = start;
  clength = SD = Ssum = 0.0;
  for (frame = start; frame < end+1; frame++){
    clengthnow = 0.0; // individual S macro value

    if ( frame != start+1 && (frame-start) % nn == 1) {   // if we're at every nn frames..
    //remainder is 1 because we are going to add S[frame] at nn+1 frame b/c S is calculated for 20->21 (and is printed for "21")
      clength = SD = 0.0;              // initialize S macro and SD (for n frames) to zero
      lastframe = frame-1;            // updating last frame so that we can print at every -1 from nth frames
      printf("frame: %d initializing clength, SD, and lastframe: %d\n", frame, lastframe);
      Ssum = 0.0;                      // initialize S micro to zero
    }


    // basically adding stats for every frame except for the very first frame (==start)
    if (frame != start) {

      //Ssum += S[frame]/double(particles[frame]);               //sum of S micro for n frames
      Ssum += S[frame]/(avgd[frame]/2.0);            //sum of S micro (normalized by <r>) for n frames
      printf("frame= %d S[frame]= %lf\n", frame, S[frame]);

      if ((frame-lastframe) % nn == 0){   //every 1 frame before n frames
                                               //Calculate MSD and store MSD, S micro, and S macro in info array
        printf("MSD calc: frame: %d, lastframe: %d\n", frame, lastframe);
        for (i=0; i<N; i++){
          dr2 = 0.0;
//          if (d[frame][i] == 0.0) continue;

          for (j=0; j<3; j++){
            dr[j] = x[frame][i][j] - x[lastframe][i][j];
           

            dr2 += dr[j] * dr[j]; 
//            if (frame==start) printf("i: %d j: %d dr[j]: %1.16lf\n", i, j, dr[j]);
          }//end of j

//            if (frame==start) printf("i: %d d[i]: %1.16lf dr2: %1.16lf\n",i, d[frame][i], dr2);
          SD += dr2;
        }//end of i

        //info[infocount][1] = SD/double(particles[frame]*particles[frame]);
        info[infocount][1] = SD/(avgd[lastframe]*avgd[lastframe]/4.0);
        info[infocount][0] = Ssum;
      printf("printing frame: %d Ssum: %lf\n", frame, Ssum);
        infocount += 1;
      }//end of if(frame-lastframe)

    }//end of if frame != start
  }//end of frame loop


  printf("infocount: %d\n", infocount);

  //print info array to final file
  FILE *otp;
  otp = fopen ("S_MSD.dat","w");
  for (i=0; i<infocount; i++){
    fprintf(otp, "%1.16e %1.16e\n", info[i][0], info[i][1]);
  }
  fclose(otp);
    



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
    //scan for the number of particles
    fscanf(inp, "%d %lf\n\n", &N, &vol);

    //then scan or the coordinates of N atoms in that frame
    for (i=0; i<N; i++){
      fscanf(inp, "H %lf %lf %lf %lf\n", &x[nf][i][0], &x[nf][i][1], &x[nf][i][2], &d[nf][i]);
    }

    nf += 1;
  }

  printf("coarsening nf : %d\n", nf);

  fclose(inp);

  inp = fopen("zwindow.dat","r");
  fscanf(inp, "%d %d", &start, &end);
  fclose(inp);


}


void read_file(){
  int dummy, dummy22, dummy24;
  double dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7, dummy8, dummy9, dummy10, dummy11, dummy13, dummy14, dummy15, dummy16, dummy17, dummy18, dummy19, dummy20, dummy21, dummy23, dummy25, dummy26, dummy27, dummy28, dummy29, dummy30;

  FILE *inp;

  inp = fopen("gradientstats.dat","r");

  nf = 0;

  printf("entered read file\n"); fflush(stdout);
  //while it is not the end of file
  while (!feof(inp)){

    //scan for various stats.. the variable of interest (S micro) is the 13th variable (dummy12)
    fscanf(inp, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %d %lf %lf %lf %lf %lf %lf %d\n", &dummy, &dummy1, &dummy2, &dummy3, &dummy4, &dummy5, &dummy6, &dummy7, &dummy8, &dummy9, &dummy10, &dummy11, &S[nf], &dummy13, &dummy14, &dummy15, &dummy16, &dummy17, &dummy18, &dummy19, &dummy20, &dummy21, &dummy22, &dummy23, &dummy24, &dummy25, &dummy26, &dummy27, &dummy28, &dummy29, &dummy30, &particles[nf]); 

//    printf("%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %d %lf %lf %lf %lf %lf %lf %d\n", dummy, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7, dummy8, dummy9, dummy10, dummy11, S[nf], dummy13, dummy14, dummy15, dummy16, dummy17, dummy18, dummy19, dummy20, dummy21, dummy22, dummy23, dummy24, dummy25, dummy26, dummy27, dummy28, dummy29, dummy30, dummy31); 
    nf += 1;
//  printf("nf: %d\n", nf); fflush(stdout);
  }

  printf("gradientstats nf : %d\n", nf);

  fclose(inp);

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



