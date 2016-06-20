//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  hurstexponent.cpp  ------------>  2dprojection.dat                      //
//     ^                              hexpomoments.dat                      //
//     |                              pdistribution.dat                     //
//     |                                                                    //
//     |---- coarseningtraj.xyz :     input(  interval of analysis  )       //
//     |---- zwindow.dat                                                    //
//                                                                          //
//                                                                          //
//         start / end = Z-Zc window                                        //
//         0 = write, 1 = append                                            //
//                                                                          //
//    *this file reads in coarseningtraj.xyz and calculates                 //
//     the hurst exponent for the input frame                               // 
//    *transverse excursion analysis                                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


#include <cstdio>
#include <cstdlib>
#include <cmath>


/////////////////////////////
// define global variables //
int N, nf, start, end, nn, frame_max, *particles, rem, printcounter, get2dprojection;
double ***x, **d, **info, *S, vol, box[3], boxh[3];

double **unit, SD, ***P, ***Q, *Pframe, meanP, variance, thirdP, fourthP, **Pmax;



/////////////////////////////
// define global functions //
void read_input ( void );
void read_file ( void );
void remove_pbc ( void );
void hurstexponent( int );

int main(int argc, char** argv){
  int i, j, k, frame, lastframe, infocount, start2;
  double clengthnow, clength, dr2, dr[3], Ssum, *avgd;


  frame_max = 5000;
  get2dprojection = 1;


  //read and store N
  FILE *inp;
  inp = fopen("coarseningtraj.xyz","r");
  fscanf(inp, "%d", &N);
  printf("N: %d\n",N);
  fclose(inp);

  if (argc < 2){
    nn=100;
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
  P = ( double*** ) calloc ( frame_max, sizeof (double**) );
  Q = ( double*** ) calloc ( frame_max, sizeof (double**) );
  Pframe = ( double* ) calloc ( frame_max, sizeof (double) );
  unit = ( double** ) calloc ( N, sizeof (double*) );
  Pmax = ( double** ) calloc ( N, sizeof (double*) );
  for (i=0; i<N; i++){
    unit[i] = ( double* ) calloc (3 , sizeof ( double ));
    Pmax[i] = ( double* ) calloc (3 , sizeof ( double ));
  }
  d = ( double** ) calloc ( frame_max, sizeof (double*) );
  S = ( double* ) calloc ( frame_max, sizeof (double) );
  particles = ( int* ) calloc ( frame_max, sizeof (int) );
  avgd = ( double* ) calloc ( frame_max, sizeof (double) );
  info = ( double** ) calloc ( frame_max, sizeof (double*) );
  for (i=0; i<frame_max; i++){
    x[i] = ( double** ) calloc ( N, sizeof (double*)) ;
    P[i] = ( double** ) calloc ( N, sizeof (double*)) ;
    Q[i] = ( double** ) calloc ( N, sizeof (double*)) ;
    d[i] = ( double* ) calloc ( N, sizeof (double)) ;
    info[i] = ( double* ) calloc ( 2, sizeof (double)) ;
    for (j=0; j<N; j++){
      x[i][j] = ( double* ) calloc ( 3, sizeof (double) );
      P[i][j] = ( double* ) calloc ( 3, sizeof (double) );
      Q[i][j] = ( double* ) calloc ( 3, sizeof (double) );
    }
  }

  read_input();

  if (get2dprojection){
    nn = end - start ;
  }



  if (start == 0 && end == 0){
    printf("Start and end are both zero, null file\nTerminating simulation\n"); 
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

  if (get2dprojection){
    start2=start;
    printf("start2: %d\n",start2);
    hurstexponent(start2);
  }
  else {
    //# of times we will do hurst exponent calculation loop/function
    rem = int( double(end - start)/double(nn));
    printf("rem: %d\n", rem);

    printcounter = 0;
    for (i=0; i<rem; i++){
      start2=start + (i*nn);
      printf("start2: %d\n",start2);
      hurstexponent(start2);

      printcounter += 1;
    }
  }
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


void hurstexponent ( int start2 ){
  int i, j ,k, Pmaxframe;
  double dr2, dr[3], sqrtSD, check, length, avgdiam, Plength, max, min;
  FILE *otp;

  SD = 0.0;
  
  //get squared distance of the interval of interest
  for (i=0; i<N; i++){
    if (d[start2][i] == 0.0) continue;
    dr2 = 0.0;
    for (j=0; j<3; j++){
      dr[j] = x[start2+nn][i][j] - x[start2][i][j];

      dr2 += dr[j]*dr[j];
    }//end of j

    SD += dr2;
  }//end of i 
  //get length
  sqrtSD = sqrt(SD);

  //calculate unit vector
  for (i=0; i<N; i++){
    if (d[start2][i] == 0.0) continue;
    for (j=0; j<3; j++){
      unit[i][j] = (x[start2+nn][i][j] - x[start2][i][j]) / sqrtSD; 
    }
  }


  double dumsum;
  //this is for looping through every frame between start and start+nn
  max = -10.0;
  for (k=1; k<nn+1; k++){


    length = 0.0;
    for (i=0; i<N; i++){
      if (d[start2][i] == 0.0) continue;
      for (j=0; j<3; j++){
        Q[k-1][i][j] = x[start2+k][i][j] - x[start2][i][j];
        length += Q[k-1][i][j]*unit[i][j];
      }//end of j
    }//end of i

    dumsum = 0.0;
    Plength = 0.0;
    for (i=0; i<N; i++){
      if (d[start2][i] == 0.0) continue;
      for (j=0; j<3; j++){
        P[k-1][i][j] = Q[k-1][i][j] - (length * unit[i][j]); 
        Plength += (Q[k-1][i][j] - (length * unit[i][j]))*(Q[k-1][i][j] - (length * unit[i][j]));
        dumsum += P[k-1][i][j]*unit[i][j];
        
      }//end of j
    }//end of i

    if (Plength > max){
      max = sqrt(Plength);
      Pmaxframe = k-1;
    }

    if (dumsum > 0.0000000001){
      printf("dot product of p and unit vector is not 0.0!\nTerminating simulation\n");
      printf("sum: %1.16e\n", dumsum);
      exit(1);
    }
  // |---> not exactly 0.0 due to numerical error accumulance

  }//end of k (frame)

  printf("Pmax frame: %d + %d, Plength: %1.16lf, end-to-end: %1.16lf\n", start2, Pmaxframe, max,sqrtSD);

  if (get2dprojection){
    otp = fopen("2dprojection.dat","w");
    double xx,yy;
    for (k=1; k<nn+1; k++){
      xx = yy = 0.0;
      for (i=0; i<N; i++){
        if (d[start2][i] == 0.0) continue;
        for (j=0; j<3; j++){
          xx += Q[k-1][i][j] * unit[i][j];
          yy += P[k-1][i][j] * (P[Pmaxframe][i][j] / max);
          if (k-1 == Pmaxframe){
            printf("xx, yy = %1.16lf, %1.16lf\n", xx,yy);
          }//end of if
        }//end of j
      }//end of i
      fprintf(otp, "%1.16lf %1.16lf\n", xx, yy);
    }//end of k
    fclose(otp);
  }//end of if get2dprojection



  //calculate moments
  double sum, sum2;

  sum = 0.0;
  // first calculate P for each frame
  for (k=0; k<nn; k++){
    for (i=0; i<N; i++){
      if (d[start2][i] == 0.0) continue;
      dr2 = 0.0;
      for (j=0; j<3; j++){
        dr2 += P[k][i][j] * P[k][i][j];
      }//end of j
      sum += dr2;
    }//end of i
   
    Pframe[k] = sqrt(sum);   //**** P and Squared Distance are not particle averaged ****//
  }//end of k


  //get <a> --------> average radius
  int dumcount=0;
  avgdiam = 0.0;
  for (i=0; i<N; i++){
    if (d[start2][i] != 0.0){
      avgdiam += d[start2][i];
      dumcount += 1;
    }
  }//end of i
  avgdiam = (avgdiam / double(dumcount))/2.0;


  // MEAN
  sum = 0.0;
  for (k=0; k<nn; k++){
    sum += Pframe[k];
  }
  meanP = sum / double(nn);
  ////////
  

  // variance 
  sum = 0.0; 
  for (k=0; k<nn; k++){
    sum = (Pframe[k] - meanP)*(Pframe[k] - meanP);
  }
  variance = sum / double(nn);
  ///////////////

  // third moment
  sum = 0.0; 
  for (k=0; k<nn; k++){
    sum = (Pframe[k] - meanP)*(Pframe[k] - meanP)*(Pframe[k] - meanP);
  }
  sum = sum / double(nn);
  thirdP = sum / pow(variance, 1.5);
  ///////////////

  // fourth moment
  sum = 0.0; 
  for (k=0; k<nn; k++){
    sum = (Pframe[k] - meanP)*(Pframe[k] - meanP)*(Pframe[k] - meanP)*(Pframe[k] - meanP);
  }
  sum = sum / double(nn);
  fourthP = sum / pow(variance, 2.0); 
  ////////////////

  FILE  *otp2;
  if (printcounter == 0 ) otp = fopen ("hexpomoments.dat","w");
  else otp = fopen("hexpomoments.dat","a");
  fprintf(otp, "%1.16lf %1.16lf %1.16lf %1.16lf %1.16lf\n", SD, meanP, variance, thirdP, fourthP);
  fclose(otp);

  if (printcounter == 0) otp2 = fopen("pdistribution.dat","w");
  else otp2 = fopen("pdistribution.dat","a");
  for (i=0; i<nn; i++){
    fprintf(otp2, "%1.16lf %1.16lf %1.16lf\n", SD, Pframe[i], avgdiam);
  }
  fclose(otp2);

    
}

