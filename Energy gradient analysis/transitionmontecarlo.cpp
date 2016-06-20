/////////////////////////////////////////////////////////////////
//                                                             //
//  transitionmontecarlo.cpp   -------> distributions..        //
//    ^                                                        //
//    |-- tpmatrix.dat                                         //
//                                                             //
//   This file:                                                //
//     reads in the transition probability matrix and          //
//     carries out a monte carlo according to the tpm          //
//                                                             //
/////////////////////////////////////////////////////////////////


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#define PI 3.141592653589793238462643383


/////////////////////////////
// define global variables //
int size, step_max, MCstep, MCcount;
long int idum;
double **tpm;


/////////////////////////////
// define global functions //
double ran2 (void);

int main (int argc, char** argv){
  int i, j, k, dummy;
  FILE *inp, *otp;
  time_t timer; //type long int

  //MC step max
  step_max = 100000;
  
  size = 0;
  //get and store current time
  time(&timer);
  idum = -timer;
  printf("seed: %ld\n", idum);



  inp = fopen("tpmatrix_14_5_11.dat","r");
  fscanf(inp, "%d", &size);
  if (size==0){
    printf("matrix file does not include the size of the matrix\nterminating simulation");
    exit(1);
  }
  fclose(inp);

  printf("tpm: %d x %d matrix\n", size, size);

  tpm = (double**) calloc (size, sizeof(double*));
  for (i=0; i<size; i++){
    tpm[i] = (double*) calloc(size, sizeof(double));
  }

  inp = fopen("tpmatrix_14_5_11.dat","r");
  // i - row
  for (i=0; i<size+1; i++){
    // j - column
    for (j=0; j<size; j++){
      if (i==0) {
        fscanf(inp, "%d\n", &dummy);
        continue;
      }
      else if (j==(size-1)){
        fscanf(inp, "%lf\n", &tpm[i-1][j]);
      }
      else {
        fscanf(inp, "%lf ", &tpm[i-1][j]);
      }
    }//end of j
  }//end of i
  fclose(inp);

/*
  //test
  FILE *test;
  test = fopen("test.dat","w");
  for (i=0; i<size; i++){
    for (j=0; j<size; j++){
      fprintf(test,"%1.16lf ", tpm[i][j]);
    }
    fprintf(test,"\n");
  }
  fclose(test);
*/


  /////////////////////////////////// 
  // Start of Monte Carlo sampling // ---> column sampling
  /////////////////////////////////// 

  //set so that it samples column (reading in row normalized matrix)

  int start, end;
  FILE *otp2, *otp3;
  double gamma, cumulant;
  otp = fopen("exitdistribution.dat","w");
  otp2 = fopen("countdistribution.dat","w");

  // INITIALLY pick a random index to start (pick a row)
  
  //makes sure that start can never be 0
  start = 0;
  while (start == 0){
    start = int(ran2() * size);
    printf("start: %d\n", start);
  }

  MCstep = 0;
  MCcount = 0;

  while (MCstep < step_max){

    cumulant = 0.0;
    
    //generate another uniform random # between [0,1]
    // this value is used to determin the exit bin
    gamma = ran2();
    printf("gamma: %lf\n", gamma);

    //value of start is the value of the row I want to scan thru to pick the exit point
    // so i loop thru the column (scanning across row)
    for (i=0; i<size; i++){
      printf("i: %d\n", i);
      cumulant += tpm[start][i]; 
      printf("cumulant: %lf\n", cumulant);
      end = i;

      //cumulant is greater than RNG, exit the matrix
      if (cumulant >= gamma) break;
    }

    //fprintf(otp, "%d\n", end);
    //fprintf(otp3, "%lf\n", cumulant);

    //if we enter the exiting column
    if (end == 0){
      start = end;
      //update MC step (which keeps track of how many "enter->exit" events happened)
      MCstep += 1;

      fprintf(otp, "%d\n", end); //print index of exit
      fprintf(otp2, "%d\n", MCcount);

      //reset MC count (which keeps track of how many MC steps occured between a single MCstep (defined above))
      MCcount = 0;

    }//end of if start == 1

    else {
      start = end;
      fprintf(otp, "%d\n", end); //print index of exit
      MCcount += 1;
    }
    //printf("start: %d\n", start);


  }//end of while loop


  fclose(otp);
  fclose(otp2);




/*
  /////////////////////////////////// 
  // Start of Monte Carlo sampling // ----> row sampling
  /////////////////////////////////// 



  //set so that it samples row (reading in column normalized matrix)

  int start, end;
  FILE *otp2, *otp3;
  double gamma, cumulant;
  otp = fopen("exitdistribution.dat","w");
  otp2 = fopen("countdistribution.dat","w");

  // INITIALLY pick a random index to start (pick a row)
  
  //makes sure that start can never be 0
  start = 0;
  while (start == 0){
    start = int(ran2() * size);
    printf("start: %d\n", start);
  }

  MCstep = 0;
  MCcount = 0;

  while (MCstep < step_max){

    cumulant = 0.0;
    
    //generate another uniform random # between [0,1]
    // this value is used to determin the exit bin
    gamma = ran2();
    printf("gamma: %lf\n", gamma);

    //value of start is the value of the column I want to scan thru to pick the exit point
    // so i loop thru the row (scanning across column)
    for (i=0; i<size; i++){
      printf("i: %d\n", i);
      cumulant += tpm[i][start]; 
      printf("cumulant: %lf\n", cumulant);
      end = i;

      //cumulant is greater than RNG, exit the matrix
      if (cumulant >= gamma) break;
    }

    //fprintf(otp, "%d\n", end);
    //fprintf(otp3, "%lf\n", cumulant);

    //if we enter the exiting column
    if (end == 0){
      start = end;
      //update MC step (which keeps track of how many "enter->exit" events happened)
      MCstep += 1;

      fprintf(otp, "%d\n", end); //print index of exit
      fprintf(otp2, "%d\n", MCcount);

      //reset MC count (which keeps track of how many MC steps occured between a single MCstep (defined above))
      MCcount = 0;

    }//end of if start == 1

    else {
      start = end;
      fprintf(otp, "%d\n", end); //print index of exit
      MCcount += 1;
    }
    //printf("start: %d\n", start);


  }//end of while loop


  fclose(otp);
  fclose(otp2);









*/



/*

sadfafsadf



  //set so that it samples row (reading in row normalized matrix)

  int start, end;
  FILE *otp2, *otp3;
  double gamma, cumulant;
  otp = fopen("exitdistribution.dat","w");
  otp2 = fopen("countdistribution.dat","w");

  // INITIALLY pick a random index to start (pick a row)
  
  //makes sure that start can never be 0
  start = 0;
  while (start == 0){
    start = int(ran2() * size);
    printf("start: %d\n", start);
  }

  MCstep = 0;
  MCcount = 0;


  while (MCstep < 500000){


    cumulant = 0.0;
    
    //generate another uniform random # between [0,1]
    // this value is used to determin the exit bin
    gamma = ran2();
    printf("gamma: %lf\n", gamma);

    //value of start is the value of the row I want to scan thru to pick the exit point
    // so i loop thru the column 
    for (i=0; i<size; i++){
      printf("i: %d\n", i);
      cumulant += tpm[i][start]; 
      printf("cumulant: %lf\n", cumulant);
      end = i;

      //cumulant is greater than RNG, exit the matrix
      if (cumulant >= gamma) break;
    }

    //fprintf(otp, "%d\n", end);
    //fprintf(otp3, "%lf\n", cumulant);

    //if we enter the exiting column
    if (start == 1){
      start = end;
      //update MC step (which keeps track of how many "enter->exit" events happened)
      MCstep += 1;

      fprintf(otp, "%d\n", end); //print index of exit
      fprintf(otp2, "%d\n", MCcount);

      //reset MC count (which keeps track of how many MC steps occured between a single MCstep (defined above))
      MCcount = 0;

    }//end of if start == 1

    else {
      start = end;
      MCcount += 1;
    }
    //printf("start: %d\n", start);


  }//end of while loop


  fclose(otp);
  fclose(otp2);
*/

  return 0;

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

//fenerates random number from [0,1]
double ran2 (void) {
  int j;
  long int k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  double tempp;
  extern long idum;

  if (idum <= 0){
    if (-(idum) < 1) idum=1;
    else idum = -(idum);
    idum2 = (idum);
    for (j=NTAB + 7; j>=0; j--){
      k=(idum)/IQ1;
      idum = IA1 * (idum - k *IQ1) - k * IR1;
      if (idum <0) idum += IM1;
      if (j<NTAB) iv[j] = idum;
    }

    iy = iv[0];
  }

  k = (idum) / IQ1;
  idum = IA1 * (idum - k * IQ1) - k * IR1;
  if (idum<0) idum += IM1;
  k = idum2 / IQ2;
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
  if (idum2 < 0) idum2 += IM2;
  j = iy/NDIV;
  iy = iv[j] - idum2;
  iv[j] = idum;
  if (iy<1) iy += IMM1;
  if ((tempp = AM*iy) > RNMX ) return RNMX;
  else return tempp;

}
