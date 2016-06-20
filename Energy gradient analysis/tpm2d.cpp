/////////////////////////////////////////////////////////////////////
//                                                                 //
//   tpm2d.cpp --------> tpmhist.dat                               //
//     |                                                           //
//     |<---- tpm.dat, Xmin, Delta, m (bin number) specified       //
//                                                                 //
//    This file reads in tpm.dat and constructs the transition     //
//    probability matrix by sorting numbers into ordered pairs     //
//    and accumulating a 2d histogram based on specifications      //
//                                                                 //
/////////////////////////////////////////////////////////////////////


#include <cstdio>
#include <cstdlib>
#include <cmath>


int main (int arc, char** argv){
  int i, j, k, count, count2, exitcount, M, index1, index2, **hist, flag, mem, mode; 
  double value1, value2, dummy, dummyc, start, end, first, last, dummyend, min, max, delta, **st, amin, amex;
  //dummy end stores the last value, to print to exit

  FILE *inp, *inp2, *inp3, *otp, *otp2, *otp3, *otp4;

  /*-- Specifications --*/
  M = 26;
  min = 0.00;
  max = 0.0545; 
  delta = (max-min) / double(M);
  mode = 1;
  //min and max boundary for <a>
  amin = 0.85;
  amex = 0.95;
  /*--------------------*/

  //read zwindow.dat to skip runs that have no windows assigned
  inp2 = fopen("zwindow.dat","r");
  fscanf(inp2, "%lf %lf", &start, &end);
  printf("start: %lf end: %lf\n", start, end);
  if (start == 0.0 && end == 0.0){
    printf("window is 0 and 0!\n"); fflush(stdout);
    exit(1);
  }
  fclose(inp2);
  int range;
  range = int((end - start)/2);
  //start = start + range;
  //end = start + range;

  if ( mode == 0 ){
    otp = fopen("tpmpairs_withoutfirstpair.dat","w");
    otp2 = fopen("tpmpairs_enter.dat","w");
    otp3 = fopen("tpmpairs_exit.dat","w");
    inp = fopen("tpm.dat","r");


    mem = end - start + 1;
    st = (double**) calloc (mem, sizeof (double*));
    for (i=0; i<mem; i++){
      st[i] = (double*) calloc (3, sizeof (double));
    }
 

    inp3 = fopen("windowstats.dat","r");
    //read in stats from windowstats.dat and calculate <a> / Eb
    int s1;
    double s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16;
    for (i=0; i<mem; i++){
      fscanf(inp3, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &s1, &s2, &s3, &s4, &s5, &s6, &s7, &s8, &s9, &s10, &s11, &s12, &s13, &s14, &s15, &s16);
      st[i][0] = double(s1); 
      st[i][1] = (s15/2.0) / (s16 / s13);
      st[i][2] = s15/2.0; 
      //printf("frame: %lf <a>/Eb: %lf\n", st[i][0], st[i][1]);
      /*
      if ((i+start) != s1) {
        printf ("Incorrect reading of windowstats\nTerminating simulation\n"); fflush(stdout);
        exit(1);
      }
      */
    }
    fclose(inp3);
    printf("read window stats\n"); fflush(stdout);

    printf("mem: %d\n", mem); fflush(stdout);
//    for (i=0; i<mem; i++){
//      printf("%lf\n", st[i][1]);
//    }


    flag = 0;
    dummyend = 0.0;
    exitcount = 0;

    //while it isn't the end of file
    while (!feof(inp)){
      fscanf(inp, "%lf ", &dummy);
      //printf("%lf\n", dummy); fflush(stdout);

      //if we read the coarsening number, don't do anything but reset count to 0
      if (fmod(dummy, 1.0) == 0.0){
        count = 0;
        count2 = 0;

        if (-dummy < start || -dummy > end) flag = 0;
        if (-dummy < start || -dummy > end) continue;
        if (-dummy >= start || -dummy <= end){
          flag = 1;
          dummyc = fabs(dummy);
        }
        //this is done at the next coarsening step
        if (exitcount > 0){
          if (flag == 1){
            if (st[int(dummyc-start)][2] >= amin && st[int(dummyc-start)][2] <= amex) {
              fprintf(otp3, "%1.16lf\n", sqrt(fabs(dummyend)/**st[int(dummyc-start)][1]*/)/*, st[int(dummyc-start)][2]*/);
            }
          }
        }
        //exit count is set to 0 when i read the very first coarsening number of file
        exitcount = 0;
      }
      else {
        count += 1;
        count2 += 1;

        if ((count%2)==1) value1 = dummy;
        //count2 !=2 condition skips the first pair
        else if ((count%2)==0 && count2 != 2){
          value2 = dummy; 
          dummyend = dummy;
          exitcount += 1;
          if (flag == 1){
            if (st[int(dummyc-start)][2] >= amin && st[int(dummyc-start)][2] <= amex ) {
              fprintf(otp,"%1.16lf %1.16lf\n", sqrt(fabs(value1)/**st[int(dummyc-start)][1]*/), sqrt(fabs(value2)/**st[int(dummyc-start)][1]*/)/*,st[int(dummyc-start)][2]*/);
            }
          }
          //onlny prints enter value of the pair after skipped pair
          if (count2 == 4){
            if (flag ==1){
              if (st[int(dummyc-start)][2] >= amin && st[int(dummyc-start)][2] <= amex ) {
                fprintf(otp2, "%1.16lf\n", sqrt(fabs(value1)/*8st[int(dummyc-start)][1]*/)/*, st[int(dummyc-start)][2]*/);
              }
            }
          }
        }
      }
    }

    fclose (inp);
    fclose (otp);
    fclose (otp2);
    fclose (otp3);

  } //end of if (mode == 0)

  /// if mode isnt 0
  else { 

    hist = (int**) calloc (M+1, sizeof (int*));
    for (i=0; i<(M+1); i++){
      hist[i] = (int*) calloc (M+1, sizeof (int));
    }



    ////////////
    //  Body  //
    ////////////
   
    inp = fopen("tpmpairs_withoutfirstpair.dat","r");
    while (!feof(inp)) {
      fscanf(inp, "%lf %lf\n", &value1, &value2);
      index1 = int((value1-min) / delta) + 1;
      index2 = int((value2-min) / delta) + 1;
      hist[index1][index2] += 1;
    }
    fclose(inp);


    /////////////
    //  Enter  //
    /////////////
    
    inp2 = fopen("tpmpairs_enter.dat","r");
    while (!feof(inp2)) {
      fscanf(inp, "%lf\n", &value2);
      index2 = int((value2-min) / delta) + 1;
      hist[0][index2] += 1;
    }
    fclose(inp2);
 

    ////////////
    //  Exit  //
    ////////////
    
    inp3 = fopen("tpmpairs_exit.dat","r");
    while (!feof(inp3)) {
      fscanf(inp, "%lf\n", &value1);
      index1 = int((value1-min) / delta) + 1;
      hist[index1][0] += 1;
    }
    fclose(inp3);


/*
    flag = 0;
    dummyend = 0.0;
    exitcount = 0;
    //while it isn't the end of file ..
    while (!feof(inp)){
      fscanf(inp, "%lf ", &dummy); 

      //if we read the coarsening number, don't do anything but reset count to 0
      // if the coarsening number we read does not fall into the zwindow range, continue
      if (fmod(dummy,1.0) == 0.0){
        count = 0;
        count2 = 0;

        if (-dummy < start || -dummy > end) flag = 0;
        if (-dummy < start || -dummy > end) continue;
        if (-dummy >= start || -dummy <= end){
          flag = 1;
          dummyc = fabs(dummy);
          printf("dummyc: %lf\n", dummyc); fflush(stdout);
        }

        //this is done at the end of every coarsening window to store the "end" values
        if (exitcount > 0){
         // fprintf(otp3, "%1.16lf\n", dummyend); 
          //because we are dealing with the exit value (which is the 1st column of the M+1 x M+1 matrix,
          // we only calculate the i value (index1) and assign j=0
          index1 = int ((sqrt(fabs(dummyend)*st[int(dummyc-start)][1])-min) / delta) + 1;
          ///////////hist[index1][0] += 1;
        }
        //exit count is set to 0 when I read the very first coarsening number of file
        exitcount = 0;
      }
      // if we are reading the gradient
      else if (flag == 1) {
        count += 1;
        count2 += 1;
        if ((count%2)==1) value1 = dummy;

        if (dummyc == start){
          printf ("value1: %lf <a>/Eb: %lf\n", value1, st[int(dummyc-start)][1]); fflush(stdout);
        }

        //count2!=2 condition skips the first pair
        else if ((count%2)==0 && count2 != 2) {
          value2 = dummy;
          dummyend = dummy;
          exitcount += 1;
          //fprintf(otp,"%1.16lf %1.16lf\n", value1, value2);
          index1 = int((sqrt(fabs(value1)*st[int(dummyc-start)][1])-min) / delta)+1;
          index2 = int((sqrt(fabs(value2)*st[int(dummyc-start)][1])-min) / delta)+1;
          //////////////////hist[index1][index2] += 1;
          //only prints enter value of the pair after skipped pair (prints once)
          if (count2 == 4){
            //printing first value
            //fprintf(otp2, "%1.16lf\n", value1);
            index2 = int((sqrt(fabs(value1)*st[int(dummyc-start)][1])-min) / delta)+1;
            //////////////hist[0][index2] += 1;
          }
        }
      }
    }
*/

    otp4 = fopen("tpmhist.dat","w");
    for (i=0; i<M+1; i++){
      for (j=0; j<M+1; j++){
        fprintf(otp4, "%d ", hist[i][j]);
      }
      fprintf(otp4, "\n");
    }
    fclose(otp4);
    
  }//end of else

  return 0;

}
