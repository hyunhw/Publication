////////////////////////////////////////////////////////
//                                                    //
//  addhist.cpp ------> tpmhist_all.dat               //
//   |                                                //
//   |<-- pools all tpmhist-${i}.dat into one array   //
//                                                    //
////////////////////////////////////////////////////////


#include <cstdio>
#include <cstdlib>
#include <cmath>


int main(int argc, char** argv){
  int i, j, k, M, folders, **hist, dummy;
  FILE *inp, *otp;
  char name[40];

  M = 20;
  folders = 50;

  hist = (int**) calloc(M+1, sizeof(int*));
  for (i=0; i<M+1; i++){
    hist[i] = (int*) calloc (M+1, sizeof(int));
  }
  
  //looping over folders
  for (i=1; i<folders+1; i++){
    sprintf (name, "tpmhist-%d.dat", i);

    inp = fopen(name,"r");
    if (inp == NULL) continue;

    for (j=0; j<M+1; j++){
      for (k=0; k<M+1; k++){
        //scan for counts
        fscanf(inp, "%d ", &dummy);
        //add to total histogram
        hist[j][k] += dummy;
      }
      fscanf(inp, "\n");
    }

    fclose(inp);
  }

  otp = fopen ("tpmhist_all.dat", "w");
  for (i=0; i<M+1; i++){
    for (j=0; j<M+1; j++){
      fprintf(otp, "%d ", hist[i][j]);
    }
    fprintf(otp,"\n");
  }
  fclose(otp);


  return 0;

}
