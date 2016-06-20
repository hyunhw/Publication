////////////////////////////////////////////////////////////
//           |---------> scatterplot.dat                  //
//           |---------> forcestats.dat                   //
//           |---------> stressmsd.dat                    //
//           |---------> stressstats.dat                  //
//           |---------> stresstraj.dat                   //
//           |---------> xyztraj.dat                      //
//           |---------> particlemsd_timeavgd.dat         //
//  vanhove.cpp -------> xyzvanhove.dat                   //
//    |                      |---> stats for largest 60   //
//    |<--- zwindow.dat                                   //
//    |<--- coarseningtraj.xyz                            //
//    |<--- particlestress.dat                            //
//                                                        //
//    This code post-processes coarseningtraj.xyz to      //
//     calculate various measurements to construct        //
//     the van Hove distribution and other tests for      //
//     assessing the heterogeneity of the system          //
//                                                        //
////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cmath>

int N, temp_frame, nf, nft, start, end, maxN, *maxbubindex, *flag, tau;
double vol, box[3], boxh[3], ***x, **d, *sd, ***stress, ***force;

void read_file(void);
double pbc_vdr(double*, double*, double*);


int main (int argc, char** argv){
  int i, j, k, z, frame, lagtime, step;
  FILE *inp, *otp, *otp2;
  double dx[3], ds[3], dr[3], dr2, dij, dij2, rij;

  /*============= system parameters ===============*/

  temp_frame = 6000; //number of frames to allocate memory for x
  maxN = 60; //number of bubbles we want to keep track of
  tau = 1;
  
  /*============= system parameters ===============*/
  

  //open file to read in N and vol
  inp = fopen("coarseningtraj.xyz","r");
  if (inp==NULL){
    printf("Input file (coarseningtraj.xyz) is missing!\nTerminating simulation\n");
    exit(1);
  }
  fscanf(inp, "%d %lf", &N, &vol);
  fclose(inp);

  //determine box size
  for (i=0; i<3; i++){
    box[i] = pow(vol, 1.0/3.0);
    boxh[i] = box[i] / 2.0;
  }

  //scan window of simulation
  inp = fopen("zwindow.dat","r");
  if (inp==NULL){
    printf("Input file (zwindow.dat) is missing!\nTerminating simulation\n");
    exit(1);
  }
  fscanf(inp, "%d %d", &start, &end);
  fclose(inp);
  if (start == 0 & end == 0){
    printf("No window specified! (start=end=0)\nTerminating simulation\n");
    exit(1);
  }
  nf = end - start;

  //allocate memory
  x = (double***) calloc (temp_frame, sizeof (double**));
  stress = (double***) calloc (temp_frame, sizeof (double**));
  force = (double***) calloc (temp_frame, sizeof (double**));
  d = (double**) calloc (temp_frame, sizeof (double*));
  maxbubindex = (int*) calloc (maxN, sizeof (int));
  flag = (int*) calloc (N, sizeof (int));
  sd = (double*) calloc (maxN, sizeof (double));
  for (i=0; i<temp_frame; i++){
    x[i] = (double**) calloc (N, sizeof (double*));
    stress[i] = (double**) calloc (N, sizeof (double*));
    force[i] = (double**) calloc (N, sizeof (double*));
    d[i] = (double*) calloc (N, sizeof (double));
    for (j=0; j<N; j++){
      x[i][j] = (double*) calloc (3, sizeof (double));
      stress[i][j] = (double*) calloc (3, sizeof (double));
      force[i][j] = (double*) calloc (3, sizeof (double));
    }
  }

  //read in diameter and trajectory
  read_file();



  //////////////
  //          //
  //   MAIN   //
  //          //
  //////////////


  /////////////////////
  // Biggest bubbles //
  /////////////////////
  
  //determine biggest 40 particles that survived at the frame = end
  
  double dummylf;
  for (i=0; i<N; i++){
    flag[i] = 0;
  }
  for (k=0; k<maxN; k++){     // k - index of 40 particles
    dummylf = 0.0;
    for (i=0; i<N; i++){
      if (d[end][i] == 0.0) continue;
      if (d[end][i] == 2.5) continue;
      if (flag[i] == 1) continue;

      if (d[end][i] > dummylf){
        maxbubindex[k] = i;
        dummylf = d[end][i];
      }
    }//end of i
    flag[maxbubindex[k]] = 1;
  }//end of k

  for (k=0; k<maxN; k++){
    printf("max bub %d and diameter: %lf\n", maxbubindex[k], d[end][maxbubindex[k]]);
  }


  /////////////////////////
  // van Hove stats calc //
  /////////////////////////

  otp = fopen("vanhovestats.dat","w");

  //printing the first three lines
  //  "-1" denotes column of time
  fprintf(otp, "-1 ");
  for (k=0; k<(maxN*3); k++){
    fprintf(otp, "%d ", tau); 
  }
  fprintf(otp, "\n");  //------> printing specified lag time
  fprintf(otp, "-1 ");
  for (k=0; k<maxN; k++){
    fprintf(otp, "%d %d %d ", maxbubindex[k], maxbubindex[k], maxbubindex[k]); 
  }
  fprintf(otp, "\n");  //------> printing specified bubble index
  fprintf(otp, "-1 ");
  for (k=0; k<maxN; k++){
    fprintf(otp, "%1.16lf %1.16lf %1.16lf ", d[end][maxbubindex[k]], d[end][maxbubindex[k]], d[end][maxbubindex[k]]); 
  }
  fprintf(otp, "\n");  //------> printing bubble diameter


  //printing actual dx dy dz
  for (frame=start; frame<end; frame++){
    fprintf(otp, "%d ", frame);

    for (k=0; k<maxN; k++){
      for (j=0; j<3; j++){
        dx[j] = x[frame][maxbubindex[k]][j] - x[frame-tau][maxbubindex[k]][j];
        if (dx[j] > boxh[j]) dx[j] -= box[j];
        else if (dx[j] < -boxh[j]) dx[j] += box[j];
      }//end of j
      fprintf(otp, "%1.16lf %1.16lf %1.16lf ", dx[0], dx[1], dx[2]); 
    }//end of k
    fprintf(otp, "\n");
  }//end of frame

  fclose(otp);
  

  ////////////////
  // Force calc //
  ////////////////
  

  otp = fopen("forcestats.dat","w");
  otp2 = fopen("totalforcestats.dat","w");
  double Fpre, Ftotal, Fframei, avgF;
  int fcount;

  //printing force stats 
  for (frame=start; frame<end; frame++){
    fprintf(otp, "%d ", frame);
    fprintf(otp2, "%d ", frame);
      fcount = 0;
      avgF = 0.0;
 
    for (k=0; k<N; k++){
      Ftotal = 0.0;
      if (d[frame][k] == 0.0) continue;
      for (j=0; j<N; j++){
        if (d[frame][j] == 0.0) continue;
        if (j <= k) continue;

        dr2 = pbc_vdr (x[frame][k], x[frame][j], dr);
        dij = ( d[frame][j] + d[frame][k] ) /2.0;
        dij2 = dij * dij;

        if (dr2 > dij2) continue;
        rij = sqrt(dr2);
        Fpre = (2.0 * (1.0/(rij*dij) - 1.0/dij2));

        for (z=0; z<3; z++){
          Ftotal += (Fpre*dr[z]) * (Fpre*dr[z]);
          force[frame][k][z] += (Fpre * dr[z]);
        }
        Fframei = sqrt(Ftotal);
        avgF+= Fframei;
      }//end of j
      fprintf(otp, "%1.16lf ", Fframei); 
      fcount +=1;
    }//end of k
    fprintf(otp2, "%1.16lf\n", avgF/double(fcount));
  }//end of frame

  fclose(otp);
  fclose(otp2);
  



  //////////////////////
  // stress traj calc //
  //////////////////////

  otp = fopen("stressstats.dat","w");

  //printing the first three lines
  //  "-1" denotes column of time
  fprintf(otp, "-1 ");
  for (k=0; k<(maxN*3); k++){
    fprintf(otp, "%d ", tau); 
  }
  fprintf(otp, "\n");  //------> printing specified lag time
  fprintf(otp, "-1 ");
  for (k=0; k<maxN; k++){
    fprintf(otp, "%d %d %d ", maxbubindex[k], maxbubindex[k], maxbubindex[k]); 
  }
  fprintf(otp, "\n");  //------> printing specified bubble index
  fprintf(otp, "-1 ");
  for (k=0; k<maxN; k++){
    fprintf(otp, "%1.16lf %1.16lf %1.16lf ", d[end][maxbubindex[k]], d[end][maxbubindex[k]], d[end][maxbubindex[k]]); 
  }
  fprintf(otp, "\n");  //------> printing bubble diameter


  //printing actual dx dy dz
  for (frame=start; frame<end; frame++){
    fprintf(otp, "%d ", frame);

    for (k=0; k<maxN; k++){
      for (j=0; j<3; j++){
        ds[j] = stress[frame][maxbubindex[k]][j] - stress[frame-tau][maxbubindex[k]][j];
        //new***  take out volume term for vanhove stress and stress traj analysis to see if anything changes
        ds[j] = ds[j] * d[frame][maxbubindex[k]] * d[frame][maxbubindex[k]] * d[frame][maxbubindex[k]] / 8.0; 
      }//end of j
      fprintf(otp, "%1.16lf %1.16lf %1.16lf ", ds[0], ds[1], ds[2]); 
    }//end of k
    fprintf(otp, "\n");
  }//end of frame

  fclose(otp);



  /////////////////////
  // stress msd calc //
  /////////////////////


  int count;
  double sum1, sum2, sum3;
  otp2 = fopen ("stressmsd.dat","w");

  for (lagtime=1 ; lagtime<nf; lagtime++){
    fprintf(otp2, "%d ", lagtime);
    count = 0;
    sum1 = sum2 = sum3 = 0.0;
    
    for (step = lagtime; step<nf; step++){
      for (k=0; k<maxN; k++){
        for (j=0; j<3; j++){
          ds[j] = fabs(stress[start+step][maxbubindex[k]][j] - stress[start+step-lagtime][maxbubindex[k]][j]);
          ds[j] = ds[j]*d[start+step][maxbubindex[k]]*d[start+step][maxbubindex[k]]*d[start+step][maxbubindex[k]]/8.0;
        }//end of j
        sum1 += ds[0]*ds[0];
        sum2 += ds[1]*ds[1];
        sum3 += ds[2]*ds[2];
      }//end of k
      count++;
    }//end of step
    sum1 = sum1 / double(count*maxN);
    sum2 = sum2 / double(count*maxN);
    sum3 = sum3 / double(count*maxN);
    fprintf(otp2, "%1.16lf %1.16lf %1.16lf\n", sum1, sum2, sum3);
  }//end of lagtime
  fclose(otp2);


  //////////////////////
  // MSD  calculation //  ---> yes time origin (prints dr2) 
  //////////////////////

  otp = fopen ("particlemsd_timeavgd.dat","w");
  
  //printing the first two lines
  //  "-1" denotes column of time
  fprintf(otp, "-1 ");
  for (k=0; k<maxN; k++){
    fprintf(otp, "%d ", maxbubindex[k]); 
  }
  fprintf(otp, "\n");  //------> printing specified bubble index
  fprintf(otp, "-1 ");
  for (k=0; k<maxN; k++){
    fprintf(otp, "%1.16lf ", d[end][maxbubindex[k]]); 
  }
  fprintf(otp, "\n");  //------> printing bubble diameter

  for (lagtime=1; lagtime<nf; lagtime++){
    fprintf (otp, "%d ", lagtime);
    
    for (k=0; k<maxN; k++){
      sd[k] = 0.0;
    }
    count = 0;

    for (step=lagtime; step<nf; step++){
      for (k=0; k<maxN; k++){
        dr2 = 0.0;
        for (j=0; j<3; j++){

          dx[j] = x[start+step][maxbubindex[k]][j] - x[start+step-lagtime][maxbubindex[k]][j]; 
          if (dx[j] > boxh[j]) dx[j] -= box[j];
          else if (dx[j] < -boxh[j]) dx[j] += box[j];

          dr2 += dx[j] * dx[j];
        }//end of j
        sd[k] += dr2;
      }//end of k
      count += 1;
    }//end of step 
    for (k=0; k<maxN; k++){
      sd[k] = sd[k] / double(count);
      fprintf(otp, "%1.16lf ", sd[k]);
    }
    fprintf(otp, "\n");
  }//end of lagtime

  fclose(otp);
  /*

  //calculate the slope of the ensemble/time averaged MSD
  inp = fopen("particlemsd_timeavgd.dat","r");
  if (inp == NULL){
    printf("input file 'particlemsd_timeavgd.dat' is missing!\nTerminating simulation\n");
    exit(1);
  }
  otp = fopen("ensemtimeMSD.dat","w");

  double sum;
  int dummyint;
  char tt[40000];
  for (i=0; i<nf+1; i++){
   
    if (i<2){
      fgets(tt, 40000, inp);
    }

    else {
      sum = 0.0;
      fscanf(inp, "%d ", &dummyint);

      //add MSD over particles
      for (k=0; k<maxN; k++){
        fscanf(inp, "%lf ", &dummylf);
        sum += dummylf;
      }//end of k
      sum = sum / double(maxN);
      fprintf(otp, "%d %1.16lf\n", dummyint, sum);
    }
  }//end of nf
  */

  /////////////////////
  // scatterplot.dat //
  /////////////////////

  FILE *otp3;
  double newsxy, newsxz, newsyz;
  otp3 = fopen("scatterplot.dat","w");
  for (frame=start; frame<end; frame++){
    for (k=0; k<maxN; k++){
      dx[0] = x[frame][maxbubindex[k]][0] - x[frame-1][maxbubindex[k]][0]; 
      if (dx[0] > boxh[0]) dx[0] -= box[0];
      else if (dx[0] < -boxh[0]) dx[0] += box[0];
      dx[1] = x[frame][maxbubindex[k]][1] - x[frame-1][maxbubindex[k]][1]; 
      if (dx[1] > boxh[1]) dx[1] -= box[1];
      else if (dx[1] < -boxh[1]) dx[1] += box[1];
      dx[2] = x[frame][maxbubindex[k]][2] - x[frame-1][maxbubindex[k]][2]; 
      if (dx[2] > boxh[2]) dx[2] -= box[2];
      else if (dx[2] < -boxh[2]) dx[2] += box[2];

      dr2 = (dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]);
      ds[0] = stress[frame][maxbubindex[k]][0] - stress[frame-1][maxbubindex[k]][0]; 
      ds[1] = stress[frame][maxbubindex[k]][0] - stress[frame-1][maxbubindex[k]][1]; 
      ds[2] = stress[frame][maxbubindex[k]][0] - stress[frame-1][maxbubindex[k]][2]; 
      newsxy = dx[0] * d[frame][maxbubindex[k]] * d[frame][maxbubindex[k]] * d[frame][maxbubindex[k]];
      newsxz = dx[1] * d[frame][maxbubindex[k]] * d[frame][maxbubindex[k]] * d[frame][maxbubindex[k]];
      newsyz = dx[2] * d[frame][maxbubindex[k]] * d[frame][maxbubindex[k]] * d[frame][maxbubindex[k]];
      //frame, diameter, stress 0, stress 1, stress 2, stress 0d3, stress 1d3, stress 2d3, delx, dely, delz, delr2
      fprintf(otp3, "%d %lf %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf %1.16lf\n", frame, d[frame][maxbubindex[k]], ds[0], ds[1], ds[2], newsxy, newsxz, newsyz, dx[0], dx[1], dx[2], dr2); 

    }//end of k
  }//end of frame
  fclose(otp3);
  return 0;
}

void read_file(){
  int i, j, k, dummyint;
  double dummylf, dummylf2, dummylf3;
  FILE *inp;

  inp = fopen("coarseningtraj.xyz","r");
  if (inp == NULL){
    printf("Input file 'coarseningtraj.xyz' is missing!\nTerminating simulation\n");
    exit(1);
  }

  nft = 0;
  while (!feof(inp)){
    fscanf(inp, "%d %lf\n\n", &dummyint, &dummylf);
    for (i=0; i<N; i++){
      fscanf(inp, "H %lf %lf %lf %lf\n", &x[nft][i][0], &x[nft][i][1], &x[nft][i][2], &d[nft][i]);
    }//end of i
    nft += 1;
  }//end of while loop
  fclose(inp);


  inp = fopen("particlestress.dat","r");
  if (inp == NULL){
    printf("Input file 'particlestress.dat' is missing!\nTerminating simulation\n");
    exit(1);
  }

  nft = 0;
  while (!feof(inp)){
    fscanf(inp, "%d\n\n", &dummyint);
    for (i=0; i<N; i++){
      fscanf(inp, "%lf %lf %lf %lf %lf %lf\n", &dummylf, &dummylf2, &dummylf3, &stress[nft][i][0], &stress[nft][i][1], &stress[nft][i][2]);
    }//end of i
    nft += 1;
  }//end of while loop
  fclose(inp);
}

double pbc_vdr (double x1[3], double x2[3], double dr[3]){
  double mdr2=0.0;
  int i;

  //PBC in all directions
  for (i=0; i<3; i++){
    dr[i] = x1[i] - x2[i];

    if (dr[i]>boxh[i]) dr[i] -= box[i];
    else if (dr[i]<-boxh[i]) dr[i] += box[i];
  
    mdr2 += dr[i]*dr[i];
  }

  return mdr2;

}



