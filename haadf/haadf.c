#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double pi=3.1415926;

/* norm distribution in 2D */
double multivariate(double cx, double cy, double sigma,
		    double x, double y){
  double c = (double)1.0/2/pi;
  double f1=exp(-0.5*(x-cx)*(x-cx)/sigma/sigma)/sigma;
  double f2=exp(-0.5*(y-cy)*(y-cy)/sigma/sigma)/sigma;
  return c*f1*f2;
}

int dot(double mat[3][3], double vec1[3], double vec2[3]){
  int i,j;
  for (i=0;i<3;i++){
    vec2[i]=0;
    for (j=0;j<3;j++)
      vec2[i]+=mat[j][i]*vec1[j];
  }
  return 0;
}

int main(){
  int i,j,k;
  int px_num, nsamples, isdepth, na;
  double px_sz, zmin, zmax;
  double proj_mat[3][3], lats[3][3], sigma;
   
  FILE *fid=fopen("HAADF.in", "r");
  fscanf(fid, "%d", &px_num);
  fscanf(fid, "%lf", &px_sz);
  fscanf(fid, "%lf", &sigma);
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      fscanf(fid, "%lf", &(proj_mat[i][j]));
  for (i=0;i<3;i++){
    zmin=0;
    for (j=0;j<3;j++)
      zmin+=proj_mat[i][j]*proj_mat[i][j];
    for (j=0;j<3;j++)
      proj_mat[i][j]/=sqrt(zmin);
  }
  for (i=0;i<3;i++){
    for (j=0;j<i;j++){
      zmin=proj_mat[i][j];
      proj_mat[i][j]=proj_mat[j][i];
      proj_mat[j][i]=zmin;
    }
  }
  for (i=0;i<3;i++){
    for (j=0;j<3;j++)
      printf("%lf ", (proj_mat[i][j]));
    printf("\n");
  }

  fscanf(fid, "%lf %lf", &zmin, &zmax);
  if (zmin>zmax) isdepth=0; else isdepth=1;
  fscanf(fid, "%d", &nsamples);

  int px_ct=(int)(px_num/2);
  int smear_sz=(int)(sigma/px_sz*3);

  double **canvas;
  canvas = (double **) malloc(sizeof(double*)*px_num);
  for (i=0;i<px_num;i++)
    canvas[i]=(double *) malloc(sizeof(double)*px_num);
  for (i=0;i<px_num;i++)
    for (j=0;j<px_num;j++)
      canvas[i][j]=0;
  
  char fxyz[100] = { '\0' };
  fgets(fxyz, sizeof(fxyz), fid);
  fgets(fxyz, sizeof(fxyz), fid);
  size_t ln = strlen(fxyz) - 1;
  if (*fxyz && fxyz[ln] == '\n') 
    fxyz[ln] = '\0';
 
  fclose(fid);
  fid = fopen(fxyz, "r");
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      fscanf(fid, "%lf", &(lats[i][j]));
  fscanf(fid, "%d", &na);
  fclose(fid);
  int lins=na+4;

  fid = fopen(fxyz, "r");
  int lino=0,ncount=0, atype, atom_id;
  double dvec[3], ivec[3], rvec[3], pvec[3];
  int pb_x, pb_y, pb_z, px_nx, px_ny, ext_x, ext_y, paint_x, paint_y;
  printf("Processing");
  while (~feof(fid)){
    if (lino % lins <=3) {
      fgets(fxyz, 100, fid);
      lino++;
      continue;
    }
    ncount++;
    printf("\r%4d out of %4d", ncount, nsamples);
    if ( ncount > nsamples ) break;
    for (atom_id=0;atom_id<na;atom_id++){
      lino++;
      fscanf(fid, "%lf %lf %lf %d", &(dvec[0]), &(dvec[1]), &(dvec[2]), &atype);
      //printf("%lf %lf %lf %d\n", (dvec[0]), (dvec[1]), (dvec[2]), atype);
      fgets(fxyz, 100, fid);
      for (pb_x=-1;pb_x<2;pb_x++)
	for (pb_y=-1;pb_y<2;pb_y++)
	  for (pb_z=-1;pb_z<2;pb_z++)
	    {
	      ivec[0]=dvec[0]+pb_x;
	      ivec[1]=dvec[1]+pb_y;
	      ivec[2]=dvec[2]+pb_z;

	      dot(lats, ivec, rvec);
	      dot(proj_mat, rvec, pvec);
	    
	      if ((isdepth) && ((pvec[2]<zmin) || (pvec[2]>zmax)))
		continue;

	      //printf("%lf %lf %lf\n", rvec[0], rvec[1], rvec[2]);
	      px_nx = round(pvec[0]/px_sz);
	      px_ny = round(pvec[1]/px_sz);
	      //printf("%d %d\n", px_nx, px_ny);
	      for (ext_x=-smear_sz; ext_x<=smear_sz;ext_x++)
		for (ext_y=-smear_sz; ext_y<=smear_sz;ext_y++){
		  paint_x=px_nx+ext_x+px_ct;
		  paint_y=px_ny+ext_y+px_ct;

		  if ( (paint_x<0) || (paint_x>=px_num))
		    continue;
		  if ( (paint_y<0) || (paint_y>=px_num))
		    continue;
		  canvas[paint_x][paint_y] += atype*atype*multivariate(pvec[0],pvec[1],sigma,(px_nx+ext_x)*px_sz, (px_ny+ext_y)*px_sz);
		  //printf("%d %d\n", paint_x, paint_y);
		  //printf("%e %e\n", (pvec[0]-(px_nx+ext_x)*px_sz)/sigma, (pvec[1]-(px_ny+ext_y)*px_sz)/sigma);
		}
	    }
    }
  }
  fclose(fid);

  double maxval=0;
  for (i=0;i<px_num;i++)
    for (j=0;j<px_num;j++)
      maxval = maxval>canvas[i][j]?maxval:canvas[i][j];
  
  fid = fopen("HAADF.out", "w");
  for (i=0;i<px_num;i++){
    for (j=0;j<px_num;j++)
      fprintf(fid, "%lf ", 1-canvas[i][j]/maxval);
    fprintf(fid, "\n");
  }
  fclose(fid);

  for (i=0;i<px_num;i++)
    free(canvas[i]);
  free(canvas);
}
