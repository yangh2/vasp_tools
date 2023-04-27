#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double pi=3.1415926;

/* norm distribution in 2D */
double Gaussian(double cx, double sigma, double x){
  double c = sqrt((double)1.0/2/pi);
  double f1=exp(-0.5*(x-cx)*(x-cx)/sigma/sigma)/sigma;
  return c*f1;
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

int vvdot(double vec[3], double vec1[3], double *prod){
  int j;
  *prod=0;
  for (j=0;j<3;j++)
    *prod+=vec[j]*vec1[j];
  return 0;
}

int main(){
  int i,j;
  int bin_num, nsamples, isdepth, na, ntype, *atom_sp, atom_sp_Hash[1000];
  double bin_sz, zmin, zmax;
  double proj_vec[3], lats[3][3], sigma;
   
  FILE *fid=fopen("zdist.in", "r");
  fscanf(fid, "%d", &ntype);
  atom_sp = (int *) malloc(sizeof(int)*ntype);
  for (i=0;i<ntype;i++){
    fscanf(fid, "%d", &(atom_sp[i]));
    atom_sp_Hash[atom_sp[i]]=i+1;
  }
  fscanf(fid, "%d", &bin_num);
  fscanf(fid, "%lf", &bin_sz);
  fscanf(fid, "%lf", &sigma);
  for (i=0;i<3;i++)
    fscanf(fid, "%lf", &(proj_vec[i]));
  zmin=0;
  for (j=0;j<3;j++)
    zmin+=proj_vec[j]*proj_vec[j];
  for (j=0;j<3;j++)
    proj_vec[j]/=sqrt(zmin);
  //printf("%lf \n" , zmin);
  /* for (j=0;j<3;j++) */
  /*   printf("%lf ", (proj_vec[j])); */
  /* printf("\n"); */

  fscanf(fid, "%lf %lf", &zmin, &zmax);
  if (zmin>zmax) isdepth=0; else isdepth=1;
  fscanf(fid, "%d", &nsamples);

  int smear_sz=(int)(sigma/bin_sz*3);

  double **canvas;
  canvas = (double **) malloc(sizeof(double*)*bin_num);
  for (i=0;i<bin_num;i++){
    canvas[i] = (double *) malloc(sizeof(double)*(ntype+1));
    for (j=0;j<ntype+1;j++)
      canvas[i][j]=0;
  }
  
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
  double dvec[3], ivec[3], rvec[3], pz;
  int pb_x, pb_y, pb_z, px_nz, ext_z, paint_z;
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
	      vvdot(proj_vec, rvec, &pz);
	    
	      if ((isdepth) && ((pz<zmin) || (pz>zmax)))
		continue;

	      //printf("%lf %lf %lf\n", rvec[0], rvec[1], rvec[2]);
	      px_nz = round((pz-zmin)/bin_sz);
	      //printf("%d %d\n", px_nx, px_ny);
	      for (ext_z=-smear_sz; ext_z<=smear_sz;ext_z++){
		paint_z=px_nz+ext_z;
		
		if ( (paint_z<0) || (paint_z>=bin_num))
		  continue;

		canvas[paint_z][0] += Gaussian(pz,sigma,(px_nz+ext_z)*bin_sz+zmin);
		canvas[paint_z][atom_sp_Hash[atype]] += Gaussian(pz,sigma,(px_nz+ext_z)*bin_sz+zmin);
		}
	    }
    }
  }
  fclose(fid);

  double sum=0;
  for (i=0;i<bin_num;i++)
    sum += canvas[i][0];
  
  fid = fopen("zdist.dat", "w");
  for (i=0;i<bin_num;i++){
    fprintf(fid, "%.8lf ", zmin+i*bin_sz);
    for (j=0;j<ntype+1;j++){
      fprintf(fid, "%.8lf ", canvas[i][j]/sum/bin_sz);
    }
    fprintf(fid, "\n");
  }
  fclose(fid);

  free(canvas);
}
