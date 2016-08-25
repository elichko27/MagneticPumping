//testSavingFile.cpp
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
  int counter;
  double *ptr_d;
  double *fMatFinalOld; 
  FILE *ptr_fp;

  fMatFinalOld = new double[10];
  
  /* Part A */
  //ptr_d = (double *)malloc(10 * sizeof(double));
  if(!fMatFinalOld)//if(!ptr_d)
    {
      printf("Memory allocation error!\n");
      exit(1);
    }else printf("Memory allocation successful.\n");
  
  /* Part B */
  for(counter = 0; counter < 10; counter++)
    fMatFinalOld[counter] = (double) rand();
  
  /* Part C */
  if((ptr_fp = fopen("test.txt", "wb")) == NULL)
    {
      printf("Unable to open file!\n");
      exit(1);
    }else printf("Opened file successfully for writing.\n");
  
  /* Part D */
  if( fwrite(fMatFinalOld, 10*sizeof(double), 1, ptr_fp) != 1)
    {
      printf("Write error!\n");
      exit(1);
    }else printf("Write was successful.\n");
  fclose(ptr_fp);
  free(fMatFinalOld);
  
  /* Part E */
  ptr_d = (double *)malloc(10 * sizeof(double));
  if(!ptr_d)
    {
      printf("Memory allocation error!\n");
      exit(1);
    }else printf("Memory allocation successful.\n");
  
  /* Part F */
  if((ptr_fp = fopen("test.txt", "rb"))==NULL)
    {
      printf("Unable to open the file!\n");
      exit(1);
    }else printf("Opened file successfully for reading.\n");
  
  /* Part G */
  if(fread(ptr_d, 10 * sizeof( double ), 1, ptr_fp) != 1)
    {
      printf( "Read error!\n" );
      exit( 1 );
    }else printf( "Read was successful.\n" );
  fclose(ptr_fp);
  
  /* Part H */
  printf("The numbers read from file are:\n");
  for(counter = 0; counter < 10; counter++)
    printf("%f ", ptr_d[counter]);
  
  /* Part I */
  free(ptr_d);
  return 0;
}
