//This program creates a structure and returns it to MATLAB.
#include "mex.h"
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    //DATA
    mxArray *mydouble,*mystring;
    double *dblptr;
    int i;
    const char *fieldnames[2]; //This will hold field names.
      //PROGRAM
      if((nlhs!=1)||(nrhs!=0))
      {
      mexErrMsgTxt("One output and no input needed");
      return;
      }
      //Create mxArray data structures to hold the data
      //to be assigned for the structure.
      mystring  = mxCreateString("This is my char");
      mydouble  = mxCreateDoubleMatrix(1,100, mxREAL);
      dblptr    = mxGetPr(mydouble);
      for(i = 0; i<100; i++)
          dblptr[i] = (double)i;
      //Assign field names
      fieldnames[0] = (char*)mxMalloc(20);
      fieldnames[1] = (char*)mxMalloc(20);
      memcpy(fieldnames[0],"Doublestuff",sizeof("Doublestuff"));
      memcpy(fieldnames[1],"Charstuff", sizeof("Charstuff"));
      //Allocate memory for the structure
      plhs[0] = mxCreateStructMatrix(1,1,2,fieldnames);
      //Deallocate memory for the fieldnames
      mxFree( fieldnames[0] );
      mxFree( fieldnames[1] );
      //Assign the field values
      mxSetFieldByNumber(plhs[0],0,0, mydouble);
      mxSetFieldByNumber(plhs[0],0,1, mystring);
      //NOTE: mxSetFieldByNumber(..) will automatically take care
      //      of allocating required memory for the fields.
  }//mexFunction(..)
