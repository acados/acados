/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */

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
