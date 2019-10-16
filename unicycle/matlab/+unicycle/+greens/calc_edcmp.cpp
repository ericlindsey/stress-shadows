/* edcmp: compute elastic Green's functions for a layered half space.
 *   Sylvain Barbot
 *   sbarbot@ntu.edu.sg
 *   Earth Observatory of Singapore
*/

#include <ctype.h>
#include <string.h>
#include "mex.h"
#include <vector>
#include <string>
#include "calc_edcmp.h"

void init(int n, char * prefix) {

  char grndir[80];

  memcpy(grndir,prefix,n);

  // load Green's functions
  mexPrintf("load Green's functions %s(.ss, .ds, .cl)\n",grndir);
  loadgrn_(grndir);  
}

void layered(int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[]){

  
  // function [ux, uy, uz] = layered_wang(slip,xs,ys,zs,L,W,strike,dip,rake,xr,yr)
  double * x, *y;
  double * ux, *uy, *uz;
  size_t nxx, nxy, ny;
  int nrec,err;

  if (11>nrhs)
     mexErrMsgIdAndTxt("MATLAB:edcmp:minrhs",
                       "not enough input arguments.");

  // test input
  double slip = mxGetScalar(prhs[0]);
  double x1 = mxGetScalar(prhs[1]);
  double x2 = mxGetScalar(prhs[2]);
  double x3 = mxGetScalar(prhs[3]);
  double length = mxGetScalar(prhs[4]);
  double width = mxGetScalar(prhs[5]);
  double strike = mxGetScalar(prhs[6]);
  double dip = mxGetScalar(prhs[7]);
  double rake = mxGetScalar(prhs[8]);

  if (3>nlhs) 
     mexErrMsgIdAndTxt("MATLAB:edcmp:minlhs",
                       "not enough output arguments.");

  /* get the length of each input vector */
  nxx = mxGetM(prhs[9]);
  ny = mxGetN(prhs[9]);
  
  if (1 != ny)
     mexErrMsgIdAndTxt("MATLAB:edcmp:dimension",
                       "x coordinates should be vectors.");

  nxy = mxGetM(prhs[10]);
  ny = mxGetN(prhs[10]);

  if (1 != ny)
     mexErrMsgIdAndTxt("MATLAB:edcmp:dimension",
                       "y coordinates should be vectors.");

  if (nxx != nxy)
     mexErrMsgIdAndTxt("MATLAB:edcmp:dimension",
                       "x and y coordinates should have the same dimension.");

  plhs[0]=mxCreateDoubleMatrix((mwSize)nxx,(mwSize)ny,mxREAL);
  plhs[1]=mxCreateDoubleMatrix((mwSize)nxx,(mwSize)ny,mxREAL);
  plhs[2]=mxCreateDoubleMatrix((mwSize)nxx,(mwSize)ny,mxREAL);

  ux=mxGetPr(plhs[0]);
  uy=mxGetPr(plhs[1]);
  uz=mxGetPr(plhs[2]);

  x=mxGetPr(prhs[9]);
  y=mxGetPr(prhs[10]);

  nrec=(int)nxx;

  // calculate displacements
  edcmp_(&slip,&x1,&x2,&x3,&length,&width,&strike,&dip,&rake,&nrec,y,x,ux,uy,uz,&err);
}

void mexFunction(int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[]) {

  if (nrhs == 0) mexErrMsgTxt("Type 'help hmmvp' for help.");
  
  // clean up command.
  int strlen = mxGetNumberOfElements(prhs[0]) + 1;
  std::vector<char> vfn(strlen);
  mxGetString(prhs[0], &vfn[0], strlen);
  for (int i = 0; i < strlen - 1; i++) vfn[i] = tolower(vfn[i]);
  std::string fn(&vfn[0]);

  // init
  if (fn[0] == 'i' || fn == "init") {
    if (nrhs < 2 || !mxIsChar(prhs[1]))
      mexErrMsgTxt("hmmvp('init',prefix)");
    
    unsigned int n = mxGetNumberOfElements(prhs[1]) + 1;
    char* prefix = new char[n+1];
    mxGetString(prhs[1], prefix, n);
    prefix[n-1]=' ';
    init(n,prefix);
  }
  
  // forward model
  if (fn[0] == 'd' || fn == "displacement") {
     if (12>nrhs)
     mexErrMsgIdAndTxt("MATLAB:edcmp:minrhs",
                       "not enough input arguments.");
    if (3>nlhs) 
     mexErrMsgIdAndTxt("MATLAB:edcmp:minlhs",
                       "not enough output arguments.");
    
    layered(nlhs,plhs,nrhs-1,&(prhs[1]));
  }
};



