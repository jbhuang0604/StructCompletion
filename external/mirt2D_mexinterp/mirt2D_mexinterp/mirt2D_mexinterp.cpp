/* 2D linear interpolation 
   Andriy Myronenko, Feb 2008, email: myron@csee.ogi.edu
   homepage: http://www.bme.ogi.edu/~myron/

ZI = mirt2D_mexinterp(Z,XI,YI) interpolates 2D image Z at the points with coordinates XI,YI.
Z is assumed to be defined at regular spaced points 1:N, 1:M, where [M,N]=size(Z).
If XI,YI values are outside the image boundaries, put NaNs in ZI.

The performance is similar to Matlab's ZI = INTERP2(Z,XI,YI,'linear',NaN).
If Z is a 3D matrix, then iteratively interpolates Z(:,:,1), Z(:,:,2),Z(:,:,3),.. etc.
This works faster than to interpolate each image independaently, such as
 ZI(:,:,1) = INTERP2(Z(:,:,1),XI,YI);
 ZI(:,:,2) = INTERP2(Z(:,:,2),XI,YI);
 ZI(:,:,3) = INTERP2(Z(:,:,3),XI,YI);

 The speed gain is from the precomputation of coefficients for interpolation, which are the same for all images.
 Interpolation of a set of images is useful, e.g. for image registration, when one has to interpolate image and its gradients
 at the same positions.
*/


#include <math.h>
#include "mex.h"


void mirt2D_mexinterp(
double* Z,
double* S,
double* T,
double* F,
int	MN,
int nrows,
int ncols,
int ndim
)
{
    int	n, in1, in2, in3, in4;
    double	t, s;
    double m1, m2, m3, m4, nan;
    int ndx, Fshift, Zshift, i, nrowsncols, fs, ft;
    
    nrowsncols=nrows*ncols;
    nan=mxGetNaN();
    
    
    for (n=0; n < MN; n++) {
        
        
        t=*(T+n);
        s=*(S+n);
        fs=(int)floor(s);
        ft=(int)floor(t);
        
        
        if (fs<1 || s>ncols || ft<1 || t>nrows){
           /* Put nans if outside*/
            for (i = 0; i < ndim; i++) F[n+i*MN]=nan;  }
        else  {
            
            ndx = ft+(fs-1)*nrows;
           
            if (s==ncols){s=s+1; ndx=ndx-nrows; }
            s=s-fs;
            if (t==nrows){t=t+1;ndx=ndx-1; }
            t=t-ft;
            
            in1=ndx-1;
            in2=ndx;
            in4=ndx+nrows;
            in3=in4-1;


			m4=t*s;
            m1=1+m4-t-s;
            m2=t-m4;
            m3=s-m4;
            
            for (i = 0; i < ndim; i++){
                Zshift=i*nrowsncols;
                F[n+i*MN]=Z[in1+Zshift]*m1+Z[in2+Zshift]*m2+Z[in3+Zshift]*m3+Z[in4+Zshift]*m4;
            }
            
        }
    }
    
    
    return;
}
// ------- Main function definitions ----------
/* Input arguments */
#define IN_Z		prhs[0]
#define IN_S		prhs[1]
#define IN_T		prhs[2]

/* Output arguments */
#define OUT_F		plhs[0]

/* Gateway routine */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *Z, *S, *T, *F;
    int  N, M, MN, nrows, ncols, vol, ndim, *newdims;
    const int *dims;
    
     /* Check for input errors */
    if (nlhs>1)
    mexErrMsgTxt("Wrong number of output parameters, usage:  Output_images = mirt2D_mexinterp(Input_images, X, Y)");
    if (nrhs!=3)
    mexErrMsgTxt("Wrong number of input parameters, usage:  Output_images = mirt2D_mexinterp(Input_images, X, Y)");

    if (!mxIsDouble(IN_Z) || !mxIsDouble(IN_S) || !mxIsDouble(IN_T))
    mexErrMsgTxt("mirt2D_mexinterp: Input arguments must be double.");
    
    if ((mxGetNumberOfDimensions(IN_S) != mxGetNumberOfDimensions(IN_T)) ||
      (mxGetNumberOfElements(IN_S) != mxGetNumberOfElements(IN_T))) mexErrMsgTxt("Inputs X, Y must have the same size");
        
    
   /* Get the sizes of each input argument */
    M = mxGetM(IN_S);
    N = mxGetN(IN_S);
    ndim = mxGetNumberOfDimensions(IN_Z);
    dims = mxGetDimensions(IN_Z);
    newdims = (int*) calloc(ndim, sizeof(int));
    
    /*Size of the array to allocate for the interpolated points*/
    *newdims=M;newdims[1]=N; MN=M*N; vol=1;
    if (ndim>2) {
        newdims[2]=dims[2];
        vol=dims[2];}
    
    /*Create the array (2D or 3D) to put the interpolated points*/
    OUT_F = mxCreateNumericArray(ndim, newdims, mxDOUBLE_CLASS, mxREAL);
    
    /* Input image size */
    nrows = *dims;
    ncols = *(dims+1);
    
    
  /* Assign pointers to the input arguments */
    Z      = mxGetPr(IN_Z);
    S       = mxGetPr(IN_S);
    T       = mxGetPr(IN_T);
    
  /* Assign pointers to the output arguments */
    F      = mxGetPr(OUT_F);
    
  /* Do the actual computations in a subroutine */
    mirt2D_mexinterp(Z, S, T, F, MN, nrows, ncols, vol);
    
    free((void*)newdims);
    return;
}


