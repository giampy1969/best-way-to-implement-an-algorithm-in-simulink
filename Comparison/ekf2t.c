/* Level-2 C file S-Function for a 2-state EKF for attitude estimation                       *
 * For more info on the algorithm, its implementation and usage, go to the root folder and   *
 * see the files README.txt, Presentation.pdf, and ekf2sim_example.slx (read the help of     *
 * the ekf2sim block therein). Copyright 2014 The MathWorks, Inc.                            */

#define S_FUNCTION_NAME ekf2t
#define S_FUNCTION_LEVEL 2

#define MAX(a,b)  ((a) > (b)  ? (a) : (b))
#define MIN(a,b)  ((a) < (b)  ? (a) : (b))
#include <math.h>

#include "simstruc.h"

/* mdlCheckParameters, check parameters, this routine is called later from mdlInitializeSizes */
#define MDL_CHECK_PARAMETERS
static void mdlCheckParameters(SimStruct *S)
{
    /* Basic check : All parameters must be real positive vectors                             */
    real_T *pr;                            

    int_T  i, el, nEls;
    for (i = 0; i < 1; i++) {
        if (mxIsEmpty(    ssGetSFcnParam(S,i)) || mxIsSparse(   ssGetSFcnParam(S,i)) ||
            mxIsComplex(  ssGetSFcnParam(S,i)) || !mxIsNumeric( ssGetSFcnParam(S,i))  )
                  { ssSetErrorStatus(S,"Parameters must be real finite vectors"); return; } 
        pr   = mxGetPr(ssGetSFcnParam(S,i));
        nEls = mxGetNumberOfElements(ssGetSFcnParam(S,i));
        for (el = 0; el < nEls; el++) {
            if (!mxIsFinite(pr[el])) 
                  { ssSetErrorStatus(S,"Parameters must be real finite vectors"); return; }
        }
    }

    /* Check number of elements in parameter: W                                               */
    if ( mxGetNumberOfElements(ssGetSFcnParam(S,0)) != 4 )
    { ssSetErrorStatus(S,"W must have 4 elements"); return; }

    /* Check number of elements in parameter: V                                               */
    if ( mxGetNumberOfElements(ssGetSFcnParam(S,1)) != 4 )
    { ssSetErrorStatus(S,"V must have 4 elements"); return; }

     /* Check number of elements in parameter: Plim                                               */
    if ( mxGetNumberOfElements(ssGetSFcnParam(S,2)) != 1 )
    { ssSetErrorStatus(S,"Plim must be a scalar"); return; }

    /* check elements  */
    pr=mxGetPr(ssGetSFcnParam(S,2));
    if ( (pr[0] < 0) )
    { ssSetErrorStatus(S,"Plim cannot be negative"); return; }

     /* Check number of elements in parameter: T                                               */
    if ( mxGetNumberOfElements(ssGetSFcnParam(S,3)) != 1 )
    { ssSetErrorStatus(S,"Sample time must be a scalar"); return; }

    /* check elements  */
    pr=mxGetPr(ssGetSFcnParam(S,3));
    if ( (pr[0] <= 0) )
    { ssSetErrorStatus(S,"T must be positive"); return; }

    /* Check number of elements in parameter: P0                                               */
    if ( mxGetNumberOfElements(ssGetSFcnParam(S,4)) != 4 )
    { ssSetErrorStatus(S,"P0 must have 4 elements"); return; }

    /* Check number of elements in parameter: x0                                               */
    if ( mxGetNumberOfElements(ssGetSFcnParam(S,5)) != 2 )
    { ssSetErrorStatus(S,"x0 must have 2 elements"); return; }

    /* Check number of elements in parameter: ve0                                               */
    if ( mxGetNumberOfElements(ssGetSFcnParam(S,6)) != 3 )
    { ssSetErrorStatus(S,"ve0 must have 3 elements"); return; }

}

/* mdlInitializeSizes - initialize the sizes array ********************************************/
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S,7);                          /* number of expected parameters        */

    /* Check the number of parameters and then calls mdlCheckParameters to see if they are ok */
    if (ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S))
    { mdlCheckParameters(S); if (ssGetErrorStatus(S) != NULL) return; } else return;

    ssSetNumContStates(S,0);                          /* number of continuous states          */
    ssSetNumDiscStates(S,9);                          /* number of discrete states            */

    if (!ssSetNumInputPorts(S,3)) return;             /* number of input ports                */
    ssSetInputPortWidth(S,0,3);                       /* first input port width               */
    ssSetInputPortWidth(S,1,3);                       /* second input port width              */
    ssSetInputPortWidth(S,2,3);                       /* third input port width               */

	ssSetInputPortDirectFeedThrough(S,0,0);           /* first port direct feedthrough flag   */
    ssSetInputPortDirectFeedThrough(S,1,0);           /* second port direct feedthrough flag  */
	ssSetInputPortDirectFeedThrough(S,2,1);           /* third port direct feedthrough flag   */

    if (!ssSetNumOutputPorts(S,1)) return;            /* number of output ports               */
    ssSetOutputPortWidth(S,0,3);					  /* first output port width              */
   
    ssSetNumSampleTimes(S,0);                         /* number of sample times               */

    ssSetNumRWork(S,33);                              /* number real work vector elements     */
    ssSetNumIWork(S,0);                               /* number int_T work vector elements    */
    ssSetNumPWork(S,0);                               /* number ptr work vector elements      */
    ssSetNumModes(S,0);                               /* number mode work vector elements     */
    ssSetNumNonsampledZCs(S,0);                       /* number of nonsampled zero crossing   */
    
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);   /* exception free code (no matlab calls)*/
    
}

/* mdlInitializeSampleTimes - initialize the sample times array *******************************/
static void mdlInitializeSampleTimes(SimStruct *S)
{
    real_T *T = mxGetPr(ssGetSFcnParam(S,3));
    ssSetSampleTime(S, 0, T[0]);
    ssSetOffsetTime(S, 0, 0);
}

/* mdlInitializeConditions - initialize the states ********************************************/
#define MDL_INITIALIZE_CONDITIONS
static void mdlInitializeConditions(SimStruct *S) {
	
    int_T     i;
	real_T  *X0  = ssGetRealDiscStates(S);

	real_T  *P0  = mxGetPr(ssGetSFcnParam(S,4));
	real_T  *x0  = mxGetPr(ssGetSFcnParam(S,5));
	real_T *ve0  = mxGetPr(ssGetSFcnParam(S,6));

	for (i=0;i<4;i++) X0[0+i]=P0[i];
	for (i=0;i<2;i++) X0[4+i]=x0[i];
	for (i=0;i<3;i++) X0[6+i]=ve0[i];
}

/* mdlStart - initialize work vectors *********************************************************/
#define MDL_START
static void mdlStart(SimStruct *S)
{
	int_T   i;
	real_T *w = ssGetRWork(S);

	for (i=0;i<33;i++) w[i]=0;
}

/* mdlUpdate - compute the states ***********************************************************/
#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid)
{

/* block parameters                                       */
real_T              *W  = mxGetPr(ssGetSFcnParam(S,0));
real_T              *V  = mxGetPr(ssGetSFcnParam(S,1));
real_T           *Plim  = mxGetPr(ssGetSFcnParam(S,2));
real_T              *T  = mxGetPr(ssGetSFcnParam(S,3));

real_T              *w  = ssGetRWork(S);


/* states                                                 */
real_T			    *X  = ssGetRealDiscStates(S);
real_T              *P  = &X[0], *x = &X[4], *vold = &X[6];

/* inputs                                                 */
InputRealPtrsType    u  = ssGetInputPortRealSignalPtrs(S,0);
InputRealPtrsType   ab  = ssGetInputPortRealSignalPtrs(S,1);
InputRealPtrsType   ve  = ssGetInputPortRealSignalPtrs(S,2);

/* vectors                                                */
real_T *agps=&w[0], *y=&w[3], *xdot=&w[5], *xm=&w[7];
real_T *A=&w[9], *Pm=&w[13], *M=&w[17], *VPm=&w[21];
real_T *iVPm=&w[25], *K=&w[29];

/* scalar variables                                       */
real_T  rx, ry, rz, sth, theta, rth, sph, phi, psi, p, q, r;
real_T  VPm_det, nrm=0;

/* integers                                               */
int_T                i, j, k;

/* calculate gps acceleration */
for (i=0;i<3;i++) agps[i]=((*ve[i])-(vold[i]))/T[0];

/* calculate heading */
psi=atan2(*ve[1],*ve[0]);

/* rotate gps acceleration */
rx=-cos(psi)*agps[0]-sin(psi)*agps[1];
ry=sin(psi)*agps[0]-cos(psi)*agps[1];
rz=9.80665-agps[2];

/* solve for theta and phi given imu and gps accelerations */
sth=(rx*(*ab[0])+rz*(sqrt(fabs(rx*rx+rz*rz-(*ab[0])*(*ab[0])))))/(2.2251e-308+rx*rx+rz*rz);
theta=atan2(sth*rx-(*ab[0]),sth*rz);
rth=rx*sin(theta)+rz*cos(theta);
sph=(ry*(*ab[1])+rth*(sqrt(fabs(ry*ry+rth*rth-(*ab[1])*(*ab[1])))))/(2.2251e-308+ry*ry+rth*rth);
phi=atan2(-sph*ry+(*ab[1]),sph*rth);

/* y = [phi;theta]; this works as the "measured" output in the ekf */
y[0]=phi;y[1]=theta;

/* states (note, phi and theta are re-initialized with the previous estimates) */
phi=x[0];
theta=x[1];

/* inputs */
p=*u[0];
q=*u[1];
r=*u[2];

/* xdot */
xdot[0]=p+sin(phi)*tan(theta)*q+cos(phi)*tan(theta)*r;
xdot[1]=cos(phi)*q-sin(phi)*r;

/* x propagation */
xm[0]=x[0]+T[0]*xdot[0];
xm[1]=x[1]+T[0]*xdot[1];

/* A matrix (A = eye(2)+T*a); */
A[0]=1+(cos(phi)*tan(theta)*q-sin(phi)*tan(theta)*r)*T[0];
A[1]=0+(-sin(phi)*q-cos(phi)*r)*T[0];
A[2]=0+(sin(phi)*(1+tan(theta)*tan(theta))*q+cos(phi)*(1+tan(theta)*tan(theta))*r)*T[0];
A[3]=1;

/* Pm=W+A*P*A': M=P*A'; */
for(i = 0; i < 2; i++)
   for(j = 0; j < 2; j++)
	   for( M[i+j*2] = 0, k = 0; k < 2; k++) M[i+j*2] += (P[i+k*2])*(A[j+k*2]);

/* Pm=W+A*P*A': Pm=A*M; */
for(i = 0; i < 2; i++)
   for(j = 0; j < 2; j++)
	   for( Pm[i+j*2] = 0, k = 0; k < 2; k++) Pm[i+j*2] += (A[i+k*2])*(M[k+j*2]);

/* Pm=W+A*P*A': Pm=Pm+W; */
for (i=0;i<4;i++) Pm[i]+=W[i];

/* VPm=V+Pm; */
for (i=0;i<4;i++) VPm[i]=V[i]+Pm[i];

/* det(V+Pm) */
VPm_det=VPm[0]*VPm[3]-VPm[2]*VPm[1];

/* inv(V+Pm) */
iVPm[0]= VPm[3]/VPm_det;
iVPm[1]=-VPm[1]/VPm_det;
iVPm[2]=-VPm[2]/VPm_det;
iVPm[3]= VPm[0]/VPm_det;

/* find K: K=Pm*iVPm; */
for(i = 0; i < 2; i++)
   for(j = 0; j < 2; j++)
	   for( K[i+j*2] = 0, k = 0; k < 2; k++) K[i+j*2] += (Pm[i+k*2])*(iVPm[k+j*2]);

/* M=I-K; */
M[0]=1-K[0];
M[1]=0-K[1];
M[2]=0-K[2];
M[3]=1-K[3];

/* P=(eye(2)-K)*Pm; */
for(i = 0; i < 2; i++)
   for(j = 0; j < 2; j++)
	   for( P[i+j*2] = 0, k = 0; k < 2; k++) P[i+j*2] += (M[i+k*2])*(Pm[k+j*2]);

/* ball norm limiter: P=P*min(nrm,Plim)/max(nrm,realmin); */
for (i=0;i<4;i++) nrm+=P[i]*P[i]; nrm=sqrt(nrm);
for (i=0;i<4;i++) P[i]*=MIN(nrm,Plim[0])/MAX(nrm,2.2251e-308);

/* EKF update equations : x=xm+K*(y-xm); */
x[0]=xm[0]+K[0]*(y[0]-xm[0])+K[2]*(y[1]-xm[1]);
x[1]=xm[1]+K[1]*(y[0]-xm[0])+K[3]*(y[1]-xm[1]);

for (i=0;i<3;i++) vold[i]=*ve[i];

}

/* mdlOutputs - compute the outputs ***********************************************************/
static void mdlOutputs(SimStruct *S, int_T tid) {
	
	real_T             *X   = ssGetRealDiscStates(S);
	real_T             *y   = ssGetOutputPortRealSignal(S,0);
    InputRealPtrsType   ve  = ssGetInputPortRealSignalPtrs(S,2);

	/* outputs */
	y[0]=X[4];
	y[1]=X[5];
	y[2]=atan2(*ve[1],*ve[0]);

}
	
/* mdlTerminate - called when the simulation is terminated ***********************************/
static void mdlTerminate(SimStruct *S) {}

/* Trailer information to set everything up for simulink usage *******************************/
#ifdef  MATLAB_MEX_FILE                      /* Is this file being compiled as a MEX-file?   */
#include "simulink.c"                        /* MEX-file interface mechanism                 */
#else
#include "cg_sfun.h"                         /* Code generation registration function        */
#endif
