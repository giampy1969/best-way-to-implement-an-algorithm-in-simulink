/* C source file for a 2-state EKF for attitude estimation implemented via legacy code tool  *
 * For more info on the algorithm, its implementation and usage, go to the root folder and   *
 * see the files README.txt, Presentation.pdf, and ekf2sim_example.slx (read the help of     *
 * the ekf2sim block therein). Copyright 2014 The MathWorks, Inc.                            */

#define MAX(a, b)  ((a) > (b)  ? (a) : (b))
#define MIN(a, b)  ((a) < (b)  ? (a) : (b))
#include <math.h>

void ekf_init(double *work1, double *work2, double *p5, double *p6, double *p7) {
    
    unsigned int i;
    
    for (i=0;i<4;i++) work1[0+i]=p5[i];
    for (i=0;i<2;i++) work1[4+i]=p6[i];
    for (i=0;i<3;i++) work1[6+i]=p7[i];
    
    for (i=0;i<33;i++) work2[i]=0;
}

void ekf_out(double *y1, double *u1, double *u2, double *u3, double *work1, double *work2, double *p1, double *p2, double p3 , double p4) {
    
    /* integers                                               */
    unsigned int  i, j, k;
    
    /* scalar variables                                       */
    double  rx, ry, rz, sth, theta, rth, sph, phi, psi, p, q, r;
    double  VPm_det, nrm=0;
    
    /* vectors (states)                                       */
    double              *P  = &work1[0], *x = &work1[4], *vold = &work1[6];
    
    /* vectors (work)                                         */
    double *agps=&work2[0], *y=&work2[3], *xdot=&work2[5], *xm=&work2[7];
    double *A=&work2[9], *Pm=&work2[13], *M=&work2[17], *VPm=&work2[21];
    double *iVPm=&work2[25], *K=&work2[29];
    
    /* calculate outputs */
    y1[0]=work1[4];
    y1[1]=work1[5];
    y1[2]=atan2(u3[1], u3[0]);
    

    /* calculate states */
        
    /* calculate gps acceleration */
    for (i=0;i<3;i++) agps[i]=((u3[i])-(vold[i]))/p4;
    
    /* calculate heading */
    psi=atan2(u3[1], u3[0]);
    
    /* rotate gps acceleration */
    rx=-cos(psi)*agps[0]-sin(psi)*agps[1];
    ry=sin(psi)*agps[0]-cos(psi)*agps[1];
    rz=9.80665-agps[2];
    
    /* solve for theta and phi given imu and gps accelerations */
    sth=(rx*(u2[0])+rz*(sqrt(fabs(rx*rx+rz*rz-(u2[0])*(u2[0])))))/(2.2251e-308+rx*rx+rz*rz);
    theta=atan2(sth*rx-(u2[0]), sth*rz);
    rth=rx*sin(theta)+rz*cos(theta);
    sph=(ry*(u2[1])+rth*(sqrt(fabs(ry*ry+rth*rth-(u2[1])*(u2[1])))))/(2.2251e-308+ry*ry+rth*rth);
    phi=atan2(-sph*ry+(u2[1]), sph*rth);
    
    /* y = [phi;theta]; this works as the "measured" output in the ekf */
    y[0]=phi;y[1]=theta;
    
    /* states (note, phi and theta are re-initialized with the previous estimates) */
    phi=x[0];
    theta=x[1];
    
    /* inputs */
    p=u1[0];
    q=u1[1];
    r=u1[2];
    
    /* xdot */
    xdot[0]=p+sin(phi)*tan(theta)*q+cos(phi)*tan(theta)*r;
    xdot[1]=cos(phi)*q-sin(phi)*r;
    
    /* x propagation */
    xm[0]=x[0]+p4*xdot[0];
    xm[1]=x[1]+p4*xdot[1];
    
    /* A matrix (A = eye(2)+T*a); */
    A[0]=1+(cos(phi)*tan(theta)*q-sin(phi)*tan(theta)*r)*p4;
    A[1]=0+(-sin(phi)*q-cos(phi)*r)*p4;
    A[2]=0+(sin(phi)*(1+tan(theta)*tan(theta))*q+cos(phi)*(1+tan(theta)*tan(theta))*r)*p4;
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
    for (i=0;i<4;i++) Pm[i]+=p1[i];
    
    /* VPm=V+Pm; */
    for (i=0;i<4;i++) VPm[i]=p2[i]+Pm[i];
    
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
    for (i=0;i<4;i++) P[i]*=MIN(nrm, p3)/MAX(nrm, 2.2251e-308);
    
    /* EKF update equations : x=xm+K*(y-xm); */
    x[0]=xm[0]+K[0]*(y[0]-xm[0])+K[2]*(y[1]-xm[1]);
    x[1]=xm[1]+K[1]*(y[0]-xm[0])+K[3]*(y[1]-xm[1]);
    
    for (i=0;i<3;i++) vold[i]=u3[i];
    
}

