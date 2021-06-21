function [x_c,P_c]=vex2(R,sensor,y)
N=1000;

 X1=repmat(y,1,N)-sample_u1(R,N);
 X=[X1(1,:).*cos(X1(2,:))+sensor(1)*ones(1,N);X1(1,:).*sin(X1(2,:))+sensor(2)*ones(1,N)];
 [u,R,factor,improv,mxv,mnv,flagstep,lamhist,var,iter,act] =minvol([X;ones(1,N)]);
x_c=X*u;
P_c=(X*diag(u)*X'-x_c*x_c')*2;
P_c=(P_c+P_c')/2;
 
 
 
 