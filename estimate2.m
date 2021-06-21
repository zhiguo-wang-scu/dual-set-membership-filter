function [x_hat,P_hat,delta]=estimate2(x0,P0,Ep,x_ver,P_ver)
n=length(x0);
a=0.001;
b=0.999;
delta1=0.001;
epsilon=0.001;

[s,phis,k,G,E]=golds1(x_ver,Ep,x0,P0,P_ver,n,a,b,delta1,epsilon);
W=(1-s)*inv(P0)+s*Ep'/P_ver*Ep;
x_hat=inv(W)*((1-s)*inv(P0)*x0+s*Ep'*inv(P_ver)*x_ver);

delta=(x_ver-Ep*x0)'/(Ep*P0*Ep'/(1-s)+P_ver/s)*(x_ver-Ep*x0);

if delta<1
    P_hat=(1-delta)*inv(W);
else
    P_hat=P0;
end
