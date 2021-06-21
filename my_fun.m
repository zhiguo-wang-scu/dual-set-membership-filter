function [f_value,delta]=my_fun(rho,x_ver,Ep,x0,P0,P_ver,n)
delta=(x_ver-Ep*x0)'/(Ep*P0/(1-rho)*Ep'+P_ver/rho)*(x_ver-Ep*x0);
f_value=(1-delta)^n*1/det((1-rho)*inv(P0)+rho*Ep'/P_ver*Ep);