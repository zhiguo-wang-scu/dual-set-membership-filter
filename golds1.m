function [s,phis,k,G,E]=golds1(x_ver,Ep,x0,P0,P_ver,n,a,b,delta,epsilon)

t=(sqrt(5)-1)/2; h=b-a;
phia=my_fun(a,x_ver,Ep,x0,P0,P_ver,n);
phib=my_fun(b,x_ver,Ep,x0,P0,P_ver,n);
p=a+(1-t)*h; q=a+t*h;
phip=my_fun(p,x_ver,Ep,x0,P0,P_ver,n);
phiq=my_fun(q,x_ver,Ep,x0,P0,P_ver,n);
k=1; G(k,:)=[a, p, q, b];
max_k=100;

while(abs(phib-phia)>epsilon)||(h>delta)
if(phip<phiq)
b=q; phib=phiq; q=p; phiq=phip;
h=b-a; p=a+(1-t)*h; phip=my_fun(p,x_ver,Ep,x0,P0,P_ver,n);
else
a=p; phia=phip; p=q; phip=phiq;
h=b-a; q=a+t*h; phiq=my_fun(q,x_ver,Ep,x0,P0,P_ver,n);
end
k=k+1; G(k,:)=[a, p, q, b];
if k>max_k
    break;
end
end
ds=abs(b-a); dphi=abs(phib-phia);
if(phip<=phiq)
s=p; phis=phip;
else
s=q; phis=phiq;
end
E=[ds,dphi];