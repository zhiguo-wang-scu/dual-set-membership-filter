clear
clc
SM=20;               % target tracking time
rN=20;               % each time has rN Monte Carlo runs
n=4;                 % the dimension of state
m=2;                 % the dimension of the measure
T=1;                 % the sampling time
f=@(x)[x(1)+T*x(3);x(2)+T*x(4);x(3);x(4)];                           % the nonlinear state function
a=420;b=420;
h=@(x)[sqrt((x(1)-a)^2+(x(2)-b)^2);atan2((x(2)-b),(x(1)-a))];        % the nonlinear measurement function
q=10;                     % power
Q=q*[T^3/3 0 T^2/2 0;0 T^3/3 0 T^2/2;T^2/2 0 T 0;0 T^2/2 0 T]; % the bound of process noise
QQ=chol(Q)';               % variance of process noise
uq=[0 0 0 0]';     % mean of process noise
dg=.5;ds=10;
sigma=1;
R=sigma*[ds^2 0;0 (dg/180*pi)^2];       % the variance of measurement noise
ur=[0 0]';                     % mean of measurement noise
error_l=zeros(rN,SM);
trace_l=zeros(rN,SM);
cputime_1=zeros(rN,SM);
for j=1:rN
x0=[50 30 5 5];           % the initial state
x=zeros(n,SM+1);          % the total state
x(:,1)=x0;


X=sample_u1(Q,SM);     % produce process noise
Y=sample_u1(R,SM);     % produce measurement

y=zeros(m,SM+1);                  % the total measure
for i=1:SM
    x(:,i+1)=f(x(:,i))+X(:,i);
    y(:,i+1)=h(x(:,i+1))+Y(:,i);
end
%%%% generate measurement ellipsoid %%%%
xk_m=zeros(4,SM+1);
xk_m(:,1)=x(:,1);
Pk_m=zeros(4,4*(SM+1));
Pk_m(1:4,1:4)=eye(4,4);

%%%% the proposed robust set membership filter %%%%%


xk_l_1=zeros(4,SM+1);
xk_l_1(:,1)=x(:,1)+sample_u1(200*eye(4,4),1);
Pk_l_1=zeros(4,4*(SM+1));
Pk_l_1(1:4,1:4)=200*eye(4,4);

delta_l=zeros(1,SM);
  FF=[1 0 T 0;0 1 0 T;0 0 1 0;0 0 0 1];
for k=1:SM
    tic
    [xk_m([1 2],k+1),Pk_m([1 2],[4*(k+1)-3,4*(k+1)-2])]=vex2(R,[a b],y(:,k+1)); % compute measurement ellipsoid
    [xk_ll_p,Pk_ll_p,~,~]...
        =predicted(FF,FF,xk_l_1(1:4,k),Pk_l_1(1:4,4*k-3:4*k),Q);
    Ep=[1 0 0 0;0 1 0 0];
    
    [xk_l_1(1:4,k+1),Pk_l_1(1:4,4*(k+1)-3:4*(k+1)),delta1]...
        =estimate2(xk_ll_p,Pk_ll_p,Ep,xk_m([1 2],k+1),Pk_m([1 2],[4*(k+1)-3,4*(k+1)-2]));
    cputime_1(j,k)=toc;
    
    delta_l(k)=(xk_l_1(:,k+1)-x(:,k+1))'/Pk_l_1(:,4*(k+1)-3:4*(k+1))*(xk_l_1(1:4,k+1)-x(1:4,k+1));
    error_l(j,k)=(norm(xk_l_1(1,k+1)-x(1,k+1)));
    trace_l(j,k)=trace(Pk_l_1(1:4,4*(k+1)-3:4*(k+1)));
    
end
end

 error1_mean=mean(error_l);
 cuptime1_mean=mean(cputime_1);
 trace1_mean=mean(trace_l);
figure
plot([trace(Pk_l_1(1:4,1:4)),trace1_mean],'-*r')
legend('DSMF')
xlabel('Time (sec)')
ylabel('Trace (m)')
axis([0 SM 450 850])
figure
plot(cuptime1_mean,'-*r')
axis([1 SM -0.5 2])
legend('DSMF')
xlabel('Time (sec)')
ylabel('Running time (sec)')