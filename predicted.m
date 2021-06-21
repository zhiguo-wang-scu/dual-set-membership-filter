function [x_hat,P_hat,h_hat,HH]=predicted(f,Jf,x0,P0,Q)
% calculate a predicted ellipsoid
% input: f nonlinear function
%        Jf Jacobian matrix
%        x0 state estimate
%        x1 state estimate of i th model
%        P0 shape matrix
%        Q process noise matrix
% output: x_hat predited state
%         P_hat predicted shape matrix
%         h_hat remainder center
%         HH    remainder shape matrix
n=length(x0);
        x_hat=Jf*x0;  % predited state
        tao(1)=sqrt(trace(Jf*P0*Jf'));
        tao(2)=sqrt(trace(Q));
        sum_tao=sum(tao);
        tao(1)=tao(1)/sum_tao;tao(2)=tao(2)/sum_tao;
        P_hat=Jf*P0*Jf'/tao(1)+Q/tao(2); % predicted shape matrix
        h_hat=zeros(n,1); HH=zeros(n,n);
end
