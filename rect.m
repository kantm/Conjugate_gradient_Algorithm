clear all
close all
clc

% d=G*m;

d=[4;5;3];
G=[2,1;1,2;1,1];

m0=[0;0];    % initial guess

m1=conjgrad(G,d,m0);
% m2=inv(G)*d;

% r= residual = -ve gradient dirn. = down slope dirn. ;
% p=conjugate gradient;

function [m]=conjgrad(G,d,m)
k=0;
dd=d;
r=G'*d;    % Residual = negative gradien dirn. = downslope dirn.
p=r;

rs_old=r'*r;
    for i=1:length(d)
            k=k+1

%         Gp=G*p;
        t=G*p;
        alpha=rs_old/(t'*t); % alpha is chosen such that r(i+1) is orthogonal to r(i) 
        m=m+alpha*p;          % update in x
        dd=dd-alpha*t;
%         r=r-alpha*Gp;         % update in r
        r=G'*dd;
        rs_new=r'*r;
        if(sqrt(rs_new)<1e-10)
            break;
        end
        p=r+(rs_new/rs_old)*p; % (rs_new/rs_old) is chosen such that p(i+1) is conjugate to p(i)
%       next search dirn. (p(i+1)) is built out of current residual & all
%       previous search dirn. (p(i)).
        rs_old=rs_new;
    end
end

% 3 steps of conjugate gradient algo:
%   r(i)=b-A*x(i);
%   alpha(i)=r(i)'*r(i)/(r(i)'*A*r(i));
%   x(i+1)=x(i)+alpha(i)*r(i);
% 3 steps of conjugate gradient algo
