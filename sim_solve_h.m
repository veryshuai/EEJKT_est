%function [nada1,err,nada2,jac,l] = sim_solve(x,bet,a,pi,n,j,ex_siz,rh,diag_Q,Q0_d,V_fail,V_succ,gam,scl)
function [err,jac,l] = sim_solve_h(x,bet,a,pi,n,ex_siz,rh,diag_Q,Q0_d,V_succ,gam,scl,cscale)
%This function creates and solves a set of ex_siz value functions for the
%search and learning Colombia Project.

nada1 = [];
nada2 = [];

c = @(x,n) ((1+x).^(1+1/bet)-1)/((1+1/bet)*n^gam*scl/cscale); %cost function

l = max((max(a*(pi + V_succ-x),0)*n^gam*scl/cscale).^bet-1,0);

denom = rh+l+diag_Q;
last_term = l.*(a*(pi+V_succ)+(1-a)*x);

err = denom.^-1.*(-c(l,n)+Q0_d*x+last_term)-x;

dl = (-a*bet*n^gam*(a*(pi+V_succ-x)*n^gam*scl/cscale).^(bet-1)).*(l>0);
dc = (1+l).^(1/bet)/(n^gam*scl/cscale) .* dl;
dlast = (a*(pi + V_succ)+(1-a)*x).*dl+(1-a)*l;
dden = denom.^-2.*dl;
derr = dden.*(-c(l,n)+Q0_d*x+last_term)+denom.^-1.*(-dc+dlast)-1;
jac = bsxfun(@times,denom.^-1,Q0_d)';
jac(1:ex_siz+1:end) = derr'; 
jac = jac';

jac(isnan(jac)==1 | abs(jac) == inf) = 0; %filter out weird stuffs

%ALTERNATIVE COST FUNCTION
% c = @(l,n) (((1+l)/n^gam).^(1+1/bet)-1)/(scl*(1+1/bet));
% 
% l = max(n^gam*(max(a*(pi + V_succ-x),0)*scl).^bet-1,0);
% 
% denom = rh+l+diag_Q;
% last_term = l.*(a*(pi+V_succ)+(1-a)*x);
% 
% err = denom.^-1.*(-c(l,n)+Q0_d*x+last_term)-x;
% 
% dl = (-a*bet*n^gam*(a*(pi+V_succ-x)*scl).^(bet-1)).*(l>0); %should be an extra scale here by my algebra, but then I fail derivative check...
% dc = ((1+l)/n^gam).^(1/bet)/(scl*n^gam) .* dl;
% dlast = (a*(pi + V_succ)+(1-a)*x).*dl+(1-a)*l;
% dden = denom.^-2.*dl;
% derr = dden.*(-c(l,n)+Q0_d*x+last_term)+denom.^-1.*(-dc+dlast)-1;
% jac = bsxfun(@times,denom.^-1,Q0_d)';
% jac(1:ex_siz+1:end) = derr'; 
% jac = jac';
% 
% jac(isnan(jac)==1 | abs(jac) == inf) = 0; %filter out weird stuffs

end