function [t,err,u,jac,l] = home_solve(x,bet,a,pi,ex_siz,rh,diag_Q,Q0_d,cost_n,gam,psi)
%function [err,jac] = end_solve(x,bet,a,pi,n,m,ex_siz,rh,diag_Q,Q0_d,cost_n,gam,psi)
%This function evaulates a set of ex_siz by cost_n value functions for the
%search and learning Colombia Project.

c = @(k,m) (1-(1-k).^(1-1/bet))./((1-1/bet)*m.^gam); %cost function

l = zeros(ex_siz*cost_n,1);
l(1:end-ex_siz) = 1-((max(a*(kron(ones(cost_n-1,1),pi) + x(ex_siz+1:end)-x(1:end-ex_siz)).*kron((1:cost_n-1)',ones(ex_siz,1)).^gam,1)).^-bet);
l(end-ex_siz+1:end) = max(1-(a*pi*cost_n^gam).^-bet,0);

denom = 1./(rh+l+kron(ones(cost_n,1),diag_Q)+kron((1:cost_n)'-1,ones(ex_siz,1))*psi);

last_term = zeros(ex_siz*cost_n,1);
last_term(1:end-ex_siz) = l(1:end-ex_siz).*(a*(kron(ones(cost_n-1,1),pi)+x(ex_siz+1:end)-x(1:end-ex_siz)));
last_term(end-ex_siz+1:end) = l(end-ex_siz+1:end).*a.*pi;

err = zeros(ex_siz*cost_n,1);
temp = Q0_d*reshape(x(ex_siz+1:end),ex_siz,cost_n-1);
err(ex_siz+1:end) = denom(ex_siz+1:end).*(-c(l(ex_siz+1:end),kron((2:cost_n)',ones(ex_siz,1)))+temp(:)+psi*kron((1:cost_n-1)',ones(ex_siz,1)).*x(1:end-ex_siz)+last_term(ex_siz+1:end))-x(ex_siz+1:end);
err(1:ex_siz) = denom(1:ex_siz).*(-c(l(1:ex_siz),ones(ex_siz,1))+(Q0_d*x(1:ex_siz))+last_term(1:ex_siz))-x(1:ex_siz);

%use linear indexing of MATLAB to do this part correctly...see help linear

jac = repmat(denom(:)',ex_siz*cost_n,1).*kron(eye(cost_n),Q0_d);
jac(1:ex_siz*cost_n+1:end) = (a*bet*(1-l(:)').*kron(1:cost_n,ones(1,ex_siz)).^gam.*(a*(kron(ones(cost_n,1),pi)'+[x(ex_siz+1:end)',x(end-ex_siz+1:end)']-x(:)')).^-1.*(l(:)>0)'.*denom(:)'.*x(:)'-a*denom(:)'.*(l(:)'+bet*(1-l(:)').*(a*(kron(ones(cost_n,1),pi)'+[x(ex_siz+1:end)',x(end-ex_siz+1:end)']-x(:)')).^-1.*(l(:)>0)'))-1; 
jac(ex_siz+1:ex_siz*cost_n+1:end-ex_siz^2*cost_n) = denom(ex_siz+1:end).*kron((1:cost_n-1)',ones(ex_siz,1))*psi;
jac(ex_siz^2*cost_n+1:ex_siz*cost_n+1:end) = a*bet*(1-l(1:end-ex_siz)').*kron(1:cost_n-1,ones(1,ex_siz)).^gam.*(a*(kron(ones(cost_n-1,1),pi)'+x(ex_siz+1:end)'-x(1:end-ex_siz)').^-1.*(l(1:end-ex_siz)>0)'.*denom(1:end-ex_siz)'.*x(1:end-ex_siz)')+a*denom(1:end-ex_siz)'.*(l(1:end-ex_siz)'+bet*(1-l(1:end-ex_siz)').*kron(1:cost_n-1,ones(1,ex_siz)).^gam.*(a*(kron(ones(cost_n-1,1),pi)'+x(ex_siz+1:end)'-x(1:end-ex_siz)').^-1.*(l(1:end-ex_siz)>0)'));
jac(isnan(jac)) = -1;
%jac(jac == inf) = 0;
%jac(jac == -inf) = 0;

t = 0;
u = zeros(size(x'));
end