function [V,l_opt] = val_loop_h(Q0,Q0_d,a,pi,mm)
%This function solves for the value function and optimal search intensity 

pi_orig = pi;

bet         = mm.b; %cost parameter
rh          = mm.r; %discount factor
dim0        = mm.dim0; %size of th_g common
dim1        = mm.dim1; %size of th_h home
net_size      = mm.net_size; %max number of network effects
gam         = mm.gam; %network effect parameter
tol         = mm.v_tolerance; %tolerance for value function loop
cscale      = mm.cscale;  %counterfactual scale, used when shutting off network effect

V = zeros(size(Q0,1),dim0,dim1,net_size+1); %[macro-prod states, theta_g,theta_h, network]
scl = max(pi_orig);
pi = pi_orig/scl;

c = @(x,n) ((1+x).^(1+1/bet)-1)/((1+1/bet)*n^gam*scl)/counter_scale; %cost function
%c = @(l,n) (((1+l)/n^gam).^(1+1/bet)-1)/((1+1/bet)*scl);


diag_Q = abs(diag(Q0));
l_opt = zeros(size(V)); %optimal search intensities

options = optimset('Display','off','Jacobian','on','GradObj','on','GradConstr','on','TolX',tol/100/scl,'TolFun',tol/scl,'DerivativeCheck','off','maxiter',10000,'Algorithm','trust-region-reflective');
tic;
    %first step is linear
    for j = 1:dim0
        for k = 1:dim1
            %display([j,k]);
            %first backwards induction step
        	l_opt(:,j,k,net_size+1) = max(((net_size+1)^gam*a(j,k)*pi*scl).^bet-1,0);
            %den = rh+l_opt(:,j,k,net_size+1)+diag_Q;
            den = rh+diag_Q;
            RHS = -c(l_opt(:,j,k,net_size+1),net_size+1)+a(j,k)*l_opt(:,j,k,net_size+1).*pi;
            Q0_diag = -Q0_d;
            %Q0_diag(1:size(Q0_d,1)+1:end) = den-l_opt(:,j,k,net_size+1);
            Q0_diag(1:size(Q0_d,1)+1:end) = den;
            V(:,j,k,net_size+1) = max(Q0_diag^-1*RHS,0);     
            
            flag = 0;
            for l = 1:net_size
                init = V(:,j,k,net_size-l+2); %use the previous iterations solution as initial guess 
                [V(:,j,k,net_size-l+1),~,flag_update] = fsolve(@(x) sim_solve_h(x,bet,a(j,k),pi,net_size+1-l,size(Q0,1),rh,diag_Q,Q0_d,V(:,j,k,net_size-l+2),gam,scl,cscale),init,options);
                if flag_update<1 %check for convergence
                    flag = flag + 1;
                end
                [~,~,l_opt(:,j,k,net_size+1-l)] = sim_solve_h(V(:,j,k,net_size+1-l),bet,a(j,k),pi,net_size+1-l,size(Q0,1),rh,diag_Q,Q0_d,V(:,j,k,net_size-l+2),gam,scl,cscale);
            end
             if flag > 0
                 display('WARNING: home value function convergence issue!');
                 %punishment = punishment + flag/10; %punishment for not converging
             end   
        end
    end
end
