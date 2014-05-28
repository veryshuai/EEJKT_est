function [V,l_opt,punishment] = val_loop_h(Q0,Q0_d,a,pi,mm)
%This function solves for the value function and optimal search intensity 

    % fix original unscaled lifetime relationship profit expectation
    pi_orig = pi;
    
    % Read in parameters
    bet         = mm.b; %cost parameter
    rh          = mm.r; %discount factor
    dim0        = mm.dim0; %size of th_g common
    dim1        = mm.dim1; %size of th_h home
    net_size      = mm.net_size; %max number of network effects
    gam         = mm.gam; %network effect parameter
    tol         = mm.v_tolerance; %tolerance for value function loop
    cscale      = mm.cs; %counterfactual scaling parameter
    
    % Initialize value function and policy function
    V = zeros(size(Q0,1),dim0,dim1,net_size+1); %[macro-prod states, theta_g,theta_h, network]
    l_opt = zeros(size(V)); %optimal search intensities
    
    % Scale the profit expectation 
    scl = max(pi_orig);
    pi = pi_orig/scl;
    
    % Scaled cost function
    c = @(x,n) ((1+x).^(1+1/bet)-1)/((1+1/bet)*n^gam*scl/cscale);
    
    % Initialize punishment for non-convergence
    punishment = 0; 
    
    % Make sure that the Q matrix (intensity matrix for events) is positive 
    diag_Q = abs(diag(Q0));
    
    % Options for solver
    options = optimset('Display','off','Jacobian','on','GradObj','on','GradConstr','on','TolX',tol/100/scl,'TolFun',tol/scl,'DerivativeCheck','off','maxiter',10000,'Algorithm','trust-region-reflective');

    % Start timer
    tic;

    % Backward induction to solve for ;value
    for j = 1:dim0
        for k = 1:dim1
            % first backwards induction 
        	l_opt(:,j,k,net_size+1) = max(((net_size+1)^gam*a(j,k)*pi*scl/cscale).^bet-1,0);
            den = rh+diag_Q;
            RHS = -c(l_opt(:,j,k,net_size+1),net_size+1)+a(j,k)*l_opt(:,j,k,net_size+1).*pi;
            Q0_diag = -Q0_d;
            Q0_diag(1:size(Q0_d,1)+1:end) = den;
            V(:,j,k,net_size+1) = max(Q0_diag^-1*RHS,0);     
            
            % subsequent steps 
            flag = 0;
            for l = 1:net_size

                % use the previous iterations solution as initial guess 
                init = V(:,j,k,net_size-l+2); 

                % solve for value using fsolve
                [V(:,j,k,net_size-l+1),~,flag_update] = fsolve(@(x) sim_solve_h(x,bet,a(j,k),pi,net_size+1-l,size(Q0,1),rh,diag_Q,Q0_d,V(:,j,k,net_size-l+2),gam,scl,cscale),init,options);

                % check for convergence
                if flag_update<1
                    flag = flag + 1;
                end

                % get the policy function
                [~,~,l_opt(:,j,k,net_size+1-l)] = sim_solve_h(V(:,j,k,net_size+1-l),bet,a(j,k),pi,net_size+1-l,size(Q0,1),rh,diag_Q,Q0_d,V(:,j,k,net_size-l+2),gam,scl,cscale);
            end

            % check for convergence problems 
            if flag > 0
                display('WARNING: home value function convergence issue!');
                punishment = punishment + flag/10; %punishment for not converging
            end
        end
    end
end
