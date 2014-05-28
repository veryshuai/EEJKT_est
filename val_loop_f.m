function [V,l_opt,punishment] = val_loop_f(Q0,Q0_d,a,pi,mm)
%calculates value function and optimal search intensity for foreign

    % fix original unscaled lifetime relationship profit expectation
    pi_orig = pi; %for help debugging.
    
    % read in parameters
    bet         = mm.b; %cost function parameter
    rh          = mm.r; %discount parameter
    dim0        = mm.dim0; %size of general th_g
    N           = mm.n_size; %maximum learning matches
    net         = mm.net_size; %maximum network effects
    tol         = mm.v_tolerance; %tolerance for loop below
    gam         = mm.gam; %network effect parameter
    cscale      = mm.cs; %counterfactual scaling term
    
    % initialize punishment for non-convergence
    punishment = 0; 
    
    % Inititialize value function and policy functions
    V = zeros(size(Q0,1),N+1,N+1,dim0,net+1); %[macro state,trials,successes,general type]
    l_opt = zeros(size(V)); %optimal search intensities
    
    %scaling
    scl = max(pi_orig);
    
    %scaled cost function
    c = @(x,n) ((1+x).^(1+1/bet)-1)/((1+1/bet)*n^gam*scl/cscale);
    
    %scaled expectation of profit from a relationship
    pi = pi_orig/scl;

    % make sure intensity matrix diagonals are positive
    diag_Q = abs(diag(Q0));
    
    % solver options
    options = optimset('Display','off','Jacobian','on','GradObj','on',...
            'GradConstr','on','TolX',tol/100/scl,'TolFun',tol/scl,...
            'DerivativeCheck','off','maxiter',10000,'Algorithm','trust-region-dogleg');

    % backward induction to solve value function
    for k = 1:dim0  
        flag= 0;
        %Solve for value function, last learning step.  j-1 is successes, k is
        for j = 1:N+1
            % solve for value when network effect maxed out
            l_opt(:,N+1,j,k,net+1) = max(((net+1)^gam*a(N+1,j,k)*pi*scl/cscale).^bet-1,0);
            Q0_diag = -Q0_d;
            den = rh+l_opt(:,N+1,j,k,net+1)+diag_Q;
            last = a(N+1,j,k)*l_opt(:,N+1,j,k,net+1).*pi-c(l_opt(:,N+1,j,k,net+1),net+1);
            Q0_diag(1:size(Q0_d,1)+1:end) = den-l_opt(:,N+1,j,k,net+1);
            V(:,N+1,j,k,net+1) = max(Q0_diag^-1*last,0);
     
            %learning maxed out
            for m = 1:net+1-j

                % initial guess
                init = V(:,N+1,j,k,net+2-m);

                % solve for value function
                [V(:,N+1,j,k,net+1-m),~,flag_update] = fsolve(@(x) sim_solve_h(x,bet,a(N+1,j,k),pi,net+1-m,size(Q0,1),rh,diag_Q,Q0_d,V(:,N+1,j,k,net+2-m),gam,scl,cscale),init,options);

                %check for convergence
                if flag_update<1 
                     flag = flag + 1;
                end

                % get policy function
                [~,~,l_opt(:,N+1,j,k,net+1-m)] = sim_solve_h(V(:,N+1,j,k,net+1-m),bet,a(N+1,j,k),pi,net+1-m,size(Q0,1),rh,diag_Q,Q0_d,V(:,N+1,j,k,net+2-m),gam,scl,cscale);
            end
        end
        
        % the rest of the value surface, 
        % N+1-m is current trial number, j is successes.
        for m = 1:N
            for j = 1:N+1-m

                %use the previous iterations solution as initial guess
                init = V(:,N+2-m,j,k,j); 

                % solve for value;
                [V(:,N+1-m,j,k,j),~,flag_update] = fsolve(@(x) sim_solve(x,bet,a(:,:,k),pi,N+1-m,j,size(Q0,1),rh,diag_Q,Q0_d,V(:,N-m+2,j,k,j),V(:,N-m+2,j+1,k,j+1),gam,scl,cscale),init,options);

                %check for convergence
                if flag_update<1 
                    flag = flag + 1;
                end

                % policy function
                [~,~,l_opt(:,N+1-m,j,k,j)] = sim_solve(V(:,N+1-m,j,k,j),bet,a(:,:,k),pi,N+1-m,j,size(Q0,1),rh,diag_Q,Q0_d,V(:,N-m+2,j,k,j),V(:,N-m+2,j+1,k,j+1),gam,scl,cscale);

            end
        end

        % add to punishment for convergence problems
        if flag > 0
            display('WARNING: foreign value function convergence issue!');
            punishment = punishment + flag/10; %punishment for not converging
        end
    end

    % end timer
    toc

end %end function

