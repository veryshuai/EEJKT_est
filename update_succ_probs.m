function [mu_f, mu_h] = update_success_probs(succ_params)
% This function draws success probability draws for a firm at home and abroad, in accordance with parameterized beta distribution

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate acceptance rates for home and foreign market

    %read_parameters
    theta0              = succ_params{1};
    theta1              = succ_params{2};
    theta2              = succ_params{3};
    bp                  = succ_params{4};
    myalpha             = succ_params{5};
    ag= bp(1); bg= bp(2); af= bp(3); bf= bp(4); ah= bp(5); bh= bp(6);
    
    %warnings off (occassional betainv did not converge warnings) 
    warning off 

    %Get random draws of theta
    th0 = betainv(rand,ag,bg);
    th1 = betainv(rand,ah,bh);
    th2 = betainv(rand,af,bf);
    th_draw   = cat(2,th0,th1,th2);
    
    %warnings on
    warning on
    
    % Read in the theta grid values 
    th    = cell(3,1);
    th{1} = theta0;
    th{2} = theta1;
    th{3} = theta2;
    
    % Find the closest grid values to each firms random draw
    indx1 = zeros(1,3);
    theta = zeros(1,3);
    for j = 1:3  % Map theta0, theta1, and theta2 draws onto grid.  
        dif            = abs(th_draw(1,j)-th{j}); 
        [~,indx1(1,j)] = min(dif);
        theta(1,j)     = th{j}(indx1(1,j));
    end
    
    % Get common vs independent components in success probabilities
    mu_h = myalpha*theta(1,1)+(1-myalpha)*theta(1,2); %true home success probability
    mu_f = myalpha*theta(1,1)+(1-myalpha)*theta(1,3); %true foreign success probability

end
