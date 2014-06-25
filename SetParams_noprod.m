% script SetParams
%This script is a handy way to deal with all the parameters of the model.
%The settings here are passed into the main program.

mm = struct();

%% technology parameters

mm.r         = 0.05;       % Rate of time preference
mm.d         = 0.03;       % Component of time preference due to exogenous death
mm.delta     = delta;      % Exogenous match separation rate
mm.b         = beta;       % Cost function parameter
mm.scale_f   = scale_f;    % Export profit function scale parameter
mm.scale_h   = scale_h;    % Domestic profit function scale parameter
mm.eta       = 5;          % Demand elasticity 
mm.gam       = gam;        % Network effect parameter
mm.cs        = cs;         % Cost scaling parameter

%% Theta distributions

mm.ah        = ah;         % Beta function, home (theta1) success parameter
mm.bh        = bh;         % Beta function, home (theta1) failure parameter
mm.af        = mm.ah;         % Beta function, foreign (theta2) success parameter (assume same as home)
mm.bf        = mm.bh;         % Beta function, foreign (theta2) success parameter (assume same as home)
mm.ag        = ag;         % Beta function, theta0 success parameter
mm.bg        = bg;        % Beta function, theta0 failure parameter
mm.F         = exp(lnF);   % cost of maintaining a client--home or foreign
mm.alpha     = alp;         %weight of "common" theta in determining match probabilities

%% Shock processes

mm.mean_phi    = 0;        % Mean productivity level
mm.rho_phi     = 0.74;     % root, idiosyncratic productivity shock (JC's estimate)
mm.sig_eps_phi = 0.46;     % std. deviation of innovation in productivity shock (JC's estimate)

mm.mean_z      = 0;        % Mean product appeal in U.S.
%mm.sig_eps_z   = sig_eps_z;% std. deviation of innovation in product appeal shock
%mm.rho_z       = rho_z;    % root, product appeal shock 

mm.mean_x_h    = 0;        % Mean home macro state 
mm.rho_x_h     = 0.961;    % root, home macro shock (JC's results, final, intermediate and cap. good expend.)
mm.sig_eps_x_h = 0.081;    % std. deviation of innovation in home macro shock (JC's results, final, intermediate and cap. good expend.)

mm.mean_x_f    = 0;        % Mean foreign macro state
mm.rho_x_f     = 0.953;    % root, foreign macro shock (Van's results, lns_US regression)
mm.sig_eps_x_f = 0.052;    % std. deviation of innovation in U.S. macro shock (Van's results, lns_US regression)

%% Discretization of state-space
mm.grid_length   = 2.5;    % number of standard deviations from mean used for discretization
mm.n_size        = 20;     % Maximum number of informative signals per firm (WAS 20)
mm.net_size      = 30;     %maximum number of network effects
mm.z_size        = 7;     % Number of discretized demand shock states (2*n+1) 
mm.phi_size      = 7;     % number of different discretized profit shocks (2*n+1)
mm.x_size        = 7;     % Number of different discretized macro shocks; same for home and foreign (2*n+1)
mm.lambda_size   = 25;     % Number of possible effort levels (WAS 40)
mm.theta_size    = 20;     % Number of possible market potential values (WAS 51)
mm.dim0          = 3;      % Number of possible theta0 values (common to both markets)
mm.dim1          = 7;      % Number of possible theta1 values (specific to home market);
mm.dim2          = 7;      % Number pf possible theta2 values (specific to foreign market);

mm.theta0           = 1/mm.dim0:1/mm.dim0:1;
mm.theta1           = 1/mm.dim1:1/mm.dim1:1;
mm.theta2           = 1/mm.dim2:1/mm.dim2:1;
mm.theta0(mm.dim0)  =  mm.theta0(mm.dim0) - 0.0001;
mm.theta1(mm.dim1)  =  mm.theta1(mm.dim1) - 0.0001;
mm.theta2(mm.dim2)  =  mm.theta2(mm.dim2) - 0.0001;

%% Solution parameters
mm.v_tolerance   = 1e-3;  % convergence tolerance, value function iterations (WAS .005)
mm.pi_tolerance  = 1e-5;  % convergence tolerance, profit function (WAS .001)
mm.T             = 50;     % horizon for calculating profit function
mm.S             = 10000;      % number of potential exporting firms to simulate (WAS 2000)
mm.burn          = 15;        %number of burn-in periods

%% Simulation restrictions
mm.maxc            = mm.net_size; %maximum number of current clients (follows old program)
mm.max_client_prod = 5000; %maximum changes in demand shock over relationship
mm.mult_match_max  = 5000; %maximum number of matches per exogenous state change interval
mm.mms             = 50000; %max event number (max matrix size)

%% Exogenous Jump Process Parameters

% Uncomment commands below to reestimate exogenous variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--Get exogenous parameters. Once these are estimated, no need to       --%
%  estimate again                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  DATA I HAVE
%     col_exp_raw = dlmread('Col_Ind_Gds_Exp_for_MLE_trial_2011_7_8.txt','\t');
%     col_exp = detrend(log(col_exp_raw(:,3))); %detrend
%     us_exp_raw = dlmread('US_man_exp_2011_8_10.txt','\t');
%     us_exp = detrend(log(us_exp_raw(:))); %detrend
% 
%     [~,sig_h,gam_h] = OU_Calibrate_ML(col_exp,1); %arguments are data and time interval, output is estimated mean, volatility, and mean reversion
%     [~,sig_f,gam_f] = OU_Calibrate_ML(us_exp,1); %arguments are data and time interval, output is estimated mean, volatility, and mean reversion
%
%     L_h = gam_h * mm.x_size; %lambda, arrival rate of shock
%     D_h = sig_h*L_h^(-.5); %delta, size of jump states
%     [Q_h,X_h] = makeq(L_h,D_h,mm.x_size);
%
%     L_f = gam_f * mm.x_size; %lambda, arrival rate of shock
%     D_f = sig_f*L_f^(-.5); %delta, size of jump states
%     [Q_f,X_f] = makeq(L_f,D_f,mm.x_size);
%
%     actual_h = zeros(size(col_exp,1),2);
%     actual_f = zeros(size(us_exp,1),2);
%     actual_h(:,1) = 1990:2007;
%     actual_f(:,1) = 1990:2009;
%     for k = 1:size(actual_h,1)
%         [~,actual_h(k,2)] = min(abs(col_exp(k)-X_h)); 
%     end
%     for k = 1:size(actual_f,1)
%         [~,actual_f(k,2)] = min(abs(us_exp(k)-X_f)); 
%     end
%
%     save('exog_est/exog.mat','L_h','D_h','Q_h','X_h','L_f','D_f','Q_f','X_f','actual_f','actual_h');

load exog_est/exog.mat
% create intensity matrices for productivity parameters
%NOTE: I choose .045 for D_z and D_p because it delivers a maximum
%difference in shipment scale of about 1500.
%2*(.045*15+4*15*.045+.3) = 7.35 (approx log(1500)) 
%D_z = .045; %delta, size of jump states
[Q_z,Z] = makeq(L_z,D_z,mm.z_size);
erg_pz = make_erg(L_z,D_z,Z); 

%D_p = .035; %delta, size of jump states
%[Q_p,Phi] = makeq(L_p,D_p,mm.phi_size);
%erg_pp = make_erg(L_p,D_p,Phi);
% normal around zero (lognormal ultimately)
erg_pp = zeros(2 * mm.phi_size,1);
for k = 1:2 * mm.phi_size + 1
    erg_pp(k) = normpdf(-3 + 3/mm.phi_size * (k-1));
end
erg_pp = erg_pp./sum(erg_pp);
Phi = (-3:3/mm.phi_size:3)' * sig_p;
Q_p = zeros(2 * mm.phi_size + 1);

% get the home and foreign aggregate intensity matrices

[Q0_h,Q0_h_d,st_h] = makebigq(Q_p,Q_h,mm.phi_size,mm.x_size,Phi,X_h);
[Q0_f,Q0_f_d,st_f] = makebigq(Q_p,Q_f,mm.phi_size,mm.x_size,Phi,X_f);

%create Q_z with zeros on the diagonal (will use this later)
Q_z_d = Q_z;
j = size(Q_z,1);
Q_z_d(1:(j+1):end) = 0;

mm.actual_h     = actual_h; %actual indexes of home macro shocks for available years
mm.actual_f     = actual_f; %actual indexes of foreign macro shocks for available years    
mm.L_b          = L_b;      %shipment hazard   

mm.Z            = Z;        %buyer productivities
mm.Phi          = Phi;      %seller productivies
mm.X_f          = X_f;      %foreign macro shocks
mm.X_h          = X_h;      %home macro shocks
mm.st_h         = st_h;     %home state (for use with intensity matrix)
mm.st_f         = st_f;     %foreign state (for use with intensity matrix)

mm.L_h          = L_h;      %arrival rate for jumps in home macro shock
mm.D_h          = D_h;      %size of jump in home macro shock
mm.Q_h          = Q_h;      %intensity matrix for home macro shock

mm.L_f          = L_f;      %arrival rate for jumps in foreign macro shock
mm.D_f          = D_f;      %size of jump in foreign macro shock
mm.Q_f          = Q_f;      %intensity matrix for foreign macro shock

mm.L_p          = L_p;      %arrival rate for jumps in own productivity
mm.D_p          = D_p;      %size of jump in own productivity
mm.Q_p          = Q_p;      %intensity matrix for own productivity
mm.erg_pp       = erg_pp;   %ergodic distribution of seller productivities

mm.L_z          = L_z;      %arrival rate for jumps in other firms productivity
mm.D_z          = D_z;       %size of jump in other firms productivity
mm.Q_z          = Q_z;      %intensity matrix for demand shocks 
mm.Q_z_d        = Q_z_d;    %with zeros on the diagonal
mm.erg_pz       = erg_pz;   %ergodic distribution of demand shocks

mm.Q0_h         = Q0_h;     %intensity matrix for home state
mm.Q0_f         = Q0_f;     %intensity matrix for foreign state
mm.Q0_h_d       = Q0_h_d;   %home with zeros on diagonal
mm.Q0_f_d       = Q0_f_d;   %foreign with zeros on diagonal
