%Read parameters into moms_nocell

burn            = mm.burn;          %number of burn in periods
if cf_num == 6 % Caculate network value
    esT         = 60 + 18; % 60 years of usable data with random macro shocks
else 
    esT         = 18 + burn;        % number of ergodic state periods to be simulated
end
scale_f         = mm.scale_f;       %US scale parameter
scale_h         = mm.scale_h;       %Colombia scale parameter
eta             = mm.eta;           %demand elasticity
myalpha         = mm.alpha;         %success function parameter
ah              = mm.ah;            %home Beta parameter
bh              = mm.bh;            %home Beta parameter  
af              = mm.af;            %foreign Beta parameter
bf              = mm.bf;            %foreign Beta parameter
S               = mm.S;             %number of firms
ag              = mm.ag;            % true common (theta0) beta parameter
bg              = mm.bg;            % true common (theta0) beta parameter
theta0          = mm.theta0;        % vector of global firm effects
theta1          = mm.theta2;        % vector of home market firm effects
theta2          = mm.theta2;        % vector of foreign market firm effects
TT              = esT;              % Number of time periods to be simulated
L_p             = mm.L_p;           %arrival rate for jumps in self productivity
L_h             = mm.L_h;
L_f             = mm.L_f;
d               = mm.d;             %exogenous death rate
maxc            = mm.maxc;          %maximum allowable number of clients
delta           = mm.delta;    %exogenous match death hazard
S               = mm.S;        %number of firms
n_size          = mm.n_size;   %number of matches which are learned from
net_size        = mm.net_size; %max number of network effects.
Z               = mm.Z;        %grid for other firm productivity
Phi             = mm.Phi;      %grid for own productivity
X_f             = mm.X_f;      %grid for foreign macro shock
X_h             = mm.X_h;      %grid for home macro shock
actual_h        = mm.actual_h; %actual home macro shock indexes (and yrs)
actual_f        = mm.actual_f; %actual foreign macro shock indexes (and yrs)
L_b             = mm.L_b;
L_z             = mm.L_z;
erg_pz          = mm.erg_pz;    %stationary distribution of buyer productivities
erg_pp          = mm.erg_pp;    %stationary distribution of seller productivites
max_client_prod = mm.max_client_prod; %maximum changes in demand shock over relationship
mult_match_max  = mm.mult_match_max; %maximum number of matches per exogenous state change interval
mms             = mm.mms; %maximum number of matrix rows (memory issues)
