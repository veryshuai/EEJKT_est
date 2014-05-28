function [pi,pi_z,c_val] = makepie(sf,st,Z,Q0,Q0_d,Q_z,Q_z_d,erg_pz,mm)

% READ IN PARAMETERS
rh              = mm.r;     %time pref
del             = mm.delta; %exogenous death rate
F               = mm.F;     %(exponentiated) fixed cost of maintaining relationship
de              = mm.eta;   %demand elasticity
pi_tol          = mm.pi_tolerance;    % stopping rule for loop below 
n_z             = 2*mm.z_size+1;         % number of demand shocks states
L_b              = mm.L_b;           %shipment hazard.


diag_Q = abs(diag(Q0)); %the diagonal of Q0 (all the same as currently set up)
h = rh + del + L_b + diag_Q(1); %hazard of "something" happening (NOTE: This only
%works if all values in diag_Q are equal)

%matrix summation: needed to evaluate continuation value (NOTE: This only
%works if all values in diag_Q are equal)
mat_tol = 1e-8; %stopping criterion
ds_trans = zeros(size(Q_z));
n = 0;
diff = 1;
while diff > mat_tol
    add_me = (Q_z^n)/(h^n);
    ds_trans = ds_trans + add_me; 
    diff = max(max(add_me));
    n = n+1;
end

% payoff to shipment
payoff = 1 / de * exp(sf) * exp((de-1)*st(1,:)+st(2,:)); %st(1,:) is productivity, st(2,:) is macro state

% get expected profits for each type of buyer 
pi_z = bsxfun(@times,ones(size(Q0,1),n_z),payoff'); %initial guess on profits
pi_z_new = ones(size(pi_z)); 
c_val = ones(size(Q0,1),n_z); %continuation value
c_val_new = c_val;
for k=1:n_z
    eps = 1;
    it = 0;
    while eps > pi_tol && it<5000   % Iterate on contraction
        it = it + 1;
        if it == 5000
            display('WARNING: makepie.m did not converge!');
        end
        c_val_new(:,k) = Q0_d*c_val(:,k)+ L_b * pi_z * ds_trans(k,:)'; %first term: likelihood of self prod or macro change times * exp lifetime value, second term: value if shipment is next event 
        c_val_new(:,k) = max(-F + (c_val_new(:,k)/h),0);
        pi_z_new(:,k) = payoff'*exp(Z(k))+c_val_new(:,k); %create new pi_z
        eps = norm(pi_z(:,k)-pi_z_new(:,k))/norm(pi_z_new(:,k)); %get eps
        pi_z(:,k) = pi_z_new(:,k); %set old pi_z equal to new pi_z 
        c_val(:,k) = c_val_new(:,k); %set old c_val equal to new c_val
    end
end 

%expectation of profits over buyer states 
pi = pi_z*erg_pz; 

end %end function
