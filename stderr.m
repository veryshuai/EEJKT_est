
clc

load results/bcov
clearvars -except bcov X

if matlabpool('size')<7
    matlabpool open 7 
end

% random seed, needed to call pop values
rng(80085);

tic;

pop = X; %use the same parameters as in the bootstrap

findif = 1e-4;

jac = zeros(41,13);

% important matlab property is that rng random seeding 
% functions different in parallel loops.
% that is why we need to call this in parfor...
W = eye(size(jac,1));
base = zeros(41,1);
parfor k = 1:1
  [~,W_temp,base_temp] = distance_noprod(pop(k,:),0,1);
  base(:,k) = base_temp;
end
W = squeeze(W);
base = squeeze(base);

% Now reevaluate the objective for a small deviation in each parameter
% and get the empirical jacobian
basis = eye(size(pop,2));

parfor k = 1:size(pop,2)
    display(k);
    ptemp = pop .* (1 + basis(k,:) * findif);
    [~,~,new] = distance_noprod(ptemp,0,1);
    jac(:,k) = (new-base) / (pop(k) * findif);
end

% Calculate the standard errors
save results/temp
temp = (jac' * W * jac)^-1 * (jac' * W * 2 * bcov * W * jac) * (jac' * W * jac)^-1;
posdefchk = min(eig((temp+temp')/2));
if posdefchk<0
    display('WARNING: Not a local maximum');
end
se = sqrt(diag(temp));
varcov = temp;

toc
matlabpool close
save  results/se_results
%diary off
