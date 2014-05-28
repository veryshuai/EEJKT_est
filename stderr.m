
clc

load bcov
clearvars -except bcov X

if matlabpool('size')<12
    matlabpool open 12
end

% random seed, needed to call pop values
rng(80085);

tic;
 addpath(genpath('/gpfs/home/dcj138/work/Colombia-Project/Colombia-Project/'));
    
 %X = [...
 %  0.0844793112802   0.1544787808223   9.7991495253733 0.5991575264434   0.0744631426408   0.2949275622659 0.9999115885627   0.1477909587384   0.3023232808900 4.6315563382334   0.0376227842532   136.1032622469932 0.4936345118741;...
 %  ];

 
pop = X; %use the same parameters as in the bootstrap

findif = 1e-3;

jac = zeros(42,13);

% important matlab property is that rng random seeding 
% functions different in parallel loops.
% that is why we need to call this in parfor...
W = zeros(42,42,1);
base = zeros(42,1);
parfor k = 1:1
  [~,W_temp,base_temp] = distance_noprod(pop(k,:));
  W(:,:,k) = W_temp;
  base(:,k) = base_temp;
end
W = squeeze(W);
base = squeeze(base);

basis = eye(size(pop,2));

parfor k = 1:size(pop,2)
    display(k);
    ptemp = pop .* (1 + basis(k,:) * findif);
    [~,~,new] = distance_noprod(ptemp);
    jac(:,k) = (new-base) / (pop(k) * findif);
end

temp = (jac' * W * jac)^-1 * (jac' * W * 2 * bcov * W * jac) * (jac' * W * jac)^-1;
posdefchk = min(eig((temp+temp')/2));
if posdefchk<0
    display('WARNING: Not a local maximum');
end
se = sqrt(diag(temp));
varcov = temp;

toc
matlabpool close
save  se-5-3-2013
%diary off
