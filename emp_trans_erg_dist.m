% This script computes the ergodic client distribution from the observed
% transition probabilities in Table 5

T = [0, 0, 0.932, 0.055, 0.009, 0.002, 0.001, 0.001, 0.000001;
    0, 0, .876, .1, .014, .008, 0, 0, 0.000001;
    .54, .08, .321, .048, .01, .002, 0, .001, 0;
    .194, .077, .375, .241, 0, .024, .009, .004, 0;
    .09, .042, .22, .271, .21, .092, 0, .027, 0;
    .059, 0, .129, .216, .215, .184, .083, .095, 0;
    0, 0, .095, .184, .181, .181, .126, .178, 0;
    0, 0, .039, .073, .089, .123, .157, .419, .073;
    0.00000001, 0.000000001, 0.000000001, 0.000000001, 0, 0, 0, .432, .526];

%Correct for missing data
for k=1:size(T,1)
    num_zeros = sum(T(k,:) == 0);
    new_val = (1 - sum(T(k,:))) / num_zeros;
    T(k,T(k,:) == 0) = new_val;
    T(k,:) = abs(T(k,:)) / sum(abs(T(k,:))); %ensure positive and sums to exactly one
end

%Get ergodic distributions
dist = ones(size(T,1)) / sum(ones(size(T,1)));
for k=1:100000
    dist = T' * dist;
end

%We only care about firms which are in the data
rel_dist = dist(3:end);
rel_dist = rel_dist / sum(rel_dist);
display(rel_dist)

