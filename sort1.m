%This snippet perfors the first large sort of the singlefirm script

%first part eliminates times greater than TT, and resets obin
obin = find((ind(1:obin,1))>0,1,'last');
obin = obin(1,1);
deathmat(obin+1:end) = 0;
scalar_temp = numel(obin+1:mms);
ind(obin+1:end,:) = repmat([0,0,0,0,-1,-1,-1,-1,zeros(1,2),-1,-1],scalar_temp,1);
[ind(1:obin,:),I] = sortrows((ind(1:obin,:)),1);
succ_prob(1:obin,:) = succ_prob(I,:);
deathmat(1:obin,1) = deathmat(I,1); 
cost_vec(1:obin,1) = cost_vec(I,1); 
for t = 2:4
    [loop_ind,~] = find(ind(2:obin,t)==0);
    for k = 1:size(loop_ind,1)
        temp = ind(loop_ind(k),t);
        ind(loop_ind(k)+1,t) = temp;
    end
end   
for t = 1:4
    [loop_ind,~] = find(succ_prob(2:obin,t)==0);
    for k = 1:size(loop_ind,1)
        temp = succ_prob(loop_ind(k),t);
        succ_prob(loop_ind(k)+1,t) = temp;
    end
end   
