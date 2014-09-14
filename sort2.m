%This snippet performs the second large sort of the singlefirm script

[ind(1:obin,:),I] = sortrows((ind(1:obin,:)),1);
succ_prob(1:obin,:) = succ_prob(I,:);
ds(1:obin,:) = ds(I,:);
deathmat(1:obin,:) = deathmat(I,:);
cost_vec(1:obin,:) = cost_vec(I,:);
ind(1,5:12)=zeros(1,8);
for t = 2:4
    [loop_ind,~] = find((ind(2:obin,t)==0));
    if isempty(loop_ind) == 0
        for k = 1:size(loop_ind,1)
            ind(loop_ind(k)+1,t) = ind(loop_ind(k),t);
        end
    end
end
for t = 5:10
    [loop_ind,~] = find((ind(2:obin,t))==-1);
    if isempty(loop_ind) == 0
        for k = 1:size(loop_ind,1)
            temp = ind(loop_ind(k),t);
            ind(loop_ind(k)+1,t) = temp;        
        end
    end
end
for t = 11:12
    [loop_ind,~] = find((ind(2:obin,t))==-1);
    if isempty(loop_ind) == 0
        for k = 1:size(loop_ind,1)
            ind(loop_ind(k)+1,t) = 0;  
        end
    end
end
for t = 1:4
    [loop_ind,~] = find((succ_prob(2:obin,t)==0));
    if isempty(loop_ind) == 0
        for k = 1:size(loop_ind,1)
            succ_prob(loop_ind(k)+1,t) = succ_prob(loop_ind(k),t);
        end
    end
end
