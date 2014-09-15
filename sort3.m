%This snippet performs the third large sort of the singlefirm script

[ind(1:obin,:),I] = sortrows((ind(1:obin,:)),1);
succ_prob(1:obin,:) = succ_prob(I,:);
ds(1:obin,:) = ds(I,:);
sh(1:obin,:) = sh(I,:);
sh_val_h(1:obin,:) = sh_val_h(I,:);
sh_val_f(1:obin,:) = sh_val_f(I,:);
deathmat(1:obin,:) = deathmat(I,:);
cost_vec(1:obin,:) = cost_vec(I,:);
    
for t = 2:4
    [loop_ind,~] = find((ind(2:obin,t))==0);
    if isempty(loop_ind)==0
        for k = 1:size(loop_ind,1)
            temp = ind(loop_ind(k),t);
            ind(loop_ind(k)+1,t) = temp;
        end
    end
end
for t = 1:4
    [loop_ind,~] = find((succ_prob(2:obin,t))==0);
    if isempty(loop_ind)==0
        for k = 1:size(loop_ind,1)
            temp = succ_prob(loop_ind(k),t);
            succ_prob(loop_ind(k)+1,t) = temp;
        end
    end
end
for t = 5:10
    [loop_ind,~] = find((ind(2:obin,t))==-1);
    for k = 1:size(loop_ind,1)
        temp = ind(loop_ind(k),t);
        ind(loop_ind(k)+1,t) = temp;        
    end
end
for t = 1:2*maxc
    loop_ind = find((ds(2:obin,t))==-1);
    for k = 1:size(loop_ind,1)
        ds(loop_ind(k)+1,t) = ds(loop_ind(k),t);
    end
end
Z_big = [-inf,Z'];
loop_ind = find(sh_val_h(:,2) ~= 0);
for k = 1:size(loop_ind,1)
    cli_ind = find(ds(loop_ind(k),:),1,'first');
    sh_val_h(loop_ind(k),1) = exp(scale_h)*exp(Phi(ind(loop_ind(k),2)))^(de-1)*exp(X_h(ind(loop_ind(k),3)))*exp(Z_big(ds(loop_ind(k),cli_ind)+1)); 
end
loop_ind = find(sh_val_f(:,2) ~= 0);
for k = 1:size(loop_ind,1)
    cli_ind = find(ds(loop_ind(k),:),1,'first');
    sh_val_f(loop_ind(k),1) = exp(scale_f)*exp(Phi(ind(loop_ind(k),2)))^(de-1)*exp(X_f(ind(loop_ind(k),3)))*exp(Z_big(ds(loop_ind(k),cli_ind)+1)); 
end

%fill in flow costs
for t = 1:2:4 
    cost_vec(1:3,1) = 0; %these are nonsense placeholders at beginning of state vector
    cost_vec(1:3,3) = 0; %these are nonsense placeholders at beginning of state vector

    % fill in
    loop_ind = find(cost_vec(2:obin,t) == -1);
    if isempty(loop_ind)==0
        for k = 1:size(loop_ind,1) %fill forward
            temp = cost_vec(loop_ind(k),t);
            cost_vec(loop_ind(k)+1,t) = temp;
        end
    end
end
