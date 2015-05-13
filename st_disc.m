function [cli_no,sales_h,sales_f,ship_f,sh_ann_f,sh_first_yr_dum,cost_h,cost_f,succ_prob,prods] = st_disc(st_ind_cont,sale_h_cont,sale_f_cont,S,TT,burn,sh,maxc,sh_val_h_cont,sh_val_f_cont,cost_vec_cont,Phi,succ_prob_cont)
%this function takes the continuous versions of state and sale vectors, and
%collapses them into aggregate annual vectors

%preallocate
cli_no = mat2cell(repmat(zeros(TT-burn,2),S,1),ones(S,1)*TT-burn,2);
sales_h = mat2cell(repmat(zeros(TT-burn,1),S,1),ones(S,1)*TT-burn,1);
sales_f = mat2cell(repmat(zeros(TT-burn,1),S,1),ones(S,1)*TT-burn,1);
ship_f = mat2cell(repmat(zeros(TT-burn,1),S,1),ones(S,1)*TT-burn,1);
cost_f = mat2cell(repmat(zeros(TT-burn,1),S,1),ones(S,1)*TT-burn,1);
cost_h = mat2cell(repmat(zeros(TT-burn,1),S,1),ones(S,1)*TT-burn,1);
sh_ann_f = cell(S,1);
sh_first_yr_dum = cell(S,1);
prods = zeros(S,1);
succ_prob = zeros(S,4);
    
for j = 1:S

    % Read in productivity and success probabilities
    if isempty(st_ind_cont{j}) == 0
        prods(j) = Phi(st_ind_cont{j}(1,2));
        succ_prob_full = full(succ_prob_cont{j}); %create full version instead of sparse
        succ_prob(j,:) = succ_prob_full(1,:); %read in first row (all rows are the same!)
    else 
        prods(j) = -1;
        succ_prob(j,:) = [-1,-1,-1,-1];
    end

    t_lag = find(st_ind_cont{j}(:,1)<burn,1,'last'); %find index of last pre burn event
    if isempty(t_lag) == 1; %deal with firms which die before burn ends
        t_lag = 0;
    end

    % get a list of foreign clients to allocate a matrix to sh_ann_f 
    t_ind = find(st_ind_cont{j}(:,1)<TT,1,'last'); %find index of last pre TT event
    fclients = unique(full(sh_val_f_cont{j}(t_lag+1:t_ind,2))); %find list of all unique clients in the 'data' period
    fclients = fclients(2:end); % get rid of the 'zero' client, which is just a placeholder 
    sh_ann_f{j} = ones(TT-burn+1,size(fclients,1))*NaN; % extra row is for NaN, will be useful a barrier for stacking later when I calculate the moments, column is the number of clients!
    sh_first_yr_dum{j} = ones(TT-burn+1,size(fclients,1))*NaN; % this holds a dummy for the first year of a relationship

    % Calculate summable flow search costs (later can easily add fixed costs)
    gaps = st_ind_cont{j}(2:end,1) - st_ind_cont{j}(1:end-1,1); %get time passed in each period
    cost_f_sumable = [0;gaps .* cost_vec_cont{j}(1:end-1,1)]; %multiply gaps by instantaneous flow cost
    cost_f_sumable = cost_f_sumable + cost_vec_cont{j}(:,2); %add in fixed costs
    cost_h_sumable = [0;gaps .* cost_vec_cont{j}(1:end-1,3)]; %multiply gaps by instantaneous flow cost
    cost_h_sumable = cost_h_sumable + cost_vec_cont{j}(:,4); %add in fixed costs


    for t = burn+1:TT
        t_ind = find(st_ind_cont{j}(:,1)<t,1,'last');
        if isempty(t_ind) == 1
            t_ind = 0;
        end
        if t_lag~=t_ind %needed so that large gaps between events don't mess things up.
            t_lag = t_lag+1;
            cli_vec_h = sh{j}(t_lag:t_ind,2*maxc+1);
            cli_vec_f = sh{j}(t_lag:t_ind,2*maxc+2);
            cli_no{j}(t-burn,1) = size(unique(cli_vec_h(cli_vec_h>0)),1); %count the number of unique clients shipped to
            cli_no{j}(t-burn,2) = size(unique(cli_vec_f(cli_vec_f>0)),1);
            ship_f{j}(t-burn) = sum(sum(sh{j}(t_lag:t_ind,maxc+1:2*maxc)));
            sales_h{j}(t-burn) = sum(sale_h_cont{j}(t_lag:t_ind,1)); %sum of sales from last period to this period
            sales_f{j}(t-burn) = sum(sale_f_cont{j}(t_lag:t_ind,1)); %sum of sales from last period to this period
            cost_f{j}(t-burn) = sum(cost_f_sumable(t_lag:t_ind,1)); %sum of search (and fixed) costs by year 
            cost_h{j}(t-burn) = sum(cost_h_sumable(t_lag:t_ind,1)); %sum of search (and fixed) costs by year 
            if sales_f{j}(t-burn)==0 && cli_no{j}(t-burn,2) > 0
                display('WARNING: in st_disc, positive clients, but zero sales!');
            end

            % Calculate annual shipments
            if isempty(fclients) ~= 1 %did the firm have any clients?
                for k = 1:size(fclients,1)
                    ind = find(sh_val_f_cont{j}(t_lag:t_ind,2) == fclients(k));
                    if isempty(ind) == 1
                        sh_ann_f{j}(t-burn,k) = NaN;
                    else
                        sh_ann_f{j}(t-burn,k) = sum(sh_val_f_cont{j}(t_lag - 1 + ind,1));
                        if sh_ann_f{j}(t-burn,k) == 0
                            display('Zeros popping up in annual shipment stuffs')
                        end
                        sh_first_yr_dum{j}(t-burn,k) = sum(sh_val_f_cont{j}(t_lag - 1 + ind,3));
                    end
                end
            end

            t_lag = t_ind;
        end
    end
end
end
