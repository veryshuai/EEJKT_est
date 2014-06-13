function res = prof_var(full_sim_dat)
% This function creates a smoothed 3D plot of profit variance conditional
% productivity phi and success probability theta

    % Read in data
    dat = full_sim_dat{10};
    prods = full(dat(:,1));
    mu_f = full(dat(:,2));
    
    % Get unique values for profits and success probabilities
    prod_list = full(unique(prods));
    mu_f_list = full(unique(mu_f));
    
    % Get foreign and domestic sales data
    sale_f_mat = full_sim_dat{4}; %foreign data
    sale_h_mat = full_sim_dat{2}; %domestic data
    
    %Preallocate
    rev_coef_var = zeros(size(prod_list,1),size(mu_f_list,1)) * inf; %preallocate
    freq_counts = zeros(size(prod_list,1),size(mu_f_list,1)); %preallocate
    
    % Get Productivity Coef of Variation
    for kind = 1:size(prod_list) %loop through productivities
        for jind = 1:size(mu_f_list) %loop through success probs 

            %only the firms with exactly that theta
            rev_temp = sale_f_mat(1,prods == prod_list(kind)...
                & mu_f == mu_f_list(jind)); 

            %domestic sales for the same firms
            dom_rev  = sale_h_mat(1,prods == prod_list(kind) ...
                & mu_f == mu_f_list(jind)); 

            %get total number of firms
            freq_counts(kind,jind) = size(rev_temp,2);

            if freq_counts(kind,jind) > 5  % remove cells with too few observations 
                if max(rev_temp(:) > 0) %coef of var doesn't exist if mean = 0
                    %only firms with either exports or domestic sales are counted 
                    rev_coef_var(kind,jind) = sqrt(var(rev_temp(rev_temp(:) > 0 | dom_rev(:) > 0)))...
                        / mean(rev_temp(rev_temp(:) > 0 | dom_rev(:) > 0));
                else
                    rev_coef_var(kind,jind) = 0; %Set coef of var to zero if no sales
                end
            end
        end
    end
    
    % Plot surface
    figure(1);
    mesh(rev_coef_var);

    % Tick and axis labeling
    set(gca,'XTick',1:7)
    set(gca,'XTickLabel',[0.1429, 0.2857, 0.4286, 0.5714, 0.7143, 0.8571, 0.9999])
    xlabel('Success Prob.')
    ylabel('Productivity')
    
end %end function
