function cf_decomposition(savename)
% This function decomposes the sales change from a counterfactualinto sales
% per client, number of active exporters, and number of clients.

    load cf_sim_results
    
    cli_no_temp     = cell2mat(cli_no'); % get client nos in a matrix like sale_f_mat
    cli_no_spc      = cli_no_temp(:,2:2:end);
    cli_no_dom      = cli_no_temp(:,1:2:end);
    spc             = sale_f_mat./cli_no_spc; %sales per client
    m_spc           = zeros(18,1); act_exp         = zeros(18,1); act_dom         = zeros(18,1);
    m_numclients    = zeros(18,1);
    numclients      = zeros(18,1);
    
    for k = 1:18
        sumable         = spc(k,:);
        cli             = cli_no_spc(k,:);
        cli_dom         = cli_no_dom(k,:);
        m_spc(k)        = sum(sumable(~isnan(sumable)))/sum(~isnan(sumable)); %mean sales per client, by year
        act_exp(k)      = sum(cli>0);
        act_dom(k)      = sum(cli_dom>0);
        numclients(k)   = sum(cli);
        m_numclients(k) = numclients(k)/sum(cli>0);
    end
    
    % decomposition
    
    Dba_sales         = log(sum(sum(sale_f_mat(10:18,:),2))) - log(sum(sum(sale_f_mat(1:9,:),2)));
    Dba_act_exp       = log(sum(act_exp(10:18)))-log(sum(act_exp(1:9)));
    Dba_m_numclients  = log(sum(numclients(10:18))/sum(act_exp(10:18))) - log(sum(numclients(1:9))/sum(act_exp(1:9)));
    Dba_m_spc         = log(sum(sum(sale_f_mat(10:18,:),2))/sum(numclients(10:18))) - log(sum(sum(sale_f_mat(1:9,:),2))/sum(numclients(1:9)));
    
    D_sales         = log(sum(sale_f_mat(2:18,:),2)) - log(sum(sale_f_mat(1:17,:),2));
    D_act_exp       = log(act_exp(2:18))-log(act_exp(1:17));
    D_m_numclients  = log(numclients(2:18)./act_exp(2:18)) - log(numclients(1:17)./act_exp(1:17));
    D_m_spc         = log(sum(sale_f_mat(2:18,:),2)./numclients(2:18)) - log(sum(sale_f_mat(1:17,:),2)./numclients(1:17));
    
    save savename
end
