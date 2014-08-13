function cf_decomposition(loadname, savename)
% This function decomposes the sales change from a counterfactualinto sales
% per client, number of active exporters, and number of clients.

    %load file
    load(loadname)

    %get file size
    sz = size(sale_f_mat,1);
    
    cli_no_temp     = cell2mat(cli_no'); % get client nos in a matrix like sale_f_mat
    cli_no_spc      = cli_no_temp(:,2:2:end);
    cli_no_dom      = cli_no_temp(:,1:2:end);
    spc             = sale_f_mat./cli_no_spc; %sales per client
    m_spc           = zeros(sz,1); 
    act_exp         = zeros(sz,1); 
    act_dom         = zeros(sz,1);
    m_numclients    = zeros(sz,1);
    numclients      = zeros(sz,1);
    
    for k = 1:sz
        sumable         = spc(k,:);
        cli             = cli_no_spc(k,:);
        cli_dom         = cli_no_dom(k,:);
        m_spc(k)        = sum(sumable(~isnan(sumable)))/sum(~isnan(sumable)); %mean sales per client, by year
        act_exp(k)      = sum(cli>0);
        act_dom(k)      = sum(cli_dom>0);
        numclients(k)   = sum(cli);
        m_numclients(k) = numclients(k)/sum(cli>0);
    end
    
    % decomposition (assume counterfactual intervention in period 10)
    
    Dba_sales         = log(sum(sum(sale_f_mat(10:sz,:),2))) - log(sum(sum(sale_f_mat(1:9,:),2)));
    Dba_act_exp       = log(sum(act_exp(10:sz)))-log(sum(act_exp(1:9)));
    Dba_m_numclients  = log(sum(numclients(10:sz))/sum(act_exp(10:sz))) - log(sum(numclients(1:9))/sum(act_exp(1:9)));
    Dba_m_spc         = log(sum(sum(sale_f_mat(10:sz,:),2))/sum(numclients(10:sz))) - log(sum(sum(sale_f_mat(1:9,:),2))/sum(numclients(1:9)));
    
    D_sales         = log(sum(sale_f_mat(2:sz,:),2)) - log(sum(sale_f_mat(1:sz-1,:),2));
    D_act_exp       = log(act_exp(2:sz))-log(act_exp(1:sz-1));
    D_m_numclients  = log(numclients(2:sz)./act_exp(2:sz)) - log(numclients(1:sz-1)./act_exp(1:sz-1));
    D_m_spc         = log(sum(sale_f_mat(2:sz,:),2)./numclients(2:sz)) - log(sum(sale_f_mat(1:sz-1,:),2)./numclients(1:sz-1));

    save(savename)
end