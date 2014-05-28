function [ac,an] = av_cost(ind,mm,lambda_f,th_ind)
%this function takes the event matrix ind, and parameter structure mm and
%returns actual average costs paid over time and firm

firm_cost = zeros(mm.S,1);
firm_net  = zeros(mm.S,1);
c = @(x,n) ((1+x).^(1+1/mm.b)-1)/((1+1/mm.b)*n^mm.gam); %cost function

for k = 1:mm.S
    if isempty(ind{k})==0
        for m = 2:size(ind{k},1)
            display([k,m]);
            l = lambda_f{min(mm.n_size,ind{k}(m-1,8))+1,min(mm.n_size,ind{k}(m-1,7))+1,th_ind(k,1),min(ind{k}(m-1,8),mm.net_size)+1}(ind{k}(m-1,2),ind{k}(m-1,3));
            cost = c(l,min(mm.net_size,ind{k}(m-1,8))+1)*(ind{k}(m,1)-ind{k}(m-1,1));
            firm_net(k) = firm_net(k) + (ind{k}(m,1)-ind{k}(m-1,1))*(min(ind{k}(m-1,8),mm.net_size)+1);
            firm_cost(k) = firm_cost(k) + cost;
        end
        firm_cost(k) = firm_cost(k)/(ind{k}(end,1)-ind{k}(1,1)); 
        firm_net(k) = firm_net(k)/(ind{k}(end,1)-ind{k}(1,1)); 
    end    
end
ac = sum(firm_cost(isnan(firm_cost)==0))/sum(isnan(firm_cost)==0);
an = sum(firm_net(isnan(firm_net)==0))/sum(isnan(firm_net)==0);
end