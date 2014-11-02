function [sale_h,sale_f] = sales_trans(scale_f_orig,scale_f_new,scale_h,de,st_cont,S,ds,sh,maxc,Z,Phi,X_h,X_f,TT)
%this function takes scale factors, elasticity of demand, and the simulated
%state transition vector and returns sales per period between events for
%each firm.

sale_h = cell(S,1);
sale_f = cell(S,1);

Z_big = [-inf,Z'];

for j = 1:S
   sale_h{j} = exp(scale_h)*exp(Phi(st_cont{j}(:,2))).^(de-1).*exp(X_h(st_cont{j}(:,3))).*sum(exp(Z_big(ds{j}(:,1:maxc)+1)).*sh{j}(:,1:maxc),2); 
   %for policy experiments, need to find the first period after the shock
   shock_time = find(st_cont{j}(:,1) < TT - 9,1,'last');
   scale_f = ones(size(st_cont{j}(:,1))) * exp(scale_f_new);
   if isempty(shock_time) == 0
       scale_f(1:shock_time) = exp(scale_f_orig);
   end
   sale_f{j} = scale_f .* exp(Phi(st_cont{j}(:,2))).^(de-1).*exp(X_f(st_cont{j}(:,4))).*sum(exp(Z_big(ds{j}(:,maxc+1:2*maxc)+1)).*sh{j}(:,maxc+1:2*maxc),2);  
    sale_f{j} = scale_f_vec .* exp(Phi(st_cont{j}(:,2))).^(de-1).*exp(X_f(st_cont{j}(:,4))).*sum(exp(Z_big(ds{j}(:,maxc+1:2*maxc)+1)).*sh{j}(:,maxc+1:2*maxc),2);  
end

end
