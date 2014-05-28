function [noe,cpe,spe,te,tc,spc] = dynam(sale_f_mat,cli_no_mat)
%This function is used to learn about the dynamics of exporting.  It should
%be run after stopping the main calibration program at the end of the file
%moms.m

cli_no = reshape(cli_no_mat(:,2),size(sale_f_mat,1),size(sale_f_mat,2));

noe = zeros(size(sale_f_mat,1),1);
te = zeros(size(noe));
cpe = zeros(size(noe));
tc = zeros(size(noe));
for k = 1:size(sale_f_mat,1)
    noe(k) = sum(sale_f_mat(k,:)>0); %number of exporters in each year
    te(k) = sum(sale_f_mat(k,sale_f_mat(k,:)>0)); %total exports 
    cpe(k) = mean(cli_no(k,sale_f_mat(k,:)>0)); %mean clients per exporter
    tc(k) = sum(cli_no(k,sale_f_mat(k,:)>0));
end
spe = te./noe;
spc = te./tc;
end

