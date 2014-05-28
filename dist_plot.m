%Presentation Plot

yax = log(histc(cli_no_mat(cli_no_mat(:,2)>0,2),1:10));
xax = log(1:10);
scatter(xax,yax(1:10),'filled');
xlabel('log clients');
ylabel('log frequencey');
title('Simulated Log-Log Client Distribution-Network Effect');

figure;
yax = histc(log(sale_f_spc(sale_f_spc(:)>0)),2:16);
scatter(2:16,log(yax),'filled');
xlabel('log sales');
ylabel('log frequencey');
title('Simulated Log-Log Total Sales Distribution-Network Effect');