% This script creates plots for counter-factuals

%orig
load orig_results;
act_exp_orig          = log(act_exp);
m_spc_orig            = log(m_spc);
m_numclients_orig     = log(m_numclients);
tot_exp_orig          = log(sum(sale_f_mat,2));

%mac_bump
load mac_bump_results;
act_exp_mac_bump      = log(act_exp);
m_spc_mac_bump        = log(m_spc);
m_numclients_mac_bump = log(m_numclients);
tot_exp_mac_bump      = log(sum(sale_f_mat,2));

%cost_dec
load cost_dec_results;
act_exp_cost_dec      = log(act_exp);
m_spc_cost_dec        = log(m_spc);
m_numclients_cost_dec = log(m_numclients);
tot_exp_cost_dec      = log(sum(sale_f_mat,2));

%red_var
load red_var_results;
act_exp_red_var      = log(act_exp);
m_spc_red_var        = log(m_spc);
m_numclients_red_var = log(m_numclients);
tot_exp_red_var      = log(sum(sale_f_mat,2));

%search_dec
load search_dec_results;
act_exp_search_dec      = log(act_exp);
m_spc_search_dec        = log(m_spc);
m_numclients_search_dec = log(m_numclients);
tot_exp_search_dec      = log(sum(sale_f_mat,2));

%no net
load no_net_results;
act_exp_no_net      = log(act_exp);
m_spc_no_net        = log(m_spc);
m_numclients_no_net = log(m_numclients);
tot_exp_no_net      = log(sum(sale_f_mat,2));

%plots
%h = figure('Position', [100, 100, 1049, 895]);
%plot(tot_exp_mac_bump - tot_exp_orig,'b','LineWidth',4);
%ylabel('Log Total Exports')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%print(h,'-depsc','mac_bump_tot_exp_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(act_exp_mac_bump - act_exp_orig,'b','LineWidth',4);
%ylabel('Log Active Exporters')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','mac_bump_act_exp_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(m_spc_mac_bump - m_spc_orig,'b','LineWidth',4);
%ylabel('Log Mean Sales per Client')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','mac_bump_mspc_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(m_numclients_mac_bump - m_numclients_orig,'b','LineWidth',4);
%ylabel('Log Mean Number of Clients')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','mac_bump_mnumcli_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(tot_exp_cost_dec - tot_exp_orig,'b','LineWidth',4);
%ylabel('Log Total Exports')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','cost_dec_tot_exp_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(act_exp_cost_dec - act_exp_orig,'b','LineWidth',4);
%ylabel('Log Active Exporters')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','cost_dec_act_exp_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(m_spc_cost_dec - m_spc_orig,'b','LineWidth',4);
%ylabel('Log Mean Sales per Client')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','cost_dec_mspc_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(m_numclients_cost_dec - m_numclients_orig,'b','LineWidth',4);
%ylabel('Log Mean Number of Clients')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','cost_dec_mnumcli_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(tot_exp_red_var - tot_exp_orig,'b','LineWidth',4);
%ylabel('Log Total Exports')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','red_var_tot_exp_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(act_exp_red_var - act_exp_orig,'b','LineWidth',4);
%ylabel('Log Active Exporters')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','red_var_act_exp_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(m_spc_red_var - m_spc_orig,'b','LineWidth',4);
%ylabel('Log Mean Sales per Client')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','red_var_mspc_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(m_numclients_red_var - m_numclients_orig,'b','LineWidth',4);
%ylabel('Log Mean Number of Clients')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','red_var_mnumcli_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(tot_exp_search_dec - tot_exp_orig,'b','LineWidth',4);
%ylabel('Log Total Exports')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','search_dec_tot_exp_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(act_exp_search_dec - act_exp_orig,'b','LineWidth',4);
%ylabel('Log Active Exporters')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','search_dec_act_exp_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(m_spc_search_dec - m_spc_orig,'b','LineWidth',4);
%ylabel('Log Mean Sales per Client')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','search_dec_mspc_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(m_numclients_search_dec - m_numclients_orig,'b','LineWidth',4);
%ylabel('Log Mean Number of Clients')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','search_dec_mnumcli_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(tot_exp_orig,'b','LineWidth',4);
%hold on;
%plot(tot_exp_no_net,'--','LineWidth',4);
%% title('No Network Effects')
%ylabel('Log Total Exports')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','no_net_tot_exp_diff.eps')

%h = figure('Position', [100, 100, 1049, 895]);
%plot(act_exp_orig,'b','LineWidth',4);
%hold on;
%plot(act_exp_no_net,'--','LineWidth',4);
%% title('No Network Effects')
%ylabel('Log Active Exporters')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','no_net_act_exp_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(m_spc_orig,'b','LineWidth',4);
%hold on;
%plot(m_spc_no_net,'--','LineWidth',4);
%% title('No Network Effects')
%ylabel('Log Mean Sales per Client')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','no_net_mspc_diff.eps')
%
%h = figure('Position', [100, 100, 1049, 895]);
%plot(m_numclients_orig,'b','LineWidth',4);
%hold on;
%plot(m_numclients_no_net,'--','LineWidth',4);
%% title('No Network Effects')
%ylabel('Log Mean Number of Clients')
%set(gca,'XTickLabel',['1992';'1994';'1996';'1998';'2000';'2002';'2004';'2006';'2008';'2010'])
%ylim([-0.1,0.5])
%hold off;
%print(h,'-depsc','no_net_mnumcli_diff.eps')

% subplots

h = figure('Position', [100, 100, 1049, 895]);

subplot(2,2,1)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(tot_exp_mac_bump - tot_exp_orig,'b','Linewidth',3);
ylabel('Log Total Exports')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

subplot(2,2,2)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(act_exp_mac_bump - act_exp_orig,'b','Linewidth',3);
ylabel('Log Active Exporters')
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

subplot(2,2,3)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(m_spc_mac_bump - m_spc_orig,'b','Linewidth',3);
ylabel('Log Mean Sales per Client')
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

subplot(2,2,4)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(m_numclients_mac_bump - m_numclients_orig,'b','Linewidth',3);
ylabel('Log Mean Number of Clients')
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

% title('Macro Shock Level Increase')
print(h,'-depsc','mac_bump_subplot_diff.eps')

h = figure('Position', [100, 100, 1049, 895]);

subplot(2,2,1)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(tot_exp_cost_dec - tot_exp_orig,'b','Linewidth',3);
ylabel('Log Total Exports')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

subplot(2,2,2)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(act_exp_cost_dec - act_exp_orig,'b','Linewidth',3);
ylabel('Log Active Exporters')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

subplot(2,2,3)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(m_spc_cost_dec - m_spc_orig,'b','Linewidth',3);
ylabel('Log Mean Sales per Client')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

subplot(2,2,4)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(m_numclients_cost_dec - m_numclients_orig,'b','Linewidth',3);
ylabel('Log Mean Number of Clients')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

% title('Fixed Cost Decrease')
print(h,'-depsc','cost_dec_subplot_diff.eps')

h = figure('Position', [100, 100, 1049, 895]);

subplot(2,2,1)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(tot_exp_red_var - tot_exp_orig,'b','Linewidth',3);
ylabel('Log Total Exports')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

subplot(2,2,2)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(act_exp_red_var - act_exp_orig,'b','Linewidth',3);
ylabel('Log Active Exporters')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

subplot(2,2,3)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(m_spc_red_var - m_spc_orig,'b','Linewidth',3);
ylabel('Log Mean Sales per Client')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

subplot(2,2,4)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(m_numclients_red_var - m_numclients_orig,'b','Linewidth',3);
ylabel('Log Mean Number of Clients')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

% title('Macro Shock Variance Decrease')
print(h,'-depsc','red_var_subplot_diff.eps')

h = figure('Position', [100, 100, 1049, 895]);

subplot(2,2,1)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(tot_exp_search_dec - tot_exp_orig,'b','Linewidth',3);
ylabel('Log Total Exports')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

subplot(2,2,2)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(act_exp_search_dec - act_exp_orig,'b','Linewidth',3);
ylabel('Log Active Exporters')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

subplot(2,2,3)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(m_spc_search_dec - m_spc_orig,'b','Linewidth',3);
ylabel('Log Mean Sales per Client')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

subplot(2,2,4)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(m_numclients_search_dec - m_numclients_orig,'b','Linewidth',3);
ylabel('Log Mean Number of Clients')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

% title('Search Cost Decrease')
print(h,'-depsc','search_dec_subplot_diff.eps')

h = figure('Position', [100, 100, 1049, 895]);

subplot(2,2,1)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(tot_exp_no_net - tot_exp_orig,'b','Linewidth',3);
ylabel('Log Total Exports')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

subplot(2,2,2)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(act_exp_no_net - act_exp_orig,'b','Linewidth',3);
ylabel('Log Active Exporters')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

subplot(2,2,3)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(m_spc_no_net - m_spc_orig,'b','Linewidth',3);
ylabel('Log Mean Sales per Client')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

subplot(2,2,4)
plot(tot_exp_cost_dec * 0, 'r--','Linewidth',1)
hold on 
plot(m_numclients_no_net - m_numclients_orig,'b','Linewidth',3);
ylabel('Log Mean Number of Clients')
L = get(gca, 'XLim');
set(gca,'XTick',linspace(L(1),L(2),5))
set(gca,'XTickLabel',['1992';'1997';'2002';'20.8';'2012'])
% ylim([-0.1,0.8])
hold off

% title('No Network Effects')
print(h,'-depsc','no_net_subplot_diff.eps')
