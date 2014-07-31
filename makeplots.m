function makeplots(loadname, savename)

% This function creates plots for a counter-factual

    %get original data
    orig = read_plot_data('orig_decomp');
    
    %get cf data
    cf = read_plot_data(loadname);
    
    % subplots
    
    h = figure('Position', [100, 100, 1049, 895]);
    
    subplot(2,2,1)
    plot(orig{1},'r','LineWidth',2);
    hold on;
    plot(cf{1},'--','LineWidth',2);
    ylabel('Log Total Exports')
    set(gca,'XTickLabel',['1992';'1997';'2002';'2007';'2012'])
    hold off;
    
    subplot(2,2,2)
    plot(orig{2},'r','LineWidth',2);
    hold on;
    plot(cf{2},'--','LineWidth',2);
    ylabel('Log Active Exporters')
    set(gca,'XTickLabel',['1992';'1997';'2002';'2007';'2012'])
    hold off;
    
    subplot(2,2,3)
    plot(orig{3},'r','LineWidth',2);
    hold on;
    plot(cf{3},'--','LineWidth',2);
    ylabel('Log Mean Sales per Client')
    set(gca,'XTickLabel',['1992';'1997';'2002';'2007';'2012'])
    hold off;
    
    subplot(2,2,4)
    plot(orig{4},'r','LineWidth',2);
    hold on;
    plot(cf{4},'--','LineWidth',2);
    ylabel('Log Mean Number of Clients')
    set(gca,'XTickLabel',['1992';'1997';'2002';'2007';'2012'])
    hold off;
    
    print(h,'-depsc',savename)
    
end
