function makeplots(loadname, savename)

% This function creates plots for a counter-factual

    %get original data
    orig = read_plot_data('results/orig_decomp');
    
    %get cf data
    cf = read_plot_data(loadname);
    
    % subplots
    
    h = figure('Position', [100, 100, 1049, 895]);
    
    subplot(2,2,1)
    plot(orig{4}(1:30),'r','LineWidth',2);
    hold on;
    plot(cf{4}(1:30),'--','LineWidth',2);
    ylabel('Log Total Exports')
    set(gca,'XTickLabel',['';'';'';'';''])
    hold off;
    
    subplot(2,2,2)
    plot(orig{1}(1:30),'r','LineWidth',2);
    hold on;
    plot(cf{1}(1:30),'--','LineWidth',2);
    ylabel('Log Active Exporters')
    set(gca,'XTickLabel',['';'';'';'';''])
    hold off;
    
    subplot(2,2,3)
    plot(orig{2}(1:30),'r','LineWidth',2);
    hold on;
    plot(cf{2}(1:30),'--','LineWidth',2);
    ylabel('Log Mean Sales per Client')
    set(gca,'XTickLabel',['';'';'';'';''])
    hold off;
    
    subplot(2,2,4)
    plot(orig{3}(1:30),'r','LineWidth',2);
    hold on;
    plot(cf{3}(1:30),'--','LineWidth',2);
    ylabel('Log Mean Number of Clients')
    set(gca,'XTickLabel',['';'';'';'';''])
    hold off;
    
    print(h,'-depsc',savename)
    
end
