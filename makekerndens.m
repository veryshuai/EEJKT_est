function [] = makekerndens(x,color)
% Make a histogram

    [y, x] = ksdensity(log(x));%,'width',0.5);
    plot(x, y,'Color',color,'LineWidth',2);
    xlabel('log sales');
    ylabel('Density');
    axis([6 22 0 0.5])

end
