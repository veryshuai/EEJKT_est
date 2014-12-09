function [] = makekerndens(x,color)
% Make a histogram

    [y, x] = ksdensity(log(x),'width',1);
    plot(x, y,'Color',color,'LineWidth',2);
    xlabel('log sales');
    ylabel('Density');
    axis([3 20 0 0.3])

end
