function [] = makekerndens(x,color)
% Make a histogram

    [y, x] = ksdensity(log(x),'width',0.9);
    plot(x, y,'Color',color,'LineWidth',2);
    xlabel('log sales');
    ylabel('Density');

end
