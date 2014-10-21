function [] = makehist(x,color)
% Make a histogram

    numOfBins = 20;
    [histFreq, histXout] = hist(log(x), numOfBins);
    %figure;
    bar(histXout, histFreq/sum(histFreq)*100,'FaceColor',color,'EdgeColor','none');
    xlabel('');
    ylabel('Frequency (percent)');
    alpha(0.7)

end
