function [] = dens_by_cohort(dat);
% This function plots densities of sales by cohort

    c_num = 3; %number of cohorts to plot
    c_max = 7; %last cohort to plot
    colors = hsv(c_max); %get colors for each histogram
    for k = 1:round(c_max/c_num):c_max
        %makehist(sale_f_coh(k,sale_f_coh(k,:)>0),colors(k,:)); % plot
        makekerndens(dat(k,dat(k,:)>0),colors(k,:));
        hold on;
    end
    hold off;

end
