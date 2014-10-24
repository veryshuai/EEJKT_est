function [] = dens_by_cohort(dat,c_num,c_max);
% This function plots densities of sales by cohort

    step = round(c_max/c_num); %cohorts between plot points
    colors = hsv(c_num); %get colors for each histogram
    
    %First sum up sales in years inbetween steps
    sum_dat = zeros(c_num,size(dat,2));
    for k = 1:size(sum_dat,1)
        sum_dat(k,:) = sum(dat(step * (k-1) + 1:step * (k-1) + step,:),1);
    end

    for k = 1:size(sum_dat,1)
        %makehist(sale_f_coh(k,sale_f_coh(k,:)>0),colors(k,:)); % plot
        makekerndens(sum_dat(k,sum_dat(k,:)>0),colors(k,:));
        hold on;
    end
    hold off;

end
