function [] = dens_by_cohort(dat,c_num,c_max);
% This function plots densities of sales by cohort

    step = round(c_max/c_num); %cohorts between plot points
    colors = hsv(c_num); %get colors for each histogram

    %remove -1's
    dat = max(dat,0);
    
    %First sum up sales in years inbetween steps
    sum_dat = zeros(c_num,size(dat,2));
    for k = 1:size(sum_dat,1)
        sum_dat(k,:) = sum(dat(step * (k-1) + 1:step * (k-1) + step,:),1);
    end

    %Now create the kernal density for each of the sums
    leg=zeros(3,2); % legend entries
    for k = 1:size(sum_dat,1)
        makekerndens(sum_dat(k,sum_dat(k,:)>0),colors(k,:)); % kernal density
        hold on; %keep all densities on the same access
        leg(k,1) = round(mean(sum_dat(k,sum_dat(k,:)>0)));
        leg(k,2) = sum(sum_dat(k,:)>0,2);
    end
    hold off;

    legend(horzcat('1-3 yrs, $', num2str(leg(1,1)),', ',num2str(leg(1,2)), ' exporters'),horzcat('4-6 yrs, $', num2str(leg(2,1)),', ',num2str(leg(2,2)),' exporters'),horzcat('7-9 yrs, $', num2str(leg(3,1)),', ',num2str(leg(3,2)),' exporters'));

end
