%in this section, we check the continuation value for a particular state to %see if it is less than zero, which is denoted by zero in our simulation.  If the continuation value is zero, the relationship is immediately endogenously dropped.

%home    
deathind_before_end = find(deathmat(:) == 1);
deathind = [deathind_before_end;obin-1];
lag = 1;
for n = 1:size(deathind,1)
    for t = 1:maxc
        loop_ind = find((ds(lag:deathind(n),end-1)) == t);
        loop_ind = loop_ind + lag - 1;
        if isempty(loop_ind) == 0 
            for o = 1:size(loop_ind,1)
                k = loop_ind(o);
            rel_time = ind(k,11); %time before exogenous separation
            p = k+1;
            dropped = 0; %dummy changes to one if relationship is dropped
            while p<deathind(n) && (ind(p,1)-ind(k,1)<rel_time-1e-12) %note: I had to take a little bit off here to deal with rounding error
                if dropped == 0 && sh(p-1,t) ~= 0; %sale is made
                    if ind(k-1,1) < shock_year % before macro policy change
                        if c_val_h_orig(ds(p-1,t),ind(p-1,2),ind(p-1,3))==0 %check for endogenous separation in the last period
                            sh(p,t) = 0; %kill any shipments in current period  
                            cost_vec(p,4) = 0;
                            ind(p,9) = ind(p,9)-1; %reduce current clients by one
                            dropped = 1;
                        end
                    else
                        if c_val_h_new(ds(p-1,t),ind(p-1,2),ind(p-1,3))==0 %check for endogenous separation in the last period
                            sh(p,t) = 0; %kill any shipments in current period  
                            cost_vec(p,4) = 0;
                            ind(p,9) = ind(p,9)-1; %reduce current clients by one
                            dropped = 1;
                        end
                    end
                elseif dropped == 1   
                    sh(p,t) = 0; %kill further shipments if relationship is over
                    cost_vec(p,4) = 0;
                    ind(p,9) = ind(p,9)-1; %reduce current clients by one
                end
                p = p+1;
            end
            temp = find((sh(1:deathind(n),t))==(ind(k,5)),1,'last');
            temp = temp(1,1);
            ind(k,11) = ind(temp,1)-ind(k,1)+1e-12;
            end
        end
    end
    sh(lag:deathind(n),2*maxc+1) = sum(sh(lag:deathind(n),1:maxc),2); %fill in penultimate column (shipment identities)
    sh(lag:deathind(n),1:maxc) = (sh(lag:deathind(n),1:maxc))>0; %make shipments binary (1 for a shipment)

    %foreign
     for t = 1:maxc
         loop_ind = find((ds(lag:deathind(n),end)) == maxc+t);
         loop_ind = loop_ind + lag - 1;
          if isempty(loop_ind) == 0 
            for o = 1:size(loop_ind,1)
                k = loop_ind(o);
            rel_time = ind(k,12); %time before exogenous separation
            p = k+1;
            dropped = 0; %dummy changes to one if relationship is dropped
            while p<deathind(n) && (ind(p,1)-ind(k,1)<rel_time-1e-12) %note: I had to take a little bit off here to deal with rounding error
                if dropped == 0 && sh(p-1,maxc+t) ~= 0; %sale is made
                    if ind(k-1,1) < shock_year % before macro policy change
                        if c_val_f_orig(ds(p-1,maxc+t),ind(p-1,2),ind(p-1,4))==0 %check for endogenous separation in the last period 
                            sh(p,maxc+t) = 0; %kill any sales in current period  
                            cost_vec(p,2) = 0;
                            ind(p,10) = ind(p,10)-1; %reduce current clients by one
                            dropped = 1;
                        end
                    else
                        if c_val_f_new(ds(p-1,maxc+t),ind(p-1,2),ind(p-1,4))==0 %check for endogenous separation in the last period 
                            sh(p,maxc+t) = 0; %kill any sales in current period  
                            cost_vec(p,2) = 0;
                            ind(p,10) = ind(p,10)-1; %reduce current clients by one
                            dropped = 1;
                        end
                    end
                elseif dropped == 1   
                    sh(p,maxc+t) = 0; %kill further sales if relationship is over
                    cost_vec(p,2) = 0;
                    ind(p,10) = ind(p,10)-1; %reduce current clients by one
                end
                p = p+1;
            end
            %display([k,t]);
            %display(ind(find(sh(:,t)==ind(k,5),1,'last'),1)-ind(k,1)+1e-12);
            temp = find((sh(1:deathind(n),maxc+t))==(ind(k,7)),1,'last');
            temp = temp(1,1);
            ind(k,12) =  ind(temp,1)-ind(k,1)+1e-12;
            end
          end
      end
    sh(lag:deathind(n),2*maxc+2) = sum(sh(lag:deathind(n),maxc+1:2*maxc),2); %fill in last column (shipment identities)
    sh(lag:deathind(n),maxc+1:2*maxc) = (sh(lag:deathind(n),maxc+1:2*maxc))>0; %make shipments binary (1 for a shipment)
    
    lag = deathind(n)+1;
end
        
%now kill the sales
sh_val_h(:,1:2) = [sh_val_h(:,1).*sum(sh(:,1:maxc),2),sh_val_h(:,2).*sum(sh(:,1:maxc),2)];
sh_val_f(:,1:2) = [sh_val_f(:,1).*sum(sh(:,maxc+1:2*maxc),2),sh_val_f(:,2).*sum(sh(:,maxc+1:2*maxc),2)];
