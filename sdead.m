function [tot_S] = sdead(S_old,t)
%this function simply puts new firms into new cells


flag = 0; tot_S = 0;

st_cont = [];

if flag == 0
    for j = 1:S_old %loop through the preset number of firm slots

        obin = 2; %current index

        % Load written files
        load(sprintf('temp_data/temp_%d_%d.mat', j, t))
        st_ind_cont_old = cind;
        ds_old = cds;
        sh_old = csh;
        deathmat = cdeathmat;
        sh_val_h_old = csh_val_h;
        sh_val_f_old = csh_val_f;
        cost_vec_old = ccost_vec;
        succ_prob_old = csucc_prob;
        
        % Find Deaths
        ind = find(deathmat==1); %find position of any death
        ind_ext = [ind;size(st_ind_cont_old,1)]; % append final slot of event matrix for firm j 
        
        % initialize new cell arrays (note old arrays called, well, old)
        S           = 1+sum(deathmat); % new number of firms (after death is included)
        tot_S       = tot_S + S; % running sum of total firm number
        st_ind_cont = cell(size(ind,1),1);
        ds          = cell(size(ind,1),1);
        sh          = cell(size(ind,1),1);
        sh_val_h    = cell(size(ind,1),1);
        sh_val_f    = cell(size(ind,1),1);
        cost_vec    = cell(size(ind,1),1);
        succ_prob   = cell(size(ind,1),1);
        

      % first firm in slot 1
      sh_val_h{1} = sh_val_h_old(1:ind_ext(1),:); %home shipment values
      sh_val_f{1} = sh_val_f_old(1:ind_ext(1),:); %foreign shipment values
      cost_vec{1} = cost_vec_old(1:ind_ext(1),:); %cost vector
      succ_prob{1} = succ_prob_old(1:ind_ext(1),:); %cost vector
      st_ind_cont{1} = st_ind_cont_old(1:ind_ext(1),:); %other firm related stuff
      ds{1} = ds_old(1:ind_ext(1),:); %demand shock (client hotel)
      sh{1} = sh_old(1:ind_ext(1),:); %shipments (client hotel)

      % loop through 2nd to nth firms in slot j
      lag = ind_ext(1)+1; %first event row of next firm in slot j
      for k = 2:size(ind) %EMPRICALLY USING IND_EXT WAS JUST MAKING AN EXTRA BLANK CELL
          sh_val_h{obin} = sh_val_h_old(lag:ind_ext(k),:);
          sh_val_f{obin} = sh_val_f_old(lag:ind_ext(k),:);
          cost_vec{obin} = cost_vec_old(lag:ind_ext(k),:);
          succ_prob{obin} = succ_prob_old(lag:ind_ext(k),:);
          st_ind_cont{obin} = st_ind_cont_old(lag:ind_ext(k),:);
          ds{obin} = ds_old(lag:ind_ext(k),:);
          sh{obin} = sh_old(lag:ind_ext(k),:);
          lag = ind_ext(k)+1; %increment first even row
          obin = obin+1; %increment firm number in slot j 
      end
      
      % Resave to disk
          save(sprintf('temp_data/temp_%d_%d.mat', j, t),'st_ind_cont' ,'ds' ,'sh' ,'sh_val_h' ,'sh_val_f' ,'cost_vec' ,'succ_prob','t','S')
      
  end
end


end
