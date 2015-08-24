function [st_cont,st_ind_cont,S,ds,sh,flag,sh_val_h,sh_val_f,cost_vec,succ_prob] = sdead(st_cont,st_ind_cont_old,S_old,ds_old,sh_old,deathmat,sh_val_h_old,sh_val_f_old,cost_vec_old,succ_prob_old)
%this function simply puts new firms into new cells
flag = 0;

st_cont = [];

new_cell_no = sum(cell2mat(deathmat));
if new_cell_no>5*S_old
    flag = 1;
    display('Too much death!');
end


% initialize new cell arrays (note old arrays called, well, old)
S           = S_old+new_cell_no; % new number of firms (after death is included)
st_ind_cont = cell(S,1);
ds          = cell(S,1);
sh          = cell(S,1); 
sh_val_h    = cell(S,1);
sh_val_f    = cell(S,1);
cost_vec    = cell(S,1);
succ_prob   = cell(S,1);

obin = S_old+1; 
if flag == 0
  for j = 1:S_old %loop through the preset number of firm slots
      ind = find(deathmat{j}==1); %find position of any death
      ind_ext = [ind;size(st_ind_cont_old{j},1)]; % append final slot of event matrix for firm j 

      % first firm in slot j
      sh_val_h{j} = sh_val_h_old{j}(1:ind_ext(1),:); %home shipment values
      sh_val_f{j} = sh_val_f_old{j}(1:ind_ext(1),:); %foreign shipment values
      cost_vec{j} = cost_vec_old{j}(1:ind_ext(1),:); %cost vector
      succ_prob{j} = succ_prob_old{j}(1:ind_ext(1),:); %cost vector
      st_ind_cont{j} = st_ind_cont_old{j}(1:ind_ext(1),:); %other firm related stuff
      ds{j} = ds_old{j}(1:ind_ext(1),:); %demand shock (client hotel)
      sh{j} = sh_old{j}(1:ind_ext(1),:); %shipments (client hotel)

      % loop through 2nd to nth firms in slot j
      lag = ind_ext(1)+1; %first event row of next firm in slot j
      for k = 2:size(ind_ext)
          sh_val_h{obin} = sh_val_h_old{j}(lag:ind_ext(k),:);
          sh_val_f{obin} = sh_val_f_old{j}(lag:ind_ext(k),:);
          cost_vec{obin} = cost_vec_old{j}(lag:ind_ext(k),:);
          succ_prob{obin} = succ_prob_old{j}(lag:ind_ext(k),:);
          st_ind_cont{obin} = st_ind_cont_old{j}(lag:ind_ext(k),:);
          ds{obin} = ds_old{j}(lag:ind_ext(k),:);
          sh{obin} = sh_old{j}(lag:ind_ext(k),:);
          lag = ind_ext(k)+1; %increment first even row
          obin = obin+1; %increment firm number in slot j 
      end
  end
end
end
