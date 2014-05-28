function [st_cont,st_ind_cont,S,ds,sh,flag,sh_val_h,sh_val_f] = sdead(st_cont,st_ind_cont_old,S_old,ds_old,sh_old,deathmat,sh_val_h_old,sh_val_f_old)
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

obin = S_old+1;
if flag == 0
  for j = 1:S_old
      ind = find(deathmat{j}==1);
      ind_ext = [ind;size(st_ind_cont_old{j},1)];
      sh_val_h{j} = sh_val_h_old{j}(1:ind_ext(1),:);
      sh_val_f{j} = sh_val_f_old{j}(1:ind_ext(1),:);
      st_ind_cont{j} = st_ind_cont_old{j}(1:ind_ext(1),:);
      ds{j} = ds_old{j}(1:ind_ext(1),:);
      sh{j} = sh_old{j}(1:ind_ext(1),:);
      lag = ind_ext(1)+1;
      for k = 2:size(ind_ext)
          sh_val_h{obin} = sh_val_h_old{j}(lag:ind_ext(k),:);
          sh_val_f{obin} = sh_val_f_old{j}(lag:ind_ext(k),:);
          st_ind_cont{obin} = st_ind_cont_old{j}(lag:ind_ext(k),:);
          ds{obin} = ds_old{j}(lag:ind_ext(k),:);
          sh{obin} = sh_old{j}(lag:ind_ext(k),:);
          lag = ind_ext(k)+1;
          obin = obin+1;
      end
  end
end
end
