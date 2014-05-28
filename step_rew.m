function [new_g,new_err] = step_rew(g,size,err)
%this function advances the counter g by one, subject to being smaller than
%size

if g < size
   new_g = g+1; 
   new_err = err;
else 
    new_g = 1;
    new_err = err + 1;
end
end