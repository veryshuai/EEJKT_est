function [new_g,new_err,vio] = step(g,size,err)
%this function advances the counter g by one, subject to being smaller than
%size

vio = 0;

if g < size
   new_g = g+1; 
   new_err = err;
else 
    new_g = g;
    new_err = err + 1;
    vio = 1;
end
end