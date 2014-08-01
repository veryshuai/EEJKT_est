function out = discount(in_vec,b)
% discounts and returns in_vec at annual rate b

out = zeros(size(in_vec,1),1); %initialize the out vector

for k=1:size(in_vec,1)
    out(k) = b^(k-1) * in_vec(k); %discount
end
