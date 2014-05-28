%topping off with Knitro

diary gradcalc.txt

%options = optimset('Display','iter','MaxIter',1);

load calib.mat

grad = zeros(14,1);

val = distance(X,shk);

for k = 1:14
    test = X;
    test(k) = X(k)+1e-8;
    grad(k) = (distance(test,shk)-val)/1e-8;
end

save grad

diary off