[base,jac] = sim_solve(init,bet,a(:,:,k),pi,N+1-m,j,size(Q0,1),rh,diag_Q,Q0_d,V(:,N-m+2,:,k,cl),cl,gam);

siz = 1e-12;
jac_num = zeros(size(Q0,1));

for u = 1:size(Q0,1)
    test = init;
    test(u) = init(u) + siz;
    new_val = sim_solve(test,bet,a(:,:,k),pi,N+1-m,j,size(Q0,1),rh,diag_Q,Q0_d,V(:,N-m+2,:,k,cl),cl,gam);

    jac_num(:,u) = (new_val-base)/siz;
end

jac_dif = sparse(jac-jac_num);