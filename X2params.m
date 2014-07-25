function X2params(X)
%This function takes an X vector passed by the optimizaiton program
%and turns it into parameters used in our simulation

    %read in parameter names from the X vector
    scale_h    =  X(3);
    ag         =  .5;
    bg         =  .5;
    lnF        =  scale_h + log(X(1));
    delta      =  X(2);
    scale_f    =  scale_h + log(X(4));
    beta       =  X(5);
    ah         =  X(6);
    bh         =  X(7);
    L_p        =  0;
    D_p        =  0;
    L_z        =  X(8);
    D_z        =  X(9);
    L_b        =  X(10);
    alp        =  0;
    gam        =  X(11)*(1+beta)/beta;
    cs         =  X(12);
    sig_p      =  X(13);

end
