function [Data, W] = target_stats()
% returns the data statistics used in the loss function

    %% TARGETS (ALL REGRESSION COEFFICIENTS, CJ system of equations May 12 COMPLETE.pdf, from email from Jim on May 13)
    cli_coefsDAT = [0.04783 0.02771 0.23173 0.24793 0.01436 0.05980];
    exp_death_coefsDAT = [0.45274 0.00235 -0.22630 0.04174 -0.01253 -0.06800];
    match_death_coefsDAT = [0.79447 0.03421 -0.03163 -0.05370 -0.02778];
    exp_sales_coefsDAT = [9.80376 -1.06980 0.28000 0.43155];
    match_ar1_coefsDAT = [1.59164 0.82637 0.32834 0.06312];
    loglog_coefsDAT = [0.02116 -1.88130 -0.05446];
    mavshipDAT = [1.21480];
    exp_dom_coefsDAT = [1.502486,0.7268059]; %regression coefficients in Marcela's email of 3-22-2013, eam_moms_out( DIAN only 21 Feb 2013).log
    dom_ar1_coefsDAT = [0.3380132,0.9764422]; %regression coefficients in Marcela's email of 3-22-2013, eam_moms_out( DIAN only 21 Feb 2013).log 

    Data = [cli_coefsDAT,exp_death_coefsDAT,match_death_coefsDAT,exp_sales_coefsDAT,match_ar1_coefsDAT,loglog_coefsDAT,mavshipDAT,exp_dom_coefsDAT,dom_ar1_coefsDAT];
    
    %% covariance matrices 
    cli_coefsCOV = ...
        [0.001701  -0.001648   -0.00026  0.00007   -0.000960  -0.000038;
        -0.001648 0.001669    0.00024   -0.00006  0.001042   -0.000023;
        -0.000264 0.000241    0.00031   -0.00013  0.000041   -0.000012;
        0.000073  -0.000063   -0.00013  0.00007   -0.000015  0.000001;
        -0.000960 0.001042    0.00004   -0.00002  0.000980   -0.000094;
        -0.000038 -0.000023  -0.00001   0.00000  -0.000094   0.000070];

   exp_death_coefsCOV = ...
       [0.00211    -.002045606  -0.00032   0.000089   -0.001209    -0.000048;
        -0.00205   0.002072826  0.00029    -0.000077  0.001313     -0.000029;
        -0.00032   0.000294071  0.00038    -0.000162  0.000051     -0.000016;
        0.00009    -.000076754  -0.00016   0.000090   -0.000019    0.000001;
        -0.00121   0.001312874  0.00005    -0.000019  0.001246     -0.000119;
        -0.00005   -.000029401  -0.00002   0.000001   -0.000119    0.000089];


    match_death_coefsCOV = ...
        [0.00043  -0.000127  -0.000026   -0.000046   -0.000034;
        -0.00013 0.000137   0.000002    0.000081    -0.000006;
        -0.00003 0.000002   0.000003    -0.000001   -0.000001;
        -0.00005 0.000081   -0.000001   0.000081    -0.000027;
        -0.00003 -0.000006  -0.000001   -0.000027   0.000043];

    exp_sales_coefsCOV = ...
        [0.00596   -0.00579   -0.00168   -0.00017;
        -0.00579  0.00577    0.00165    0.00009;
        -0.00168  0.00165    0.00097    -0.00008;
        -0.00017  0.00009    -0.00008   0.00015];

    match_ar1_coefsCOV = ...  
        [0.00236   -0.00016  -0.00048   -0.000285;
        -0.00016  0.00001   0.00001    -0.000005;
        -0.00048  0.00001   0.00033    0.000124;
        -0.00028  -0.00000  0.00012    0.000193];

    loglog_coefsCOV = ...
        [0.022690  -0.01531  0.002421;
         -0.015310 0.01261   -0.002269;
         0.002421  -0.00227  0.000443];

    mavshipCOV = [2.27402];

    exp_dom_coefsCOV = ...
        [0.1014345 0;
        0 0.0064868];

    dom_ar1_coefsCOV = ...
        [0.0126289 0;
        0 0.0008622];

    W = blkdiag(cli_coefsCOV,exp_death_coefsCOV,match_death_coefsCOV,exp_sales_coefsCOV,match_ar1_coefsCOV,loglog_coefsCOV,mavshipCOV,exp_dom_coefsCOV,dom_ar1_coefsCOV);
    
end
