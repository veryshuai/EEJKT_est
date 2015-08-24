function [Data, W] = target_stats()
% returns the data statistics used in the loss function

    %% TARGETS (ALL REGRESSION COEFFICIENTS, CJ system of equations May 12 COMPLETE.pdf, from email from Jim on May 13)
    cli_coefsDAT = [0.04783 0.02771 0.23173 0.24793 0.01436 0.05980];
    exp_death_coefsDAT = [0.45274 0.00235 -0.22630 0.04174 -0.01253 -0.06800];
    match_death_coefsDAT = [0.79447 0.03421 -0.03163 -0.05370 -0.02778];
    exp_sales_coefsDAT = [9.80376 -1.06980 0.28000 0.43155];
    match_ar1_coefsDAT = [1.59164 0.82637 0.32834 0.06312];
    loglog_coefsDAT = [0.02116 -1.88130 -0.05446];
  % mavshipDAT = 10.59674540; % average in levels
    mavshipDAT = 0.9706402010; % average in logs
    exp_dom_coefsDAT = [1.502486,0.7268059]; %regression coefficients in Marcela's email of 3-22-2013, eam_moms_out( DIAN only 21 Feb 2013).log
    dom_ar1_coefsDAT = [0.3380132,0.9764422]; %regression coefficients in Marcela's email of 3-22-2013, eam_moms_out( DIAN only 21 Feb 2013).log 

    Data = [cli_coefsDAT,exp_death_coefsDAT,match_death_coefsDAT,exp_sales_coefsDAT,match_ar1_coefsDAT,loglog_coefsDAT,mavshipDAT,exp_dom_coefsDAT,dom_ar1_coefsDAT];
    
    %% covariance matrices 

cli_coefsCOV = ...   
   [0.001701382136  -0.001647552086 -0.000264459058  0.000073487156 -0.000959542049 -0.000038182353;
    -.001647552086   0.001668967452  0.000241358631 -0.000063003478  0.001042437496 -0.000022979405;
    -.000264459058   0.000241358631  0.000314554178 -0.000132500191  0.000040625185 -0.000012385983;
    0.000073487156  -0.000063003478 -0.000132500191  0.000073241876 -0.000015391800  0.000000882907;
    -.000959542049   0.001042437496  0.000040625185 -0.000015391800  0.000979801955 -0.000093769282;
    -.000038182353  -0.000022979405 -0.000012385983  0.000000882907 -0.000093769282  0.000069673070];


exp_death_coefsCOV = ...
   [0.002111871537  -0.002045605521  -0.000321397341  0.000089417063 -0.001208558470 -0.000047812999;
   -0.002045605521   0.002072826086   0.000294070954 -0.000076753577  0.001312874423 -0.000029401386;
   -0.000321397341   0.000294070954   0.000383671241 -0.000161878633  0.000050920706 -0.000016052533;
    0.000089417063  -0.000076753577  -0.000161878633  0.000089679064 -0.000019381918  0.000001153981;
   -0.001208558470   0.001312874423   0.000050920706 -0.000019381918  0.001246198277 -0.000118914785;
   -0.000047812999  -0.000029401386  -0.000016052533  0.000001153981 -0.000118914785  0.000088682676];
 
match_death_coefsCOV = ...    
    [0.000427699832  -0.000127448542  -0.000026253311  -0.000045560267 -0.000033872676; 
    -0.000127448542   0.000137342677   0.000001753036   0.000081273853 -0.000005694758; 
    -0.000026253311   0.000001753036   0.000002501209  -0.000000672050 -0.000000605574;
    -0.000045560267   0.000081273853  -0.000000672050   0.000081302304 -0.000027033495;
    -0.000033872676  -0.000005694758  -0.000000605574  -0.000027033495  0.000042624529];
    
exp_sales_coefsCOV = ...
    [0.005961707962 -0.005792913998 -0.001678065759 -0.000174738189; 
    -0.005792913998  0.005773713389  0.001647870126  0.000090992680;
    -0.001678065759  0.001647870126  0.000967730696 -0.000082167676; 
    -0.000174738189  0.000090992680 -0.000082167676  0.000147385590];
 
match_ar1_coefsCOV = ...
     [0.002362705419 -0.000156827604 -0.000481368001 -0.000284962242;
     -0.000156827604  0.000014721579  0.000014300618 -0.000004772804;
     -0.000481368001  0.000014300618  0.000332671396  0.000124115947;
     -0.000284962242 -0.000004772804  0.000124115947  0.000193084788];

loglog_coefsCOV = ...
  [0.022689638548  -0.015310017854   0.002420753146;
  -0.015310017854   0.012614810554  -0.002268716169;             
   0.002420753146  -0.002268716169   0.000443232645];            

    mavshipCOV = 0.00415553^2;  % note that I squared the standard error here

    exp_dom_coefsCOV = ...
        [0.1014345 0;
        0 0.0064868];

    dom_ar1_coefsCOV = ...
        [0.0126289 0;
        0 0.0008622];

 %   W = blkdiag(cli_coefsCOV,exp_death_coefsCOV,match_death_coefsCOV,exp_sales_coefsCOV,match_ar1_coefsCOV,loglog_coefsCOV,mavshipCOV,exp_dom_coefsCOV,dom_ar1_coefsCOV);
    
    % The following block is an ad hoc attempt to put the regressions on a more 
    % equal footing in terms of their weight in the fit metric.
    
    ncli_coefsCOV         = cli_coefsCOV./trace(cli_coefsCOV);
    nexp_death_coefsCOV   = exp_death_coefsCOV./trace(exp_death_coefsCOV);
    nmatch_death_coefsCOV = match_death_coefsCOV./trace(match_death_coefsCOV);
    nexp_sales_coefsCOV   = exp_sales_coefsCOV./trace(exp_sales_coefsCOV);
    nmatch_ar1_coefsCOV   = match_ar1_coefsCOV./trace(match_ar1_coefsCOV);
    nloglog_coefsCOV      = loglog_coefsCOV./trace(loglog_coefsCOV);
    nmavshipCOV           = mavshipCOV./trace(mavshipCOV);
    nexp_dom_coefsCOV     = exp_dom_coefsCOV./trace(exp_dom_coefsCOV);
    ndom_ar1_coefsCOV     =  dom_ar1_coefsCOV./trace( dom_ar1_coefsCOV);
    
   W = blkdiag(ncli_coefsCOV,nexp_death_coefsCOV,nmatch_death_coefsCOV,nexp_sales_coefsCOV,...
       nmatch_ar1_coefsCOV,nloglog_coefsCOV,nmavshipCOV,nexp_dom_coefsCOV,ndom_ar1_coefsCOV);
    
    
end
