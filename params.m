% This file holds parameter values for use in running serach and learning model code

% Link between X parameter values and parameters used in the model
    % lnF         =  scale_h+log(X(1));
    % delta       =  X(2);
    % scale_h     =  X(3);
    % scale_f     =  scale_h + log(X(4));
    % beta        =  X(5);
    % a           =  X(6);
    % b           =  X(7);
    % L_z         =  X(8);
    % D_z         =  X(9);
    % L_b         =  X(10);
    % gam         =  X(11)*(1+beta)/beta;
    % cost scalar =  X(12);
    % sig p       =  X(13);
    
% for scripts that take a single parameter vector
X = [0.0338,0.267,11.344,0.512,0.087,0.716,3.161,0.532,0.087,8.36,0.298,111.334585731,0.650]; 

% for scripts that take multiple parameter vectors
pop = [...  
    0.0338,0.267,11.344,0.512,0.087,0.716,3.161,0.532,0.087,8.36,0.298,111.334585731,0.650;...
    0.0338,0.267,11.344,0.512,0.087,0.716,3.161,0.532,0.087,8.36,0.298,111.334585731,0.650;...
    0.0338,0.267,11.344,0.512,0.087,0.716,3.161,0.532,0.087,8.36,0.298,111.334585731,0.650;...
   ];
