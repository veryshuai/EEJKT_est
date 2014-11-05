% This inserted script returns policy function

% First half policy 

case_str = 'non';

% Get parameters
SetParams_noprod;

%reset random seed
if seed == 1
    rng(80085);
end

% Get policy and value functions
[lambda_f_orig,lambda_h_orig,pi_tilda_h_orig,pi_tilda_f_orig,c_val_h_orig,c_val_f_orig,punishment_orig] = solve(mm);

% Second half policy

switch cf_num
    case 0
        case_str = 'est';
        increase = 1; % no counterfactual

        % Get policy and value functions
        lambda_f_new   =  lambda_f_orig;
        lambda_h_new   =  lambda_h_orig;
        pi_tilda_h_new =  pi_tilda_h_orig;
        pi_tilda_f_new =  pi_tilda_f_orig;
        c_val_h_new    =  c_val_h_orig;
        c_val_f_new    =  c_val_f_orig;
        punishment_new =  punishment_orig;

    case 1
        case_str = 'sim';
        increase = 1; % no counterfactual

        % Get parameters (firm number)
        SetParams_noprod;

        % Get policy and value functions
        lambda_f_new   =  lambda_f_orig;
        lambda_h_new   =  lambda_h_orig;
        pi_tilda_h_new =  pi_tilda_h_orig;
        pi_tilda_f_new =  pi_tilda_f_orig;
        c_val_h_new    =  c_val_h_orig;
        c_val_f_new    =  c_val_f_orig;
        punishment_new =  punishment_orig;

    case 2
      %cost_dec_trans(1);

    case 3

        case_str = 'mac';

        % Set counterfactual
        increase = 1.2; % 20% jump in macro shock
        scale_f    =  scale_f + log(increase);

        % Get parameters
        SetParams_noprod;
        
        %reset random seed
        if seed
            rng(80085);
        end
    
        % Get policy and value functions
        [lambda_f_new,lambda_h_new,pi_tilda_h_new,pi_tilda_f_new,c_val_h_new,c_val_f_new,punishment_new] = solve(mm);

    case 4
      %red_var_trans(1);

    case 5
      %search_dec_trans(1);

    case 6 % calculate value of network

        case_str = 'val';
        increase = 1; % no counterfactual

        % Get parameters
        SetParams_noprod;

        % Get policy and value functions
        lambda_f_new   =  lambda_f_orig;
        lambda_h_new   =  lambda_h_orig;
        pi_tilda_h_new =  pi_tilda_h_orig;
        pi_tilda_f_new =  pi_tilda_f_orig;
        c_val_h_new    =  c_val_h_orig;
        c_val_f_new    =  c_val_f_orig;
        punishment_new =  punishment_orig;

    case 7 % plot policies 

        increase = 1; % no counterfactual
        run_policy_stuff;

    otherwise

        case_str = 'non';
        increase = 1; % no counterfactual

        % Get policy and value functions
        lambda_f_new   =  lambda_f_orig;
        lambda_h_new   =  lambda_h_orig;
        pi_tilda_h_new =  pi_tilda_h_orig;
        pi_tilda_f_new =  pi_tilda_f_orig;
        c_val_h_new    =  c_val_h_orig;
        c_val_f_new    =  c_val_f_orig;
        punishment_new =  punishment_orig;

end

% Make punishment variable
try
    punishment = punishment_orig + punishment_new;
catch e
    punishment = punishment_orig;
end

