% This inserted script returns policy function

% First half policy 

case_str = 'non';


% Get parameters
SetParams_noprod;

%reset random seed
rng(80085);

% Get policy and value functions
[lambda_f_orig,lambda_h_orig,pi_tilda_h_orig,pi_tilda_f_orig,c_val_h_orig,c_val_f_orig,punishment_orig] = solve(mm);

% Second half policy

switch cf_num
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
        rng(80085);
    
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

        case_str = 'non';
        increase = 1; %no counterfactual
        policy_plot('baseline',lambda_f_orig,0,lambda_f_orig);
        return

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

