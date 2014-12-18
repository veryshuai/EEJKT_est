function resp = query(quest,ans)
    %This function asks the user a request, checks for allowable responses, and returns the users response

    % Query user for desired task
    acceptable = 'FALSE'; %Is user input interpretable?
    while acceptable == 'FALSE'
        resp = input(quest, 's');
    
        display([char(10) 'Your input was ' num2str(resp) '.' char(10)]);
    
        %Check for validity
        for k = 1:size(ans,2)
            if resp == ans{k}
                acceptable = 'TRUE ';
            end
        end
        if acceptable == 'FALSE'
            display(['Sorry, I do not understand.  Try again.' char(10)]);
        end
    end
