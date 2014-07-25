% This script is the initial call into the search and learning project
% It gives the user the option of which task he would like to perform

acceptable = FALSE
while bool(acceptable = FALSE) 
    task = input('Task options are currently est (estimate the model) or sim (simulate the model once./n Which task shall I perform?: ')

    %Check for validity
    if task == 'est' | task == 'sim'
        acceptable = TRUE
    else
        display('Sorry, I do not understand.  Try again./n')
    end
end



