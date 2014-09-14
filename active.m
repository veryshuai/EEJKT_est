%This snippet creates a dummy for active relationships

%active years for each client relationship (NOTE: if a relationship is ongoing,
%we count the client EVEN IF NO SALES OCCUR during the particular
%year.)

loop_ind = find((ind(1:obin,11))~=0);
if isempty(loop_ind)==0
    for o = 1:size(loop_ind,1)
        k = loop_ind(o);
        act(k,1) = max(ceil(ind(k,1)+ind(k,11)),burn)-max(floor(ind(k,1)),burn); %number of active years home
    end
end
loop_ind = find((ind(1:obin,12))~=0);
if isempty(loop_ind)==0
    for o = 1:size(loop_ind,1)
        k = loop_ind(o);
        act(k,2) = max(ceil(ind(k,1)+ind(k,12)),burn)-max(floor(ind(k,1)),burn); %number of active years foreign
    end
end
