% This code snippet fills in the client hotel, so that each row of a current relationship contains the demand shock (until this point, only the first row of a relationship contains the demand shock

%so far we have just recorded the demand shock in the period in which they change.  To
%calculate profits, we need to fill in the ds matrix.
%home
for t = 1:maxc
     loop_ind = find((ds(:,end-1)) == t);
      if isempty(loop_ind) == 0 
        for o = 1:size(loop_ind,1)
            k = loop_ind(o);
        rel_time = ind(k,11); %time before exogenous separation
        p = k+1;
        while p<size(ind,1) && (ind(p,1)-ind(k,1)<rel_time-1e-12) %note: I had to take a little bit off here to deal with rounding error
            if ds(p,t) == 0; 
                temp = ds(p-1,t);
                ds(p,t) = temp;    
            end
            p = p+1;
        end
        end
      end
end
%foreign
for t = 1:maxc
    loop_ind = find((ds(:,end)) == maxc+t);
    if isempty(loop_ind) == 0
        for o = 1:size(loop_ind,1)
            k = loop_ind(o);
        rel_time = ind(k,12); %time before exogenous separation
        p = k+1;
        while p<size(ind,1) && (ind(p,1)-ind(k,1)<rel_time-1e-12) %note: I had to take a little bit off here to deal with rounding error
            if ds(p,maxc+t) == 0; 
                temp = ds(p-1,maxc+t); 
                ds(p,maxc+t)=temp;    
            end
            p = p+1;
        end
        end
    end
end

