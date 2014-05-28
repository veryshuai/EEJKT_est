function [avdur,qtile] = quant_dur(dur,initex,tiles)

%% 
cutoff    = 0:1/tiles:1;
cutoff    = 100*cutoff;
newmatch  = find(initex>0);
[nr,nc]   = size(newmatch);
sizedis   = reshape(initex(newmatch),nr*nc,1);
qtile     = prctile(sizedis,cutoff);

%% 
dur2   = dur(newmatch);
init2  = initex(newmatch);
avdur  = zeros(1,tiles-1);
  for qq = 1:tiles
    qqindx    = init2 <= qtile(qq+1) & init2 > qtile(qq);
    tot       = sum(sum(qqindx));
    avdur(qq) = sum(sum(dur2.*qqindx))/tot;
  end

end