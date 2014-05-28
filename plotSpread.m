function out = plotSpread(varargin)
%PLOTSPREAD plots distributions of points by spreading them around the y-axis
%
% SYNOPSIS: ph = plotSpread(data,binWidth,spreadFcn,xNames,showMM)
%           ph = plotSpread(ah,...)
%
% INPUT data: cell array of distributions or nDatapoints-by-mDistributions array
%       binWidth : (opt) width of bins that control which data points are
%           considered close enough to be spread. Default: 0.1
%       spreadFcn : (opt) cell array of length 2 with {name,param}
%           if name is 'lin', the spread goes linear with the number of
%             points inside the bin, until it reaches the maximum of 0.9 at
%             n==param.
%           if name is 'xp', the spread increases as 1-exp(log(0.9)*x).
%             param is empty
%           Default {'xp',[]}
%       xNames : (opt) cell array of length nDistributions containing x-tick names
%               (instead of the default '1,2,3')
%       showMM : (opt) if 1, mean and median are shown as red crosses and
%                green squares, respectively. Default: 0
%                2: only mean
%                3: only median
%                4: mean +/- standard error of the mean (no median)
%                5: mean +/- standard deviation (no median)
%       ah  : handles of axes into which to plot
%
% OUTPUT ph: plot handle
%
% REMARKS: plotSpread is useful for distributions with a small number of
%          data points. For larger amounts of data, distributionPlot is
%          more suited.
%
% EXAMPLE: data = {randn(25,1),randn(100,1),randn(300,1)};
%          figure,plotSpread(data,[],[],{'25 pts','100 pts','300 pts'})
%
% created with MATLAB ver.: 7.9.0.3470 (R2009b) on Mac OS X  Version: 10.5.7 Build: 9J61
%
% created by: jonas
% DATE: 11-Jul-2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def_binWidth = 0.1;
def_spreadFcn = {'xp',[]};
def_xNames = [];
def_showMM = false;

%% CHECK INPUT

% check for axes handle
if ~iscell(varargin{1}) && length(varargin{1}) == 1 && ...
        ishandle(varargin{1}) && strcmp(get(varargin{1},'Type'),'axes')
    ah = varargin{1};
    data = varargin{2};
    varargin(1:2) = [];
    newAx = false;
else
    ah = gca;
    data = varargin{1};
    varargin(1) = [];
    newAx = true;
end

% check data. If not cell, convert
if ~iscell(data)
    [nPoints,nData] = size(data);
    data = mat2cell(data,nPoints,ones(nData,1));
else
    % get nData
    data = data(:);
    nData = length(data);
    % make sure all are vectors
    badCol = ~cellfun(@isvector,data) & ~cellfun(@isempty,data);
    if any(badCol)
        nCols = cellfun(@(x)(size(x,2)),data(badCol));
        warning('PLOTSPREAD:AUTORESHAPE',...
            'Elements %s of the cell array are not vectors. They will be reshaped automatically',...
            num2str(find(badCol)'));
        data(badCol) = cellfun(@(x)(x(:)),data(badCol),'UniformOutput',false);
    end
end

if length(varargin) > 0 && ~isempty(varargin{1})
    binWidth = varargin{1};
else
    binWidth = def_binWidth;
end
if length(varargin) > 1 && ~isempty(varargin{2})
    spreadFcn = varargin{2};
else
    spreadFcn = def_spreadFcn;
end
if length(varargin) > 2 && ~isempty(varargin{3})
    xNames = varargin{3};
else
    xNames = def_xNames;
end
if length(varargin) > 3 && ~isempty(varargin{4})
    showMM = varargin{4};
else
    showMM = def_showMM;
end


%% TRANSFORM DATA
% Here, I try to estimate what the aspect ratio of the data is going to be
fh = figure('Visible','off');
goodData = find(~cellfun(@isempty,data));
plot([0.5;nData+0.5],[min(data{goodData(1)});max(data{goodData(end)})],'o');
aspectRatio = get(gca,'DataAspectRatio');
close(fh);

tFact = aspectRatio(2)/aspectRatio(1);

%% SPREAD POINTS
[m,md,sem,sd] = deal(nan(nData,1));
for iData = 1:nData
    currentData = data{iData}(:);
    goodIdx = isfinite(currentData);
    currentData(~goodIdx) = [];
    
    if ~isempty(currentData)
        
        % transform and sort
        currentData = currentData / tFact;
        %currentData = sort(currentData);
        
        % add x
        currentData = [ones(size(currentData))*iData,currentData];
        
        % step through the data in 0.1 increments. If there are multiple
        % entries, spread along x
        for y = min(currentData(:,2)):binWidth:max(currentData(:,2))
            % find values
            valIdx = find(currentData(:,2) >= y & currentData(:,2) < y+binWidth);
            nVal = length(valIdx);
            if nVal > 1
                % spread
                switch spreadFcn{1}
                    case 'xp'
                        spreadWidth = 0.9*(1-exp(log(0.9)*(nVal-1)));
                    case 'lin'
                        spreadWidth = 0.9*min(nVal-1,spreadFcn{2})/spreadFcn{2};
                end
                spreadDist = spreadWidth / (nVal - 1);
                if isEven(nVal)
                    offset = spreadDist / 2;
                else
                    offset = eps;
                end
                for v = 1:nVal
                    currentData(valIdx(v),1) = iData + offset;
                    % update offset
                    offset = offset - sign(offset) * spreadDist * v;
                end
            end
        end
        
        % update data
        currentData(:,2) = data{iData}(goodIdx);
        data{iData} = currentData;
        
        if showMM > 0
            m(iData) = nanmean(currentData(:,2));
            md(iData) = nanmedian(currentData(:,2));
            sd(iData) = nanstd(currentData(:,2));
            sem(iData) = sd(iData)/sqrt(sum(isfinite(currentData(:,2))));
        end
    end % test isempty
end


%% plot
set(ah,'NextPlot','add')
try
allData = cat(1,data{:});
catch me
    if strcmp(me.identifier,'MATLAB:catenate:dimensionMismatch')
        % attempt to save that if there are 'strange' empties
        emptyIdx = cellfun(@isempty,data);
        [data{emptyIdx}] = deal([]);
        allData = cat(1,data{:});
    end
end
ph = plot(ah,allData(:,1),allData(:,2),'.');

% if ~empty, use xNames
set(ah,'XTick',1:nData);
if ~isempty(xNames)
    set(ah,'XTickLabel',xNames)
end
% have plot start/end properly
xlim([0,nData+1])

% add mean/median
mh = [];mdh=[];
if showMM
    % plot mean, median. Mean is filled red circle, median is green square
    if any(showMM==[1,2])
        mh = plot(1:nData,m,'+r','Color','r','MarkerSize',12);
    end
    if any(showMM==[1,3])
        mdh = plot(1:nData,md,'sg','MarkerSize',12);
    end
    if showMM == 4
        mh = plot(1:nData,m,'+r','Color','r','MarkerSize',12);
        mdh = myErrorbar(1:nData,m,sem);
    end
    if showMM == 5
        mh = plot(1:nData,m,'+r','Color','r','MarkerSize',12);
        mdh = myErrorbar(1:nData,m,sd);
    end
        
end

if nargout > 0
    out = ph;
end


