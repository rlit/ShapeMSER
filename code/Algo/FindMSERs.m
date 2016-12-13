function [mserIdxs, mserProps] = FindMSERs(ADJ, val, area, filters)

% Default filter values
if nargin < 4, filters = struct; end
if ~isfield(filters, 'MinArea'),      filters.MinArea =      -Inf; end
if ~isfield(filters, 'MaxArea'),      filters.MaxArea =      +Inf; end
if ~isfield(filters, 'MaxScore'),     filters.MaxScore =     +Inf; end
if ~isfield(filters, 'MinDiversity'), filters.MinDiversity = -Inf; end
if ~isfield(filters, 'MinPeakClimb'), filters.MinPeakClimb = +Inf; end
if ~isfield(filters, 'KeepStable'),   filters.KeepStable =  false; end

TT = ADJ - ADJ';
leafs = find(~any(TT>0,2));

[row,col] = find(TT==-1);
parentsIdxs = zeros(length(val),1);
parentsIdxs(row) = col;

nTimesClimbed = zeros(size(val),'single');
nTimesPeak    = zeros(size(val),'single');

% Detect MSERs and assign instability scores Q.
% Q = Inf - not an MSER.
Q = Inf(length(val),1);
parentMser = zeros(length(val), 1);
for k=1:length(leafs)
    
    % find branch
    branch = Climb(parentsIdxs, leafs(k));
    nTimesClimbed(branch) = nTimesClimbed(branch) + 1;
    
    % calc instability
    v = val(branch);
    a = area(branch);
    
    dv = abs( v(1:end-2) - v(3:end) );
    da = (    a(1:end-2) - a(3:end) ) ./ a(2:end-1);
    q  = da ./ dv;
    
    % find minima
    if length(q) < 3,  continue; end
    minimaLoc = 1 + find(...% the "+1" is cuz q is sampled from 2 to end-1
        q(2:end-1) < q(1:end-2) & ...
        q(2:end-1) < q(3:end));
    q = q(minimaLoc);
    minimaLoc = branch(minimaLoc);
    
    if isempty(minimaLoc), continue; end
    nTimesPeak(minimaLoc) = nTimesPeak(minimaLoc) + 1;
    Q(minimaLoc) = min(Q(minimaLoc), q(:));
    
    % set MSER hierarchy
    if length(q) > 1
        parentMser(minimaLoc(2:end)) = minimaLoc(1:end-1);
    end
    
end
peakToClimbRatio = nTimesPeak ./ nTimesClimbed;

% Extract msers
[q,mserIdxs] = sort(Q, 'ascend');
isMser = (q<Inf);
q = q(isMser);
mserIdxs = mserIdxs(isMser);
val = val(mserIdxs);
area = area(mserIdxs);
peakToClimbRatio = peakToClimbRatio(mserIdxs);

% Renumber parents
LUT = zeros(length(Q)+1,1);
LUT(mserIdxs+1) = 1:length(mserIdxs);
parentMser = LUT(parentMser(mserIdxs)+1);

% Filter MSERs
if nargin >= 4,
    
    % Scan bottom-up
    [v,idx] = sort(val, 'descend'); %#ok<ASGLU>
    for p = idx(:)',
        if ... % kill region if one of filters is "active"
                peakToClimbRatio(p) > filters.MinPeakClimb || ...
                area(p) < filters.MinArea || ...
                area(p) > filters.MaxArea || ...
                q(p)    > filters.MaxScore
            q(p) = Inf;
            continue% <-- no need to check overlap
        end
        
        u = parentMser(p);
        if u > 0
            % parent exist -> check overlap
            nonoverlap = (area(u) - area(p)) / area(u);
            if nonoverlap < filters.MinDiversity
                % overlab too big -> %kill less stable region, or child region
                if filters.KeepStable && q(p) < q(u)
                    q(u) = Inf;
                else
                    q(p) = Inf;
                end                
            end            
        end
    end
    
end

isMser                   = (q < Inf);
mserProps.instability    = q(               isMser);
mserProps.val            = val(             isMser);
mserProps.area           = area(            isMser);
mserProps.peakClimbRatio = peakToClimbRatio(isMser);
mserIdxs                 = mserIdxs(        isMser);
