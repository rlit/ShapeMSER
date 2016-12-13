function result = FloodFill(ADJ, f, p)
nVertices = length(f);
nFloods   = numel(p);

result = false(nVertices,nFloods);

if iscell(ADJ)
    adjLut = ADJ;
else
    adjLut = Adj2Lut(ADJ);
end

for iFlood = 1:nFloods
    currPoints = p(iFlood);
    currThreshold = f(currPoints);
    
    isInCurrFlood = false(nVertices,1);
    while ~isempty(currPoints)
        isInCurrFlood(currPoints) = true;
        
        % check for currPoints' adjLut
        currNeighbors = unique([adjLut{currPoints}]);  
        
        % see if currNeighbors should become a point in next iteration
        isNextPoints = ...
            f(currNeighbors) >= currThreshold & ...
            isInCurrFlood(currNeighbors) == 0;
        currPoints = currNeighbors(isNextPoints);
    end
    
    result(:,iFlood) = isInCurrFlood;
end

end