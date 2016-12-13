function weights = GetVertexWeights(vw_fun,evals,evecs,timeVals,siHksDesc,siHksNorm)

assert(ischar(vw_fun),'GetVertexWeights:BadVwFunType','first input must be a string')

switch vw_fun
    case 'CT'
        RemoveZeroEval();
        weights = evecs.^2 * (1./sqrt(evals(:)));% commute-time
        
    case 'CT_3D'
        RemoveZeroEval();
        weights = evecs.^2 * (1./sqrt(evals(:).^3));% commute-time for volumes

    case 'heat_kernel'
        RemoveZeroEval();
        weights = evecs.^2 * exp(-evals(:) * timeVals(1));% diffusion
        
    case 'siHksNorm'
        weights = siHksNorm;%
        
    case 'siHks_at0'
        weights = siHksDesc(:,1);%
        
    otherwise
        error('GetVertexWeights:BadVwFunName','unknown VW fun name "%s"',vw_fun)
end
%%
    function RemoveZeroEval()
        isZeroVal = evals < eps(max(evals));
        evals = evals(~isZeroVal);
        evecs = evecs(:,~isZeroVal);
    end

end