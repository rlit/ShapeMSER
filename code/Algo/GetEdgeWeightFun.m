function ewFun = GetEdgeWeightFun(distFun,shape,evals,evecs,timeVals,siHksDesc,siHksNorm)

assert(ischar(distFun),'GetVertexWeights:BadVwFunType','first input musat be a string')


switch distFun
    case 'SIHKS_at0'
        ewFun = @(x,y) (siHksDesc(x,1) - siHksDesc(y,1)).^2;
        
    case 'siHksNorm'
        ewFun = @(x,y) (siHksNorm(x)   - siHksNorm(y)).^2;
        
    case 'geodesic'
        shapePoints = @(x) [shape.X(x(:)) shape.Y(x(:)) shape.Y(x(:))];
        ewFun = @(x,y) sum((shapePoints(x) - shapePoints(y)).^2,2);
        
    otherwise
        ewFun = @(x,y) HeatKernelDistance(x,y);
        
end

%%
    function distRes = HeatKernelDistance(v1,v2)
        nDist = numel(v1,1);
        distRes = zeros(nDist,1);
        
        hks = @(x,y,t,evals) (x.*y) * exp(-evals(:) * t(:)');
        
        switch distFun
            case 'scalar_HK'
                RemoveZeroEval()
                eVecDiff = abs(evecs(v1,:).^2 - evecs(v2,:).^2);
                eValsFactor = exp(- evals(:) * timeVals(1));
                distRes = eVecDiff * eValsFactor;
                
            case 'scalar_CT'
                RemoveZeroEval()
                eVecDiff = abs(evecs(v1,:).^2 - evecs(v2,:).^2);
                distRes = eVecDiff * (1./sqrt(evals(:)));% commute-time
                
            case 'time_integral'
                RemoveZeroEval()
                valsSum   = meshgrid(evals);
                valsSum   = valsSum + valsSum';
                
                expVals = exp(-valsSum * timeVals(1)) - exp(-valsSum * timeVals(2));
                
                eValsFactor = expVals ./ valsSum;
                %
                eVecDiff = evecs(v1,:).^2 - evecs(v2,:).^2;
                for iDist = 1:nDist
                    currVecDiff = eVecDiff(iDist,:);
                    distRes(iDist) = currVecDiff * eValsFactor * currVecDiff';
                end
                
            case 'diffusion_distance'
                RemoveZeroEval()
                eVecDiff = (evecs(v1,:) - evecs(v2,:));
                distRes  = hks(eVecDiff,eVecDiff,2*timeVals(1),evals) ;
                
            case 'inverse_CT'
                RemoveZeroEval()
                c_xy = evecs(v1,:) .* evecs(v2,:) * (1./evals);
                distRes = 1 ./ c_xy;
                
            case 'exponent_CT'
                RemoveZeroEval()
                c_xy = evecs(v1,:) .* evecs(v2,:) * (1./evals);
                distRes = exp(-c_xy);
                
            case 'inverse_HKS'
                RemoveZeroEval()
                distRes = 1 ./ hks(evecs(v1,:),evecs(v2,:),timeVals(1),evals);
                
            case 'exponent_HKS'
                RemoveZeroEval()
                distRes = exp(-hks(evecs(v1,:),evecs(v2,:),timeVals(1),evals));
                
            case 'inverse_SIHKS_diff'
                hks1 = hks(evecs(v1,:),evecs(v1,:),timeVals,evals);
                hks2 = hks(evecs(v2,:),evecs(v2,:),timeVals,evals);
                lastNonZero = find(all(hks1 & hks2),1, 'last');
                hks1 = hks1(:,1:lastNonZero);
                hks2 = hks2(:,1:lastNonZero);
                
                % remove sacling by derivative of a log
                hksDiff = (log(hks1) - log(hks2));
                
                % remove shift by taking abs of ftt
                hksDiff = abs(fft(hksDiff,[],2));
                
                distRes = 1./sum(hksDiff.^2,2);
                
            case 'inverse_SIHKS'
                hksVals = hks(evecs(v1,:),evecs(v2,:),timeVals,evals);
                lastNonZero = find(all(hksVals),1, 'last');
                hksVals = hksVals(:,1:lastNonZero);
                
                % remove sacling by derivative of a log
                hksVals = diff(log(hksVals),1,2);
                
                % remove shift by taking abs of ftt
                hksVals = abs(fft(hksVals,[],2));
                
                distRes = 1 ./ sum(hksVals.^2,2);
                
                %%
                %     case 'inverse_HKS_norm'
                %         %%
                %         distRes = sqrt(...
                %             hks(evecs(v1,:),evecs(v1,:),timeVals(1)) .* ...
                %             hks(evecs(v2,:),evecs(v2,:),timeVals(1))) ./ ...
                %             hks(evecs(v1,:),evecs(v2,:),timeVals(1));
            otherwise
                error('GetEdgeWeightFun:BadVwFunName','unknown EW fun name "%s"',distFun)
        end
        
    end
%%
    function RemoveZeroEval()
        isZeroVal = evals < eps(max(evals));
        evals = evals(~isZeroVal);
        evecs = evecs(:,~isZeroVal);
    end
end