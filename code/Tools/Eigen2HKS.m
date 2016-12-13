function desc = Eigen2HKS(evals,evecs,timeSamples)
expVal = exp(-bsxfun(@times, evals(:)          , timeSamples(:)'));
desc   =      bsxfun(@times,(evecs.^2) * expVal, timeSamples(:)');
desc   = normalize(desc, 'l2', 2);