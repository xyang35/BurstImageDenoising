function sigma2 = computeSigma2FromDiffVector(diffVector)
    %sigma2 = mean(diffVector.^2);
    sigma2 = var(diffVector);