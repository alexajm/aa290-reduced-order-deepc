function [sMEV, ciEV, cumE, resL2, nBasis] = HornParallelAnalysis(data, K)
% Function HornParallelAnalysis.m simulates a distribution of eigenvalues
% by resampling a set of random variables of the real data size, and
% compares the eigenvalues of the real data and the distribution of
% eigenvalues from simulation; then the number of retained basis factors is
% decided by keeping those who are bigger than 95% of simulated
% distribution of eigenvalues. This is to implement the parallel analysis
% approach proposed by Horn (1965) and developed by Ledesma et al. (2007).
% 
% Input argument: data - standardized real data 
%                      (subtract mean and divide by std for each variable);
%                 K    - number of resampling.
% Output argument: sMEV - mean eigenvalues from simulation;
%                  ciEV - confidence intervals of eigenvalues from
%                         simulation;
%                  cumE - cumulative energy content / total energy content;
%                  resL2 - L2 norm of residual for each eigenvalue of real
%                          data;
%                  nBasis - number of retained basis factors.
%
% [Version 0] written by Lanya T. Cai on June 22 2016. Contact
% lanyavikins[at]gmail.com if there's any question.

    dataSize   = size(data);
    nEntry     = dataSize(1);
    nVar       = dataSize(2);
    covMtx     = data' * data ./ (nEntry-1); % Compute the covariance matrix
    [~, S, ~]  = svd(covMtx);                % singular value decomposition
    lambda     = diag(S);                    % the eigenvalues of covMtx
    cumE       = zeros(nVar,1);              % the cumulative energy content
    resL2      = zeros(nVar,1);       % L2 norm of residuals when retaining 
                                  % certain number of principal components.
    for iEV = 1:nVar,
        cumE(iEV)  = sum(lambda(1:iEV))/sum(lambda);
        resL2(iEV) = sqrt(sum(sum((pcares(data,iEV)).^2)));
    end

    simLatentRoots = zeros(nVar,K);
    for iSim = 1:K,
        simData    = randn(dataSize);
        covMtx_sim = simData' * simData ./ (length(simData)-1); % Compute the covariance matrix of simData
        [~, simS, ~]     = svd(covMtx_sim);                     % singular value decomposition
        simLatentRoots(:,iSim) = diag(simS);                    % the eigenvalues of covMtx of simData
    end
    sMEV = mean(simLatentRoots,2); % means of eigenvalue distributions from simulation 
    stdLatentRoots = std(simLatentRoots,0,2);
    lowCI_LR       = zeros(nVar,1);
    highCI_LR      = zeros(nVar,1);
    for iEV = 1:nVar,
        avgLR = sMEV(iEV);
        stdLR = stdLatentRoots(iEV);
        lowCI_LR(iEV)  = norminv(0.025, avgLR, stdLR);
        highCI_LR(iEV) = norminv(0.975, avgLR, stdLR);
    end
    ciEV = [lowCI_LR highCI_LR]; % size nVar-by-2
    nBasis = sum(lambda>highCI_LR); % Determine number of basis factors to retain.
end % The end of function HornParallelAnalysis.m