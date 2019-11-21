function [nullDist,nullPhase] = generateCoherenceNullDist(nPerms,continuousData,pointProcessData,movingWin,params)
% GENERATECOHERENCENULLDIST generates a null distribution for spike-field
%   coherence.
%
%   This function is basically a wrapper for coherencycpt, looping nPerms
%   times over random selection of continuous and point process data in
%   order to generate a null distribution for use with cluster statistics.
%
%   Null distributions should have the same first two dimensions as the
%   coherogram. The third dimension is trials, and the fourth is each
%   permutation.
%
%   see also: COHGRAMCPT
%


% author: EHS::20160414


% checking trial dimensions.
if isequal(length(pointProcessData),size(continuousData,2))
    nTrials = size(continuousData,2);
else
    error('continuous and point process data have different trial dimensions')
end

% makaing sure to average over each random trial permutation.
params.trialave = 1;

display('running shuffle test to determine null distribution for spike-field coherence')
parfor prm = 1:nPerms
    
    % updating user every 50 permutations.
    if mod(prm,50)==0
        display(sprintf('      - calculated coherence for permutation %d out of %d.',prm,nPerms))
    end
    
    % permutation trials.
    randIdxCont = round((nTrials-1).*rand(nTrials,1) + 1);
    randIdxPt = round((nTrials-1).*rand(nTrials,1) + 1);
    
    % calculating coherence.
    [nullDist(:,:,prm),nullPhase(:,:,prm)] = cohgramcpt(continuousData(:,randIdxCont),pointProcessData(randIdxPt), movingWin, params);
    
end


return


