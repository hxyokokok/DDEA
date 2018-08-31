function DDEA(Global)
% <algorithm> <A-G>
% dynamic decomposition based evolution algorithm 

%--------------------------------------------------------------------------
% The copyright of the PlatEMO belongs to the BIMK Group. You are free to
% use the PlatEMO for research purposes. All publications which use this
% platform or any code in the platform should acknowledge the use of
% "PlatEMO" and reference "Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu
% Jin, PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective
% Optimization, 2016".
%--------------------------------------------------------------------------

% Copyright (c) 2016-2017 BIMK Group
    
	w = ones(Global.M) * 1e6; % weights for ASF function
	w(logical(eye(Global.M))) = 1;
	N = Global.N;
	M = Global.M;
    Population = Global.Initialization();
    %% Optimization
    fit_ = zeros(N,1);
    while Global.NotTermination(Population)
    	MatingPool = binaryTournamentSelection(-fit_,N);
        Offspring  = Global.Variation(Population(MatingPool),N);
        All = [Population Offspring];
	    nFit = All.objs;
        % normalize-1
	    Z_star = min(nFit,[],1)-1e-6;
	    Z_nad = max(nFit,[],1);
	    nFit = bsxfun(@rdivide,bsxfun(@minus,nFit,Z_star),Z_nad - Z_star);
% 	    nFit(:,Z_nad-Z_star<1e-10) = 1e-10;
%         flag = sum(nFit,2)<1e-10;
%         nFit(flag,:) = 1e-10;
        % sort by 1-norm
	    [~,sortedIdx] = sort(sum(nFit,2));
	    All = All(sortedIdx);
	    nFit = nFit(sortedIdx,:);
        % construct the set of reference points
		Lambda = bsxfun(@rdivide,nFit,sum(nFit,2));
		Lambda(Lambda<1e-6) = 1e-6;
	    % find extreme points
	    extremeIdx = [];
	    for i = 1 : M
	    	[~,e_idx] = min(max(bsxfun(@times,nFit,w(i,:)),[],2));
	    	extremeIdx = [extremeIdx e_idx];
	    end
	    extremeIdx = unique(extremeIdx);
		% move the extreme points to the next population
		candidateIdx = 1 : 2*N;
		candidateIdx(extremeIdx) = [];
		nextIdx = extremeIdx;
		% distance matrix
		DM = pdist2(Lambda,Lambda);
		for j = 1 : 2*N
			DM(j,j) = 0;
		end
	    nearestDst2nextPop = min(DM(candidateIdx,nextIdx),[],2);

	    % aggregation matrix
		AM = zeros(2*N,2*N);
	    for j = 1 : M
	    	AM = max(AM, bsxfun(@rdivide,nFit(:,j),Lambda(:,j)'));
	    end
		while length(nextIdx) < N
			% select a solution d distant to the existed solution
			[~,d_idx] = max(nearestDst2nextPop);
			d_idx = candidateIdx(d_idx);
			% find the associated solutions
			associatedIdx = candidateIdx(DM(candidateIdx,d_idx)<=nearestDst2nextPop);
			[~,s_idx] = min(AM(associatedIdx,d_idx));
			s_idx = associatedIdx(s_idx);
			% update the nearest distance to the next population
			nearestDst2nextPop = min(DM(candidateIdx,s_idx),nearestDst2nextPop);
			nearestDst2nextPop(candidateIdx==s_idx) = [];
			% move this solution to the next population
			candidateIdx(candidateIdx==s_idx) = [];
			nextIdx = [nextIdx s_idx];
		end
        Population = All(nextIdx);
        fit_ = [zeros(1,length(extremeIdx)) 1:(N-length(extremeIdx))]';
    end
end