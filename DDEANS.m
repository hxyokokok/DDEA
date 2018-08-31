function DDEANS(Global)
% <algorithm> <A-G>
% DDEA + non-dominated sorting

%--------------------------------------------------------------------------
% The copyright of the PlatEMO belongs to the BIMK Group. You are free to
% use the PlatEMO for research purposes. All publications which use this
% platform or any code in the platform should acknowledge the use of
% "PlatEMO" and reference "Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu
% Jin, PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective
% Optimization, 2016".
%--------------------------------------------------------------------------

% Copyright (c) 2016-2017 BIMK Group

    Population = Global.Initialization();
    N = Global.N;
    FrontNo = zeros(1,N);
    rankings = zeros(1,N);
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,length(Population),FrontNo,rankings);
        Offspring  = Global.Variation(Population(MatingPool));
        All = [Population Offspring];
	    % Non-dominated sorting
	    [FrontNo,MaxFNo] = NDSort(All.objs,N);
	    % ranking 
	    rankings = inf(size(FrontNo));
	    for level = 1 : (MaxFNo-1)
		    rankings(FrontNo==level) = DDRanking(All(FrontNo==level).objs,sum(FrontNo==level));
	    end
	    rankings(FrontNo==MaxFNo) = DDRanking(All(FrontNo==MaxFNo).objs,N-sum(FrontNo<MaxFNo));
	    Next = ~isinf(rankings);
	    Population = All(Next);
        FrontNo = FrontNo(Next);
	    rankings = rankings(Next);
    end
end

% ranking the top K solutions in Fitness 
function rankings = DDRanking(Fitness,K)
    [N,M] = size(Fitness);
    w = ones(M) * 1e6; % weights for ASF function
	w(logical(eye(M))) = 1;
	rankings = inf(N,1);
    % normalize-1
    Z_star = min(Fitness,[],1)-1e-6;
    Z_nad = max(Fitness,[],1);
    nFit = bsxfun(@rdivide,bsxfun(@minus,Fitness,Z_star),Z_nad-Z_star);
    % nFit(:,Z_nad-Z_star<1e-10) = 1e-10;

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
    if K <= length(extremeIdx)
    	extremeIdx = extremeIdx(randperm(length(extremeIdx),K));
    	rankings(extremeIdx) = 0;
    	return;
    end
	% move the extreme points to the next population
	candidateIdx = 1:N;
	candidateIdx(extremeIdx) = [];
	nextIdx = extremeIdx;
	% distance matrix
	DM = pdist2(Lambda,Lambda);
	for j = 1 : N 
		DM(j,j) = 0;
	end
    nearestDst2nextPop = min(DM(candidateIdx,nextIdx),[],2);
    % aggregation matrix
	AM = zeros(N,N);
    for j = 1 : M
    	AM = max(AM, bsxfun(@rdivide,nFit(:,j),Lambda(:,j)'));
    end
	while length(nextIdx) < K
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
    rankings(nextIdx) = [zeros(1,length(extremeIdx)) 1:(K-length(extremeIdx))]';
end