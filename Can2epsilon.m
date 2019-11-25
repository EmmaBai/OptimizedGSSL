function epsilon_res = Can2epsilon(X,CandidateValues,TrIndex,TeIndex)

n_can = length(CandidateValues);

if n_can == 1
    epsilon_res = CandidateValues(1);
else
    n_median = ceil(0.5*n_can);
    epsilon = CandidateValues(n_median);
    
    W = KMsparse(X,epsilon);
    
    FlagT = TestConnection(TrIndex,TeIndex,W);
    
    Index = union(TrIndex,TeIndex);
    temp = setdiff(Index,FlagT);

    if isempty(temp)
        CandidateV = CandidateValues(1:n_median);
    else
        CandidateV = CandidateValues((n_median+1):n_can);
    end
    
    epsilon_res = Can2epsilon(X,CandidateV,TrIndex,TeIndex);
end
    