function epsilon = DscBandW(X,TrIndex,TeIndex)

X_L = X(TrIndex,:);
X_U = X(TeIndex,:);


nl = length(TrIndex);
nu = length(TeIndex);


%% calculate and sort candidate values
Candi = zeros(nu,nl);

for i = 1:nu
    for j = 1:nl
        Candi(i,j)=sum((X_U(i,:)-X_L(j,:)).^2);
    end
end

 


CandidateValues = unique(sort(min(Candi,[],2)));


epsilon = Can2epsilon(X,CandidateValues,TrIndex,TeIndex);


