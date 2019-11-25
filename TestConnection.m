function FlagT = TestConnection(TrIndex,TeIndex,W)

Index = union(TrIndex,TeIndex);
FlagT = TrIndex;
Ltr = length(TrIndex);
Lte = length(TeIndex);

temp = zeros(Lte,1);

for i = 1:Ltr
    for j=1:Lte
        if W(TrIndex(i),TeIndex(j))>0
            temp(j)=1;
        end
    end
end
 

if sum(temp)>0 && sum(temp)<Lte
    temp_TeIndex = find(temp>0);
    temp_TeIndex = TeIndex(temp_TeIndex);
    UTrIndex = union(FlagT,temp_TeIndex);
    UTeIndex = setdiff(Index,UTrIndex);
    FlagT = TestConnection(UTrIndex,UTeIndex,W);
elseif sum(temp)==0
    FlagT = TrIndex;
elseif sum(temp)==Lte
    temp_TeIndex = find(temp>0);
    temp_TeIndex = TeIndex(temp_TeIndex);
    FlagT = union(FlagT,temp_TeIndex);
end


