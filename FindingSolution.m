%
% (C) Copyright 2004.-, HyunJung (Helen) Shin (2004-12-16).
%
function [beta, Z] = FindingSolution(beta, C, C0, L, Y, TrIndex, TeIndex, const)
    
% Combined Laplacian
% L is a <1x5 cell> array. Each cell is a <3588x3588 double> Laplacian matrix
K   = size(L,2);
n   = size(L{1},1);
nTr = length(TrIndex);
Y(TeIndex)=0;
 
% Class balancing matrix with Different Costs
E = ones(n,1);
ClassBalance = ones(1,nTr);
% Find out-of-class indexes from training set
idxC1 = find(Y(TrIndex)==-1);
% Find in-class indexes from training set
idxC2 = setdiff([1:nTr]', idxC1);
% Proportions of out-of-class and in-class samples in the training set
p1= length(idxC1)/nTr;
p2= length(idxC2)/nTr;
p = max(p1,p2);
% Fill ClassBalance array TODO
ClassBalance(idxC1) = 1 + const*sign(p2-p1)*p;
ClassBalance(idxC2) = 1 + const*sign(p1-p2)*p;

%Fill E TODO
E(TrIndex) = ClassBalance;

% Build matrix E from its diagonal
E = sparse(diag(E));

% Weight each element of Y using E
Y = E*Y;

% Define options for the optimization subroutine
options = optimset('Display',...
                   'iter',...
                   'GradObj',...
                   'on',...
                   'LargeScale',...
                   'off',...
                   'MaxFunEvals',...
                   70);

% Perform optimization (Find minimum of constrained nonlinear multivariable function)               
% @dualfune is the dual function defined in dualfune.m, with arguments (L,Y,E)
beta = fmincon(@dualfune,... % fun
               zeros(K,1),...% x0 solution initial point
               ones(1,K),... % A matrix for inequality constraints
               C,... % b right side of inequality constraints
               [],...% Aeq
               [],...% beq
               ones(K,1)*0.05*C,...% lb lower-bound for solution
               ones(K,1)*C0,...%ub upper-bound for solution
               [],... %nonlcon
               options,...
               L,... % @fun argument 
               Y,...% @fun argument 
               E);% @fun argument 
           
CombinedL = E;

% For each data set (5 in total)
for k=1:K 
  % Add the weighted Laplacian matrix 
  CombinedL = CombinedL + beta(k)*L{k};
end 

% Z seems to be the score-vector f from the paper
Z=sparse(CombinedL\Y);

return
