%
% (C) Copyright 2004.-, HyunJung (Helen) Shin (2004-12-16).
%
function [ROCscore,CpuTime,gBeta] = starter

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Information
%--------------------------------------------------------------------
% data loading.... 
% "L{.}": Normalized Laplacian, 5 x [ 3588 x 3588 ]
% "W{.}": Similarity Matrix, 5 x [ 3588 x 3588 ]
% "yMat"   : output [3588 x 13]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: NOTE: W is not used, just the Laplacian matrices L
fileHead  = 'multi_biograph';

% Create result file in res/ subfolder
rltFileName = sprintf('res/%s_res.txt', fileHead);

% Load the multi_biograph.mat file in the root
load( sprintf('%s.mat', fileHead) ); 

CV_fold = 5;                 % The number of folds used for Cross-Validation 
KindOfClass = [1:13];        % Define the Output class 

% K is the number of datasets (5 here) 
K = 5;    
% const is used for class balancing within FindingSolution function
const = 0.7;
% ratio is used to calculate C0  = ratio*C
ratio  = 0.4;

% Loop over each class. "cls" is the class index for Output class
for cls=1: length(KindOfClass) % 
  
  switch cls  % C determined by cross-validation
    case {1, 2, 7, 8},      C=5.0; 
    case {5, 6, 9, 10},     C=10.0; 
    case {3, 4, 13},        C=25;  
    case {11},              C=100; 
    case {12},              C=2.5;
  end
  
  % TODO: check C0
  C0  = ratio*C; 

  % yMat is a <3588x13 doubles> (n.b nb of proteins = 3588)
  % yMat(x,y) equals 1 if prot. x is of class y, -1 else
  Y = yMat(:,KindOfClass(cls));
  
  % Check for Class Imbalance
  NData  = length(Y);
  NClass = hist(Y,2);
  info=sprintf('\n\n Class( 1): %d\n Class(-1): %d\n', NClass(2), NClass(1));
  fprintf('%s',info);    

  % Data Partition for Cross-Validation
  % find prot. indices that are not in the current cls class
  IndexC1 = find(Y==-1);
  % Set difference of two arrays
  % C = setdiff(A,B) returns the data in A that is not in B.
  IndexC2 = setdiff([1:NData]', IndexC1);
 
  % Get the number of elements whithin each of the k-fold subsamples
  LengthOfAfoldC1 = floor(length(IndexC1)/CV_fold);
  LengthOfAfoldC2 = floor(length(IndexC2)/CV_fold);

  % For each k-fold subsamples
  for cv=1:CV_fold    
    
    % Get the test indexes  
    TeIndex = [IndexC1(LengthOfAfoldC1*(cv-1)+1 : LengthOfAfoldC1*cv);...
      IndexC2(LengthOfAfoldC2*(cv-1)+1 : LengthOfAfoldC2*cv)];
    % Get the training indexes
    TrIndex = setdiff([1:NData]', TeIndex);
    
    % beta is output/input of FindingSolution function
    % it the weights vector of each dataset
    beta  = zeros(K,1);

    cputime1  = cputime;
    % Call of the main function. 
    % L is the graph laplacian matrix (read from multi_biograph.mat file) 
    [beta, Z] = FindingSolution(beta, C, C0, L, Y, TrIndex, TeIndex, const);
    cputime2  = cputime;
    
    %  PERFORMANCE EVALUATION 
    % Evaluate confusion matrices
    TrConfMat = ConfusionMatrix(Y(TrIndex), Z(TrIndex));
    TeConfMat = ConfusionMatrix(Y(TeIndex), Z(TeIndex));

    TrErr(cv) = (1-trace(TrConfMat)/length(TrIndex))*100; 
    TeErr(cv) = (1-trace(TeConfMat)/length(TeIndex))*100; 
      
    CpuTime(cv)= cputime2-cputime1;
    ROCscore(cv) = calcrocscore(Z(TeIndex)',Y(TeIndex)');      
    
    gZ{cv}    = Z;
    gBeta{cv} = beta;
    
    fprintf('\n\n\n%s: class=%d (%d)', fileHead, cls, cv);
    fprintf('\nTrErr: %-10.4f'    , TrErr(cv)      );
    fprintf('\nTeErr: %-10.4f'    , TeErr(cv)      );
    fprintf('\nROCscore: %-10.4f' , ROCscore(cv)   );
    fprintf('\ncputime: %-10.4f'  , CpuTime(cv)    );
    
  end %end of For(cv=1:CV_fold)

  %  SAVE-RESULT for one class
  save(sprintf('res/%s_F%d.mat', fileHead, cls), 'gZ','gBeta','TrErr','TeErr','ROCscore');
  
  %  REPORT in file
  Fp = fopen(rltFileName,'a');
      
  str='';
  tempBeta = [];
  for v=1:CV_fold
    str      = [str sprintf('%s', '%10.4f')]; 
    tempBeta = [tempBeta gBeta{v}];
  end 
    
  fprintf(Fp, '\n\n\n======%s  Class [%d] : %s ======', fileHead, cls, datestr(now));
  fprintf(Fp, '\n%s',info);
  fprintf(Fp, '\n%s ', sprintf('%-10s',' ')),        fprintf(Fp, '%10.0d', [1:CV_fold]); fprintf(Fp, '%10s', 'avg');
  fprintf(Fp, '\n%s:', sprintf('%-10s','TrErr')),    fprintf(Fp, str, TrErr );    fprintf(Fp, ' : %10.4f', mean(TrErr) );
  fprintf(Fp, '\n%s:', sprintf('%-10s','TeErr')),    fprintf(Fp, str, TeErr );    fprintf(Fp, ' : %10.4f', mean(TeErr) ); 
  fprintf(Fp, '\n%s:', sprintf('%-10s','ROCscore')), fprintf(Fp, str, ROCscore ); fprintf(Fp, ' : %10.4f', mean(ROCscore) ); 
  fprintf(Fp, '\n%s:', sprintf('%-10s','CPUtime')),  fprintf(Fp, str, CpuTime );  fprintf(Fp, ' : %10.4f', mean(CpuTime) ); 
  fprintf(Fp, '\n\n%s:\n', sprintf('%-10s','Beta')), str = [sprintf('%-11s','') str sprintf('%s',' : %10.4f\n')]; 
  fprintf(Fp, str, [tempBeta' ; mean(tempBeta')]);        
    
  fclose(Fp);  
     
end 




