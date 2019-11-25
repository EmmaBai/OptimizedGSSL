
function [ROCscore,CpuTime,gBeta] = starter(C,Y)

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Information
%--------------------------------------------------------------------
% data loading.... 
% "L{.}": Normalized Laplacian
% "W{.}": Similarity Matrix,
% "yMat"   : output
%

fileHead  = 'measurement';

% Create result file in res/ subfolder
rltFileName = sprintf('res/%s_res.txt', fileHead);

% Load the multi_biograph.mat file in the root
load( sprintf('%s.mat', fileHead) ); 

CV_fold = 10;  


% K is the number of datasets
K = 3;    
const = 0.7;
ratio  = 0.4;

  
  C0  = ratio*C; 

  
  % Check for Class Imbalance
  NData  = length(Y);
  NClass = hist(Y,2);
  info=sprintf('\n\n Class( 1): %d\n Class(-1): %d\n', NClass(2), NClass(1));
  fprintf('%s',info);    


  IndexC1 = find(Y==-1);
  IndexC2 = setdiff([1:NData]', IndexC1);
 
  LengthOfAfoldC1 = floor(length(IndexC1)/CV_fold);
  LengthOfAfoldC2 = floor(length(IndexC2)/CV_fold);

  % For each k-fold subsamples
  for cv=1:CV_fold    
    

    TeIndex = [IndexC1(LengthOfAfoldC1*(cv-1)+1 : LengthOfAfoldC1*cv);...
      IndexC2(LengthOfAfoldC2*(cv-1)+1 : LengthOfAfoldC2*cv)];
    % Get the training indexes
    TrIndex = setdiff([1:NData]', TeIndex);
    

    beta  = zeros(K,1);

    cputime1  = cputime;
   
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
     





