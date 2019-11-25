function [false_alarms,hits] = calcroc(output,LTE)
% [false_alarms,hits] = calcroc(output,LTE)

%[ld,idx] = sort(output); 
%i=0 ;

%hits= 1-cumsum(LTE(idx)>0)/sum(LTE > 0) ;
%false_alarms=1-cumsum(LTE(idx)<0)/sum(LTE < 0) ;

%if nargout<2,
%  false_alarms = [false_alarms;hits];
%end ;

%return 

% old

%ROC
[ld,ii] = sort(output); 
lte = LTE(ii);
lk = ones(1,length(lte));

hits=zeros(1,length(lte)) ;
false_alarms=zeros(1,length(lte)) ;
np=sum(lte > 0) ;
nn=sum(lte < 0) ;
for i=1:length(lte)
  lk(i) = -1;
  hits(i) = sum((lk > 0) & (lte >0))/np;
  false_alarms(i) = sum((lk > 0) & (lte < 0))/nn ;
end






