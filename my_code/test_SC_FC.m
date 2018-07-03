load('./Connectivity_significant_matrix.mat');
C_read = ConSig;
C_emp_read = ConSigEmp;
C_emp_read(isinf(C_emp_read))=1;   
for i=1:68

%      C_normalize(:,:,i) = normalize_data(C_read(:,:,i),1);
%      C = C_normalize(:,:,i);
%          [m nAreas] = size(C);
%     C(1:nAreas+1:nAreas*nAreas) = 0;
%          x = (C-min(min(C)))./(max(max(C))-min(min(C)))*1+0;
%          
%                C_emp = C_emp_read(:,:,i);
%       y = (C_emp-min(min(C_emp)))./(max(max(C_emp))-min(min(C_emp)))*1+0;

     x = C_read(:,:,i);
     y = C_emp_read(:,:,i);
% 
%     FC_sim=atanh(FC(find(tril(FC,-1))));
% %     FC_sim=zscore(FC(find(tril(FC,-1))));
%     FC_sim(isinf(FC_sim))=1;     %replace infinte values by 1
%     
%     %correlation coeficient between simulated functional connectivity and
%     %empirical functional connectivity
%     FC_emp = atanh(C_emp(find(tril(C_emp,-1))));
% %     FC_emp = zscore(C_emp(find(tril(C_emp,-1))));
%     FC_emp(isinf(FC_emp))=1;
%     correleation = corrcoef(FC_sim(1:end),FC_emp);
%     
%     correleation(i) = correleation(1,2);
    
    nAreas = size(x,1);
        IsubDiag = find(tril(ones(nAreas),-1)); % indices of lower diagonal elements of SC matrix
%     xyCorr=corrcoef(x(IsubDiag),y(IsubDiag));
  xyCorr=corrcoef(zscore(x(IsubDiag)),zscore(y(IsubDiag)));
    Corr(i,1) = xyCorr(1,2);
end


