clear all
load('./Connectivity_significant_matrix.mat');
C = ConSig(:,:,1);
U=unique(ConSig(:,:,1));
R=sort(randn(1,length(C)));
Z=zeros(105,105);

% for i=1:1:105*105
%     
%     minMatrixC = min(C(:));
%     [rowC,colC] = find(C==minMatrixC);
%    
%     minMatrixR = min(R(:));
%     [rowR,colR] = find(R==minMatrixR);
%     
%     Z(rowC,colC)=R(rowR,colR);
%     C(rowC,colC)=10^5;
%     R(rowR,colR)=10^5;
%     Z(rowC,colC);
%     
% end 



for i=1:length(U)
C(C==U(i))=R(i); 
end 

X = (C-min(min(C)))/(max(max(C))-min(min(C)))*0.2+0.4;
