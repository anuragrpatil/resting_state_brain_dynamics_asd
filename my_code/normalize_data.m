function X = normalize_data(C_read,flag)

% 
% % load('./Connectivity_significant_matrix.mat');
% % C = ConSig(:,:,1);
% % U=unique(ConSig(:,:,1));
% % R=sort(randn(1,length(C)));
% 
% C = C_read;
% U = unique(C_read);
% R=sort(randn(1,length(C)));
% 
% % for i=1:length(U)
% % C(C==U(i))=R(i); 
% % end
% 
% %  X = (C-min(min(C)))/(max(max(C))-min(min(C)))*0.2+0.4;
% %  X = zscore(C_read);
% %   X = C_read(:,:,no)/norm(C_read(:,:,no));
% 
% Z = reshape(C,[105*105 1]);
% if flag==1
% X = 0.5+ (C-mean(Z))*(0.1/std(Z));
% elseif flag==2
% X = 0+ (C-mean(Z))*1/std(Z);
load('Distances_between_centers_sig.mat')
Dist= Dist;
DistEmp=DistEmp; 

X=(C_read./Dist);
% X = zscore(X)

end


