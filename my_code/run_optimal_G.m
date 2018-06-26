clc
clear all 
close all

%%%%%%%%%%%%%%%%%%%load data%%%%%%%%%%%%%

load('./Connectivity_significant_matrix.mat');
C_read = ConSig;
C_emp_read = ConSigEmp;
C_emp_read(isinf(C_emp_read))=1;

% C_read = hdf5read('Human_68.hdf5','/C');
% C_emp_read = hdf5read('Human_68.hdf5','/CC');


numSubs = 1;                        %number of subjects
C_normalize = zeros(105,105,numSubs);
C= zeros(105,105);
C_emp= zeros(105,105);


% C= zeros(68,68);
% C_emp= zeros(68,68);


%%%%%%%%%%%%%%%%%%%%% Parameters that are varied%%%%%%%%%%%%%
%  G = 0:0.5:4;
%    G=0:0.05:1.5;
  G=0.5;
simTime = 20*60*1000; % in 1s of ms 
dt = 0.1;
% noiseAmp=0.001:0.001:0.009;
noiseAmp=0.001;
corr = cell(numSubs,1);
maxfrNsubs = cell(numSubs,1);
FC = cell(numSubs,1);
bds = cell(numSubs,1);

% parpool(20);
currSub = C_read(:,:,9);
 parfor no = 1:numSubs

%      C_normalize(:,:,no) = C_read(:,:,no)/norm(C_read(:,:,no));
%         C_normalize(:,:,no) = (currSub-mean(currSub(:)))/std(currSub(:));
     C_normalize(:,:,no) = normalize_data(C_read(:,:,9),1);
     C = C_normalize(:,:,no);

%     C=C_read(:,:,no);
    [m nAreas] = size(C);
    C(1:nAreas+1:nAreas*nAreas) = 0;
%     C=zscore(C);
    C = (C-min(min(C)))./(max(max(C))-min(min(C)))*1+0;
%     C_emp = normalize_data(C_emp_read(:,:,no),2);
      C_emp = C_emp_read(:,:,9);
      C_emp = (C_emp-min(min(C_emp)))./(max(max(C_emp))-min(min(C_emp)))*1+0;

   

    [corr{no} maxfrNsubs{no}  ]= optimal_noiseAmp(C,C_emp,simTime,dt,G,noiseAmp);
    


 end
 
%  for i=1:9
% maxFr(i)=max(maxfrNsubs{1, 1}{1, i});
% end
 
save('done_for_ASD_EIC','corr','maxfrNsubs');
 
   


    



