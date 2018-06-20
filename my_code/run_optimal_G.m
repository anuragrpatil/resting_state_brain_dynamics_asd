clc
clear all 
close all

%%%%%%%%%%%%%%%%%%%load data%%%%%%%%%%%%%

% load('./Connectivity_significant_matrix.mat');
% C_read = ConSig;
% C_emp_read = ConSigEmp;

C_read = hdf5read('Human_68.hdf5','/C');
C_emp_read = hdf5read('Human_68.hdf5','/CC');


numSubs = 1;                        %number of subjects
% C_normalize = zeros(105,105,numSubs);
% C= zeros(105,105);
% C_emp= zeros(105,105);


C= zeros(68,68);
C_emp= zeros(68,68);


%%%%%%%%%%%%%%%%%%%%% Parameters that are varied%%%%%%%%%%%%%
G = 0.5:0.5:0.5;
simTime = 20*60*100; % in 10s of ms 
dt = 0.1;
noiseAmp=0.001;
corr = cell(length(numSubs),1);
maxfrNsubs = cell(length(numSubs),1);

% parpool(20);

 parfor no = 1:numSubs

%     C_normalize(:,:,no) = C_read(:,:,no)/norm(C_read(:,:,no));
%     C = C_normalize(:,:,no);
    C=C_read(:,:,no);
    [m nAreas] = size(C);
    C(1:nAreas+1:nAreas*nAreas) = 0;

    C_emp = C_emp_read(:,:,no);
    

    [corr{no} maxfrNsubs{no}]= optimal_G(C,C_emp,simTime,dt,G,noiseAmp)
    ;


 end
 
 save('results','corr','maxfrNsubs');
 
   


    



