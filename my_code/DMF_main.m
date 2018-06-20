function [correleation maxfr]= DMF_main(C,C_emp,simTime,dt,G,noiseAmp)
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Load file%%%%%%%%%%%%%%%%%%%%
% C = sc_td;

% load('./Connectivity_significant_matrix.mat');
% C_read = ConSig;
% C_normalize = C_read(:,:,9)/norm(C_read(:,:,9));
% C = C_normalize(:,:);
%  [m nAreas] = size(C);
% C(1:nAreas+1:nAreas*nAreas) = 0;
% 
% C_emp_read = ConSigEmp;
% C_emp_read(isinf(C_emp_read))=1;
% C_emp = C_emp_read(:,:,9)/norm(C_emp_read(:,:,9));

% load('./Connectivity_significant_matrix.mat');

% C_read = ConSig;
% C_normalize = zscore(C_read(:,:,9));
% C = C_normalize(:,:);
% [m nAreas] = size(C);
% C(1:nAreas+1:nAreas*nAreas) = 0;
% 
% C_emp_read = ConSigEmp;
% C_emp_read(isinf(C_emp_read))=1;
% C_emp = zscore(C_read(:,:,9));


% 
% C_read = hdf5read('Human_68.hdf5','/C');

% C = C_read;
%  [m nAreas] = size(C);
% C(1:nAreas+1:nAreas*nAreas) = 0;

% C_emp_read = hdf5read('Human_68.hdf5','/CC');
% C_emp = C_emp_read;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dt = 0.1;
    tend = simTime;
%     tend = 120000;
    T = 0:dt:tend;
%     G = 0.5;

S = zeros(tend*(1/dt),length(C),length(G));
S_persec = zeros(tend*dt,length(C),length(G));
% B = zeros(tend*dt,length(C),length(G)); 
% bds = zeros(tend*dt+2,length(C),length(G)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Compute DMF%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    [S H x] = DMF_eulers_explicit(T,dt,C,G,noiseAmp);



     S_persec = S(1:10:end,:);
     frN=H(1:10:end,:);
     maxfr=max(frN);
   
    % savetofile('Hvalues',H,G);
    %plot(T,S(:,1)');

    %bold signal transformation
    disp('synaptic activity has been calculated')
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Compute Bold%%%%%%%%%%%%%%%%%%%%

     nn = length(T(1,1:1/dt:end));  %nn = no.of 100 of milliseconds
     dtt = 0.001;
    T_bold = nn*dtt;  %no of seconds
    B = zeros(T_bold/dtt-1,length(C),length(G));

    parfor i =1:nAreas 

        B(:,i) = BOLD(T_bold,S_persec(:,i));

    end
% 


  %---------- Downsampling and reordering removing the first 500ms
    ds   = 2000; 
    bds = B(20500:ds:end,:);

    
    disp('Bold signal has been found')
    %%%%%%%%%%%%%%%%%%%Correleation coefficient%%%%%%%%%%%%%%%%%%%%%%%
    FC = corrcoef(bds);          % BOLD correlation matrix = Simulated Functional Connectivity
    

    FC_sim=atanh(FC(find(tril(FC,-1))));
    FC_sim(isinf(FC_sim))=1;     %replace infinte values by 1
    %correlation coeficient between simulated functional connectivity and
    %empirical functional connectivity
    FC_emp = atanh(C_emp(find(tril(C_emp,-1))));
    FC_emp(isinf(FC_emp))=1;
    correleation = corrcoef(FC_sim,FC_emp);
    
    correleation = correleation(1,2);
    
    %%%%%%%%%%%%%%%%%%%
%     maxH = max(max(H));
      
% save('workspace_june20_2');    
end

