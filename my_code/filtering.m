bds = cell2mat(bds);
% static_FCemp_filt = ConSigEmp(:,:,1);
static_FCemp_filt = fcPath;
sim_ts = corrcoef(bds,'rows','pairwise');
     [b,a]       =       butter(1, [0.1]/((1/1.94)/2));
 
   % Compute low-pass filtered version of simulated fMRI
        for node = 1:68,
            simts_filt(:,node)    =   filtfilt(b,a,sim_ts(:,node));
        end
        
        Isubdiag = find(tril(ones(68),-1));
        
             % same for filtered data
        static_FC_sim_filt                  =   corrcoef(simts_filt,'rows','pairwise');        
        static_FC_cc_filt = corrcoef(static_FC_sim_filt(Isubdiag),static_FCemp_filt(Isubdiag));