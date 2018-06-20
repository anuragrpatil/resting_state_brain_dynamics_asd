function [fcCorrs] = compute_optimal_G(sc_td, fc_td, startG,endG,incG,simTime,dt,noiseAmp,saveFigPath)
    fc_emp = fc_td;
%     nAreas = size(fc_emp,1);
    G = startG:incG:endG;
%     nIters = 1;
    fcCorrs = zeros(length(G),1);
    fcs = cell(length(G),1);
    maxfrNs = cell(length(G),1);
%     avgFC = zeros(nAreas,nAreas);
    parfor i = 1:length(G)
        Gcurr = G(i);
%         parfor j = 1:nIters
%             [fcs{j}, frN] = DMF_dev(scPath,fic,ffi,simTime,dt,G,noiseAmp,lesionAreas,ficWeights,setFICWeights);
%         end
         [fcs{i}, frN] = DMF_excitatory(sc_td, simTime,dt,Gcurr,noiseAmp);
%        [fcs{i}, frN]=DMF_dev(scPath,fic,ffi,simTime,dt,Gcurr,noiseAmp,[],[],false);
%         for k = 1:nIters
%            avgFC = avgFC + fcs{k}; 
%         end
%         avgFC = avgFC/nIters;
        maxfrNs{i} = max(frN);
        fcCorrs(i) = find_corr(fcs{i},fc_emp);
    end
    figure('Visible','off');
    plot(G,fcCorrs,'x');
    print('-djpeg',saveFigPath);
    close all;
    save('./results/optimalG_264_dt1_sigma0.001_8mins.mat','fcs', 'maxfrNs', 'fcCorrs');
    %    save('/results/68/optimalG_68_fic_dt1_sigma0.001_8mins_EI.mat','fcs','maxfrNs','fcCorrs');
end
