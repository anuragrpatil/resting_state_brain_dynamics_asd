function [corr maxfrNs FC bds]= optimal_G(C,C_emp,simTime,dt,G,noiseAmp)


corr = zeros(1,length(G));
maxfrNs = cell(1,length(G));
FC = cell(1,length(G));
bds= cell(1,length(G));

parfor g=1:length(G)
    
   
  [corr(1,g) maxfrNs{g}  ]=DMF_main(C,C_emp,simTime,dt,G(g),noiseAmp); 
  
  
end

end