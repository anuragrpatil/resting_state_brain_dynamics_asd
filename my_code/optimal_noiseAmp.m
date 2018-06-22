function [corr maxfrNs FC bds]= optimal_noiseAmp(C,C_emp,simTime,dt,G,noiseAmp)


corr = cell(1,length(noiseAmp));
maxfrNs = cell(1,length(noiseAmp));
FC = cell(1,length(noiseAmp));
bds= cell(1,length(noiseAmp));

parfor na=1:length(noiseAmp)
    
   
  [corr{na} maxfrNs{na} FC{na} bds{na}]=optimal_G(C,C_emp,simTime,dt,G,noiseAmp(na)); 
  
  
end

end