clear all

%scPath = hdf5read('Human_68.hdf5','/C');
%fcPath =  hdf5read('Human_68.hdf5','/CC');
 load('./Connectivity_significant_matrix.mat');
 C_read = ConSig;
 C_normalize = C_read(:,:,9);
 C_normalize = (C_normalize - mean(C_normalize(:)))/std(C_normalize(:));
 C = C_normalize(:,:);
[m nAreas] = size(C);
C(1:nAreas+1:nAreas*nAreas) = 0;
scPath = C;
% 
 C_emp_read = ConSigEmp;
% 
 C_emp = C_emp_read(:,:,9);
% C_emp(isinf(C_emp))=1;
fcPath = C_emp;

startG = 0.5;
endG = .5;
incG = 0.5;
fic = false;
ffi = false;
simTime = 20*60*100;
dt = 0.1;
noiseAmp = 0.001;
saveFigPath = 'G_20mins_test.jpeg';
saveCorrsPath = 'hello_test.hdf5';
% p = parpool(4);
[fcCorrs maxfrNs bds] = compute_optimal_G(scPath,fcPath,startG,endG,incG,simTime,dt,noiseAmp,saveFigPath);
% p.delete();
% clear p;
rmCmd = strcat(['rm ',saveCorrsPath]);
if(~system(strcat(['test -e ',saveCorrsPath])))
    system(rmCmd);
end
%h5create(saveCorrsPath,'/FCCorrelations',length(startG:incG:endG));
%h5write(saveCorrsPath,'/FCCorrelations',fcCorrs);
