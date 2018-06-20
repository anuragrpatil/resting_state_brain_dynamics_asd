clear all

scPath = hdf5read('Human_68.hdf5','/C');
fcPath =  hdf5read('Human_68.hdf5','/CC');
% load('./Connectivity_significant_matrix.mat');
% C_read = ConSig;
% C_normalize = C_read(:,:,9)/norm(C_read(:,:,9));
% C = C_normalize(:,:);
%  [m nAreas] = size(C);
% C(1:nAreas+1:nAreas*nAreas) = 0;
% scPath = C;
% 
% C_emp_read = ConSigEmp;
% C_emp_read(isinf(C_emp_read))=1;
% C_emp = C_emp_read(:,:,9)/norm(C_emp_read(:,:,9));
% fcPath = C_emp;

startG = 0.5;
endG = 0.5;
incG = 0.5;
fic = false;
ffi = false;
simTime = 20*60*100;
dt = 0.1;
noiseAmp = 0.001;
saveFigPath = 'G_20mins_test.jpeg';
saveCorrsPath = 'hello_test.hdf5';
% p = parpool(4);
[fcCorrs] = compute_optimal_G(scPath,fcPath,startG,endG,incG,simTime,dt,noiseAmp,saveFigPath);
% p.delete();
% clear p;
rmCmd = strcat(['rm ',saveCorrsPath]);
if(~system(strcat(['test -e ',saveCorrsPath])))
    system(rmCmd);
end
h5create(saveCorrsPath,'/FCCorrelations',length(startG:incG:endG));
h5write(saveCorrsPath,'/FCCorrelations',fcCorrs);
