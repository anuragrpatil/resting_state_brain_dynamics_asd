load('./Connectivity_significant_matrix.mat');
C_read = ConSig;
C_emp_read = ConSigEmp;
C_emp_read(isinf(C_emp_read))=1;
C_emp_read_TD = zeros(105,105,2);
C_read_TD = zeros(105,105,2);
%%%%%%%Averaging out the data%%%

for i = 1:34
    C_emp_read_TD(:,:,1) = C_emp_read_TD(:,:,1) + C_emp_read(:,:,i);
    C_read_TD(:,:,1) = C_read_TD(:,:,1) + C_read(:,:,i);
end
C_emp_read_TD(:,:,1) = C_emp_read_TD(:,:,1)/34;
C_read_TD(:,:,1) = C_read_TD(:,:,1)/34;

for i = 34:68
    C_emp_read_TD(:,:,2) = C_emp_read_TD(:,:,2) + C_emp_read(:,:,i);
    C_read_TD(:,:,2) = C_read_TD(:,:,2) + C_read(:,:,i);
end
C_emp_read_TD(:,:,2) = C_emp_read_TD(:,:,2)/34;
C_read_TD(:,:,2) = C_read_TD(:,:,2)/34;

ConSig=C_read_TD;
ConSigEmp = C_emp_read_TD;

save('Connectivity_significant_matrix_AVG.mat','ConSig','ConSigEmp');

