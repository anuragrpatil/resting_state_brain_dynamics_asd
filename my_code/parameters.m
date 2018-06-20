function parameter = parameters(C)
global w 
global a 
global b
global d
global Jn
global gamma
global I0
global sigma
global tau_s
global n 

w = 0.9 ;
a= 270 ;
b= 108;
d = 0.154; 
Jn = 0.2609;
gamma = 0.641*10^-6;
I0 = 0.3;
sigma = 0.001;
tau_s = 100*10^-3;
[m,n] = size(C);

