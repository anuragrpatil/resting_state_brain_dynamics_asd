function [S_persec frN] = DMF_eulers_explicit_EIC(T,dt,C,G,noiseAmp)

tmax=T(end);
tspan1=0:dt:tmax;
[m nAreas] = size(C)

taon=100;
taog=10;
gamma=0.641;
sigma=noiseAmp;
% JN=0.15;
JN = 0.2609;
J=0.9*ones(nAreas,1);
% I0=0.382;
I0=0.3;
Jexte=1.0;
Jexti=0.7;
w=1.4;
%Iext=0.02;
curr=zeros(tmax,nAreas);
curr_new=zeros(tmax,nAreas);
neuro_act=zeros(tmax,nAreas);
S_persec=zeros(tmax,nAreas);
S_persec=zeros(tmax,nAreas);
ix=1;
we=G;   %we=G

a= 270 ;
b= 108;
d = 0.154; 
g=0.087;
I=177.;
c=615.; 

f_Hn = @(x) (a*x-b)./(1-exp(-1*d*(a*x-b)));
f_Hg = @(x) (c*x-I)./(1-exp(-1*g*(c*x-I)));

k=1;
sn=0.001*ones(nAreas,1);
sg=0.001*ones(nAreas,1);
Iext=zeros(nAreas,1);
delta=0.02*ones(nAreas,1);

for k=1:500
 sn=0.001*ones(nAreas,1);
 sg=0.001*ones(nAreas,1);
 nn=1;
 j=0;
for i=2:1:length(tspan1)
  xn=I0*Jexte+w*JN*sn+we*JN*C*sn-J.*sg;
  xg=I0*Jexti+JN*sn-sg;
  rn=f_Hn(xn);
  rg=f_Hg(xg);
  sn=sn+dt*(-sn/taon+(1-sn)*gamma.*rn./1000.)+sqrt(dt)*sigma*randn(nAreas,1);
  sn(sn>1) = 1;  
  sn(sn<0) = 0;      
  sg=sg+dt*(-sg/taog+rg./1000.)+sqrt(dt)*sigma*randn(nAreas,1);
  sg(sg>1) = 1;        
  sg(sg<0) = 0;
  j=j+1;
  if j==10
   curr(nn,:)=xn'-125/310;
   S_persec(nn,:)=rn';
   frN(nn,:) = rn';
   nn=nn+1;
   j=0;
   
  end
 end

 currm=mean(curr(1000:end,:),1);
 flag=0;
 for n=1:1:nAreas
  %if (n==2&3) && abs(currm(n)+0.026)>0.005 
  if abs(currm(n)+0.026)>0.005
   if currm(n)<-0.026 
    J(n)=J(n)-delta(n);
    delta(n)=delta(n)-0.001;
    if delta(n)<0.001;
       delta(n)=0.001;
    end
   else 
    J(n)=J(n)+delta(n);
   end
  else
   flag=flag+1;
  end
 end

 if flag==nAreas
  break;
 end
end
