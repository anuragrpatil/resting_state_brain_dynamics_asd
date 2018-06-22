function [S_persec frN] = DMF_eulers_explicit(T,dt,C,G,noiseAmp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Parameters for DMF%%%%%%%%%%%%%%%%%%
w = 0.65 ;          % wEE or w+ synpatics weight for intra area excitatory-excitatory
%  w = 1.4
a= 270 ;
b= 108;
d = 0.154; 
 JN = 0.2609; % synaptic coupling for NMDA synapses
%  JN = 0.15
gamma = 0.641;
I0 = 0.3;
% I0 = 0.24;          % external input current
% I0 = 0.382
sigma = noiseAmp;
taon = 100;        % time constant for NMDA synapses
[m nAreas] = size(C);

% taog=10;
% g=0.087;
% I=177.;
% c=615.; 
% 
% JN=0.15;
% J=0.9*ones(nAreas,1);
% 
% Jexte=1.0;
% Jexti=0.7;




%%%%%%%%%%%%%%%%%%%%%%%%%Initialization%%%%%%%%%%%%%%%%%%%%
Sn = 0.001*ones(1,nAreas);
Sg = 0.001*ones(1,nAreas);
xn = zeros(1,nAreas);
xg = zeros(1,nAreas);
Hn = zeros(1,nAreas);
Hg = zeros(1,nAreas);
frN = zeros(T(end),nAreas);
S_persec = zeros(T(end),nAreas);
j = 0;
nn = 0;


%%%%%%%%%%%%%%%%%%%%%%%%% Model %%%%%%%%%%%%%%%%%%%%%%
 
f_Sn = @(S,H) -1*S/taon + (gamma/1000)*(1-S).*H;
f_Hn = @(x) (a*x-b)./(1-exp(-1*d*(a*x-b)));
f_Hg = @(x) (c*x-I)./(1-exp(-1*g*(c*x-I)));

for t = 1:1:(length(T)-1)
    
%     for i = 1:1:nAreas
%         
%         xn(t,i)= w*JN*Sn(t,i)+G*JN*sum(C(i,:).*Sn(t,1:nAreas))+I0;
%         
%         Hn(t,i)= f_Hn(xn(t,i));
%         Sn(t+1,i) = Sn(t,i)+dt*f_Sn(Sn(t,i),Hn(t,i))+sqrt(dt)*sigma*randn;
%         
%         if Sn(t+1,i)>1
%          Sn(t+1,i)= 1; 
%         elseif Sn(t+1,i)<0
%             Sn(t+1,i)=0;
%         end  
%         
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%vectorized code without inhibition%%%%%%%%%%

      xn = w*JN*Sn + G*JN*(C*Sn')' + I0; 
      Hn = f_Hn(xn); 
      Sn = Sn + dt*f_Sn(Sn,Hn) + sqrt(dt)*sigma*randn(1,nAreas);
      Sn(Sn>1) = 1;
      Sn(Sn<0) = 0;
      j = j + 1;
      if(j == 1/dt)
        nn = nn + 1;
        j = 0;
        S_persec(nn,:) = Sn;
        frN(nn,:) = Hn;
      end




%%%%%%%%%%%%%%%%%%%inhibition control added%%%%%%%%%%%%%%%%%%%        
%         xn(t,:)= I0*Jexte+w*JN*Sn(t,:)+G*JN*(C*Sn(t,1:nAreas)')'+I0;
%         xg(t,:)= I0*Jexti+JN*Sg(t,:)+JN*Sn(t,:)-Sg(t,:);
% 
%         Hn(t,:)= f_Hn(xn(t,:));
%         Hg(t,:)= f_Hg(xg(t,:));
% 
%         
%         Sn(t+1,:) = Sn(t,:)+dt*f_Sn(Sn(t,:),Hn(t,:))+dt^(1/2)*sigma*randn(1,nAreas);
%         if Sn(t+1,:)>1
%          Sn(Sn(t+1,:)>1)= 1; 
%         elseif Sn(t+1,:)<0
%             Sn(Sn(t+1,:)<0)=0;
%         end 
% 
%         Sg(t+1,:)=Sg(t,:)+dt*(-Sg(t,:)/taog+Hg(t,:)./1000.)+sqrt(dt)*sigma*randn(1,nAreas);
%         if Sg(t+1,:)>1
%          Sg(Sg(t+1,:)>1)= 1; 
%         elseif Sg(t+1,:)<0
%             Sg(Sg(t+1,:)<0)=0;
%         end 
%        
         
        
        
    
  
   
 
end

    
   
