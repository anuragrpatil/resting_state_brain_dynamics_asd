function [Sn Hn xn] = DMF_eulers_explicit(T,dt,C,G,noiseAmp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Parameters for DMF%%%%%%%%%%%%%%%%%%
w = 0.9 ;          % wEE or w+ synpatics weight for intra area excitatory-excitatory
% w = 1.4
a= 270 ;
b= 108;
d = 0.154; 
JN = 0.2609; % synaptic coupling for NMDA synapses
% JN = 0.15
gamma = 0.641;
I0 = 0.3;          % external input current
% I0 = 0.382
sigma = noiseAmp;
taon = 100;        % time constant for NMDA synapses
[m nAreas] = size(C);

% taog=10;
% 
% 
% 
% g=0.087;
% I=177.;
% c=615. 
% 
% JN=0.15;
% J=0.9*ones(nAreas,1);
% 
% Jexte=1.0;
% Jexti=0.7;




%%%%%%%%%%%%%%%%%%%%%%%%%Initialization%%%%%%%%%%%%%%%%%%%%
Sn = 0.001*ones(length(T)-1,nAreas);
Sg = 0.001*ones(length(T)-1,nAreas);
xn = zeros(length(T)-1,nAreas);
xg = zeros(length(T)-1,nAreas);
Hn = zeros(length(T)-1,nAreas);
Hg = zeros(length(T)-1,nAreas);



%%%%%%%%%%%%%%%%%%%%%%%%% Model %%%%%%%%%%%%%%%%%%%%%%
 
f_Sn = @(S,H) -1*S/taon + (gamma/1000)*(1-S).*H;
f_Hn = @(x) (a*x-b)./(1-exp(-1*d*(a*x-b)));
f_Hg = @(x) (c*x-I)./(1-exp(-1*g*(c*x-I)));

for t = 1:1:(length(T)-1)
%     for i = 1:1:nAreas
%         
%         x(t,i)= w*JN*S(t,i)+G*JN*sum(C(i,:).*S(t,1:nAreas))+I0;
%         
%         H(t,i)= f_H(x(t,i));
%         S(t+1,i) = S(t,i)+dt*f_S(S(t,i),H(t,i))+sqrt(dt)*sigma*randn;
%         
%         if S(t+1,i)>1
%          S(t+1,i)= 1; 
%         elseif S(t+1,i)<0
%             S(t+1,i)=0;
%         end  
%         
%     end

        xn(t,:)= w*JN*Sn(t,:)+G*JN*(C*Sn(t,1:nAreas)')'+I0;


        Hn(t,:)= f_Hn(xn(t,:));

        Sn(t+1,:) = Sn(t,:)+dt*f_Sn(Sn(t,:),Hn(t,:))+dt^(1/2)*sigma*randn(1,nAreas);
        
        if Sn(t+1,:)>1
         Sn(Sn(t+1,:)>1)= 1; 
        elseif Sn(t+1,:)<0
            Sn(Sn(t+1,:)<0)=0;
        end 



        
%         xn(t,:)= I0*Jexte+w*JN*Sn(t,:)+G*JN*(C*Sn(t,1:nAreas)')'+I0;
%         xg(t,:)= I0*Jexti+JN*Sg(t,:)+JN*Sn(t,:)-Sg(t,:);
% 
%         Hn(t,:)= f_Hn(xn(t,:));
%         Hg(t,:)= f_Hg(xg(t,:));
% 
%         
%         Sn(t+1,:) = Sn(t,:)+dt*f_Sn(Sn(t,:),Hn(t,:))+dt^(1/2)*sigma*randn(1,nAreas);
%         Sn(Sn>1)=1;
%         Sn(Sn<0)=0;
% 
%         Sg(t+1,:)=Sg(t,:)+dt*(-Sg(t,:)/taog+Hg(t,:)./1000.)+sqrt(dt)*sigma*randn(1,nAreas);
%         Sg(Sg>1)=1;
%         Sg(Sg<0)=0;
       
         
        
        
    
  
   
 
end

    
   
