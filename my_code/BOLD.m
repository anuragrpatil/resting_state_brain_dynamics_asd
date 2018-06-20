% based on Gustavo's code implementation of Hemodynamic model
 
function [b] = BOLD(T,S)

% The Hemodynamic model with one simplified neural activity
% 
% BOLD(f0,T,n_fig,suite)
%
% T       : total time (s)

% global itaus itauf itauo ialpha Eo dt

ch_int = 0;         % 0: Euler, 1: ode45

dt  = 0.001;        % (s)
t0  = (0:dt:T)';
n_t = (length(t0)-2);

t_min = 20;
n_min = round(t_min/dt);

r_max = max(S);

% BOLD model parameters

taus   = 0.65; % 0.8;    % time unit (s)
tauf   = 0.41; % 0.4;    % time unit (s)
tauo   = 0.98; % 1;      % mean transit time (s)
alpha  = 0.32; % 0.2;    % a stiffness exponent
itaus  = 1/taus;
itauf  = 1/tauf;
itauo  = 1/tauo;
ialpha = 1/alpha;
Eo     = 0.34; % 0.8;    % resting oxygen extraction fraction
vo     = 0.02;
k1     = 7*Eo; 
k2     = 2; 
k3     = 2*Eo-0.2;

% Initial conditions

x0  = [0 1 1 1];


% Euler method

    t      = t0;
    x      = zeros(n_t,4);
    x(1,:) = x0;
    for t = 1:n_t-1
        x(t+1,1) = x(t,1) + dt*(S(t,1)-taus*x(t,1)-tauf*(x(t,2)-1) );
        x(t+1,2) = x(t,2) + dt*x(t,1);
        x(t+1,3) = x(t,3) + dt*itauo*(x(t,2)-x(t,3)^ialpha);
        x(t+1,4) = x(t,4) + dt*itauo*(x(t,2)*(1-(1-Eo)^(1/x(t,2)))/Eo - (x(t,3)^ialpha)*x(t,4)/x(t,3));
    end



% t  = t(n_min:end);
% s  = x(n_min:end,1);
% fi = x(n_min:end,2);
% v  = x(n_min:end,3);
% q  = x(n_min:end,4);

t  = t(:);
s  = x(:,1);
fi = x(:,2);
v  = x(:,3);
q  = x(:,4);


b  = vo*( k1.*(1-q) + k2.*(1-q./v) + k3.*(1-v) );

clear x;

