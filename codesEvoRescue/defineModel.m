%%% Define the model for drug resistance
%%% Teemu Kuosmanen

%%% IMPORTANT: RUN THIS FIRST TO DEFINE THE PARAMETERS TO WORKSPACE

% PARAMETERS
K=1; % scaled carrying capacity
rs=0.5; % growth rate for S
rr=0.8*rs; % growth rate R
m0=10^-6; % control independent, baseline mutation rate
dmax=2*rs; % maximal rate of death as u -> umax
umax=1000; % constraint 1: maximal rate of control
alfa=1*umax^-1*10^-5; % linear coefficient for control dependent mutation rate 
h=40; % rate of control which yields 50 % of max. effect dmax
k=2.3; % steepness of sigmoidal function, hill-function coefficient 
S0=0.99; % initial condition for S, close to K to obtain control law u(S)
R0=0; % initial condition for R

% INITIALIZE
dt=0.001; % time step, make smaller if necessary
Tf=35; % end time
T=0:dt:Tf; % time vector
n=length(T);
S=zeros(1,n); % sensitive pop size as function of time
R=zeros(1,n); % resistent pop size as function of time
S(1)=S0;
R(1)=R0;
lambda1=zeros(1,n); % multiplier function for S
lambda2=zeros(1,n); % multiplier function for R
u=zeros(1,n); % optimal control strategy as function of time

% COST FUNCTION: C(u)= \int_{0}^{Tf} S(t,u(t))*mu(u(t))*dt

% Boundary conditions for the multipliers
lambda1(end)=0;
lambda2(end)=0;

% HAMILTONIAN: H(u) = S*mu + lambda1*fS +lambda2*fR

% Define the model DYNAMICS using function handles
d=@(t,u) -dmax/(1+(u/h)^k)+dmax; % Hill-type sigmoidal dose-response, d(u)
dd=@(u) ((dmax)*(k/h^k)*u^(k-1))/((1+(u/h)^k)^2); % derivative of dose response, d'(u)
mu=@(t,u) m0+alfa*u; % linear dose-dependent mutation rate, \mu(u)

% % independent logistic growth dynamics for S and R
% fS=@(t,S,R,u) rs*S*(1-S/K)-d(t,u)*S-mu(t,u)*S; % logistic growth dynamics for sensitive cells dS/dt; try also gompertz: rs*S*log(K/(S+R))
% fR=@(t,S,R,u) rr*R*(1-R/K)+mu(t,u)*S; % logistic growth dynamics for resistant cells dR/dt
% L1=@(t,lambda1,lambda2,S,R,u) -mu(t,u)-lambda1*(rs*(1-(2*S)/K)-d(t,u)-mu(t,u))-lambda2*mu(t,u); % sensitive multiplier dynamics, dlambda1/dt  
% L2=@(t,lambda1,lambda2,S,R,u) -lambda2*rr*(1-(2*R)/K); % resistant multiplier dynamics, dlambda2/dt

% common carrying capacity
fS=@(t,S,R,u) rs*S*(1-(S+R)/K)-d(t,u)*S-mu(t,u)*S; % logistic growth dynamics for sensitive cells dS/dt; try also gompertz: rs*S*log(K/(S+R))
fR=@(t,S,R,u) rr*R*(1-(S+R)/K)+mu(t,u)*S; % logistic growth dynamics for resistant cells dR/dt
L1=@(t,lambda1,lambda2,S,R,u) -mu(t,u)-lambda1*(rs*(1-(2*S+R)/K)-d(t,u)-mu(t,u))-lambda2*(mu(t,u)-rr*R/K); % sensitive multiplier dynamics, dlambda1/dt  
L2=@(t,lambda1,lambda2,S,R,u) lambda1*rs*S/K-lambda2*rr*(1-(2*R+S)/K); % resistant multiplier dynamics, dlambda2/dt

%% Plot pharmacodynamics and dose-dependent mutation rate
U=logspace(0,3,n);  % Create a log-spaced vector from 0 to umax
D=zeros(1,n); % d(U)
DD=zeros(1,n); % d'(U)
for i=1:n
    t=T(i);
    D(i)=d(t,U(i));
    DD(i)=dd(U(i));
end

% Figure 1B
figure(1)
subplot(1,2,1)
semilogx(U,rs-D,'LineWidth',5)
title('Pharmacodynamics', 'fontsize', 20)
xlabel('Drug concentration u', 'fontsize', 18)
ylabel('Cell Growth rate', 'fontsize', 18)
set(gca,'FontSize', 18, 'FontWeight', 'bold')
xline(h,'--')

subplot(1,2,2)
semilogx(U,mu(0,U)/m0,'LineWidth',5)
title('Linear dose-dependent mutation rate', 'fontsize', 20)
xlabel('Drug concentration u', 'fontsize', 18)
ylabel('Fold change to mutation rate', 'fontsize', 18)
set(gca,'FontSize', 18, 'FontWeight', 'bold')
xline(h,'--')























