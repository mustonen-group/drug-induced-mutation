%%% Plot control laws
%%% Teemu Kuosmanen

%%% Run the files defineModel.m, solveOptControl.m, solveDiscounted.m first
%%% to get the optimal time-dependent controls u(t) and corresponding trajectories S(t)

%% Use your own data or use the examples below
u1=load('u(t).mat','u');
u1=u1.u; % rename as u1
u2=load('u(t)_discounted.mat','u');
u2=u2.u; % rename the discount solution as u2
S1=load('S(t).mat','S');
S1=S1.S; % rename as u1
S2=load('S(t)_discounted.mat','S');
S2=S2.S; % rename the discount solution as u2

%% The inverse function method
TS=linspace(0,1,10000); % discretize S
US=linspace(0,umax,10000); % discretize U
nn=length(TS);
inds=zeros(1,nn);
U1=zeros(1,nn);
U2=zeros(1,nn);
for i=1:nn
    % no discount
    [~,ind]=min(abs(S1-TS(i))); % find the time when the pop size was TS for each element 
    inds(i)=ind;
    U1(i)=u1(ind); % record what was the control used at this time
    % discount
    [~,ind]=min(abs(S2-TS(i))); % find the time when the pop size was TS for each element 
    inds(i)=ind;
    U2(i)=u2(ind);
end
figure(6)
plot(TS,U1,'LineWidth',3,'color',[0.4940, 0.1840, 0.5560])
hold on
plot(TS,U2,'LineWidth',3,'color',[0.9290, 0.6940, 0.1250])
yline(U1(2),'--','LineWidth',3,'color',[0.4940, 0.1840, 0.5560])
yline(U2(end),'--','LineWidth',3,'color',[0.9290, 0.6940, 0.1250])
xlabel('Population size S')
ylabel('Optimal control u(S)')
legend('no discount', 'discounted')
title('Control laws')
set(gca,'FontSize', 18, 'FontWeight', 'bold')

