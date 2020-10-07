%%% Solve control law u(S) from the stationary profile (EQUATION 5)
%%% Teemu Kuosmanen

%%% IMPORTANT: Run first defineModel.m, solveOptControl.m

% The inverse function method
TS=linspace(0,1,10000); % discretize S
US=linspace(0,umax,10000); % discretize U
nn=length(TS);
inds=zeros(1,nn);
U=zeros(1,nn);
for i=1:nn
    [~,ind]=min(abs(S-TS(i))); % find the time when the pop size was TS for each element 
    inds(i)=ind;
    U(i)=u(ind); % record what was the control used at this time
end
figure(6)
plot(TS,U,'LineWidth',3)
xlabel('Population size S')
ylabel('Optimal control u(S)')
hold on

% The stationary profile equation (does not apply to discounted problem)
% EQ: S*mu(u)+fS*alfa/(d'(u)+alfa)=0
UU=zeros(1,nn);
D=zeros(1,nn);
DD=zeros(1,nn);
for i=1:nn
    D(i)=d(0,US(i)); 
    DD(i)=dd(US(i));
end

for i=1:nn
    FS=mu(0,US)+alfa./(DD+alfa).*(rs*(1-TS(i))-D-mu(0,US));
    [m,ind]=min(abs(FS));
    UU(i)=US(ind);
end
plot(TS,UU,'Linewidth',3)
title('Control law','fontsize',18)
hold on
yline(UU(1),'--', 'linewidth',3)
legend('Control law from numerical solution', 'Control law from Eq. 5','Optimal constant dose u(0)')
set(gca,'FontSize', 18, 'FontWeight', 'bold')

