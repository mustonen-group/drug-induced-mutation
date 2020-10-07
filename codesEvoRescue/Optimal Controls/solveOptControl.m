%%% Solve the optimal control strategy using Pontryagin's minimum principle
%%% Teemu Kuosmanen

%%% IMPORTANT: Run first defineModel.m 

%%% If the algorithm does not converge, try adjusting the initial proposals
%%% and plot the switching function W at different time points to see if
%%% the root searching interval for fzero should be adjusted. 
%%% 
%%% If you change the model dynamics (Hamiltonian), make sure that 
%%% switching function (fun and W at lines 105 and 108) corresponds 
%%% to dH/du.

% Initialize simulation variables
test=-1; % convergence variable
iter=1; % record number of iterations
delta=0.001; % convergence threshold
max_iter=1000; % maximum number of iterations
us=zeros(max_iter,n); % save optimal controls for diagnostic purposes
tests=zeros(max_iter,1);

% Initial educated guesses
u=500*ones(1,n); % finetune the proposal if necessary
S=zeros(1,n); % sensitive pop size as function of time
R=zeros(1,n); % resistent pop size as function of time
% lambdas must not be equal or zero!
lambda1=0.01*ones(1,n); % multiplier function for S
lambda2=0.1*ones(1,n); % multiplier function for R

while(test < 0)
    if iter > max_iter
        disp('Out of iterations')
        break
    else 
        iter
    end
    
    oldu=u;
    oldS=S;
    oldR=R;
    oldlambda1=lambda1;
    oldlambda2=lambda2;
    
    % Solve S and R FORWARD in time with Runge-Kutta and initial guesses
    S(1)=S0;
    R(1)=R0;
    for i=2:n
        tk=T(i-1);
        sk=S(i-1);
        rk=R(i-1);
        uk=u(i-1);
        % S
        k1=dt*fS(tk,sk,rk,uk);
        k2=dt*fS(tk+k1/2,sk+k1/2,rk+k1/2,uk+k1/2);
        k3=dt*fS(tk+k2/2,sk+k2/2,rk+k2/2,uk+k2/2);
        k4=dt*fS(tk+k3,sk+k3,rk+k3,uk+k3);
        S(i)=sk+1/6*(k1+2*(k2+k3)+k4);
        % R
        k1=dt*fR(tk,sk,rk,uk);
        k2=dt*fR(tk+k1/2,sk+k1/2,rk+k1/2,uk+k1/2);
        k3=dt*fR(tk+k2/2,sk+k2/2,rk+k2/2,uk+k2/2);
        k4=dt*fR(tk+k3,sk+k3,rk+k3,uk+k3);
        R(i)=rk+1/6*(k1+2*(k2+k3)+k4);
    end
    S=real(S);
    R=real(R);
    
    % Solve the multipliers BACKWARDS in time with the obtained S and R
    lambda1(end)=0;
    lambda2(end)=0;
    for i=n-1:-1:1
        tk=T(i+1);
        Lk=lambda1(i+1);
        Lk2=lambda2(i+1);
        sk=S(i+1);
        rk=R(i+1);
        uk=u(i+1);
        % solve lambda1
        k1=dt*L1(tk,Lk,Lk2,sk,rk,uk);
        k2=dt*L1(tk,Lk+k1/2,Lk2+k1/2,sk+k1/2,rk+k1/2,uk+k1/2);
        k3=dt*L1(tk,Lk+k2/2,Lk2+k2/2,sk+k2/2,rk+k2/2,uk+k2/2);
        k4=dt*L1(tk,Lk+k3,Lk2+k3,sk+k3,rk+k3,uk+k3);
        lambda1(i)=Lk-1/6*(k1+2*(k2+k3)+k4);
        % solve lambda2
        k1=dt*L2(tk,Lk,Lk2,sk,rk,uk);
        k2=dt*L2(tk,Lk+k1/2,Lk2+k1/2,sk+k1/2,rk+k1/2,uk+k1/2);
        k3=dt*L2(tk,Lk+k2/2,Lk2+k2/2,sk+k2/2,rk+k2/2,uk+k2/2);
        k4=dt*L2(tk,Lk+k3,Lk2+k3,sk+k3,rk+k3,uk+k3);
        lambda2(i)=Lk2-1/6*(k1+2*k2+2*k3+k4);    
    end
    lambda1=real(lambda1);
    lambda2=real(lambda2);
    
    % UPDATE u based on the updated trajectories and multipliers
    U=logspace(1,3,n);
    %U=linspace(0,umax,n);
    DD=zeros(1,n);
    for i=1:n
        DD(i)=dd(U(i)); % Derivative of the dose response d'(u)
    end
    
    % Search the roots of the switching function H_u
    Z1=zeros(1,n); % the first root is never optimal (maximizing), prove using Legendre-Clebsch
    Z2=zeros(1,n); % the second root is the optimal (minimizing)
    %W=zeros(n,n);
    for i=1:n
        % Numerically search the 2nd root of the switching function
        fun=@(u) S(i)*(alfa-lambda1(i)*(dd(u)+alfa))+lambda2(i)*S(i)*alfa; % add lambda3 for isoperimetric constraint
        % Plot W to visualize the roots and adjust the root finding
        % algorithm if necessary
        W=S(i)*(alfa-lambda1(i)*(DD+alfa))+lambda2(i)*S(i)*alfa;
%         [x_end,ind]=min(abs(W));
%         Z2(i)=U(ind);
        [x_end,ind]=min(W);

        if x_end < 0 && fun(U(end)) > 0
             Z2(i)=fzero(fun,[U(ind),U(end)]);
        elseif x_end > 0
            Z2(i)=0;
        else
            Z2(i)=U(end);
        end

%         % Plot for illustrating root finding
%         if sum(i==[1, 5001, 10001, 20001])>0
%             plot(U,W,'LineWidth',2)
%             xlim([0 125])
%             %leg(plot_iter)=['t=',num2str(T(i))];
%             hold on
%             %plot_iter=plot_iter+1;
%         end
        
%       % Alternative way to try...
%         if x_end < 0 && fun(U(1)) > 0
%             Z1(i)=fzero(fun,[U(1),U(ind)]);
%         elseif x_end < 0 && fun(U(end)) > 0
%             Z2(i)=fzero(fun,[U(ind),U(end)]);
%         elseif x_end > 0
%             Z1(i)=0;
%             Z2(i)=0;
%         else 
%             Z1(i)=U(end);
%             Z2(i)=U(end);
%         end
    end
    us(iter,:)=Z2; % collect updated u for each iteration
    u=Z2;
    
    % check for convergence
    temp1=delta*sum(abs(u))-sum(abs(oldu-u));
    temp2=delta*sum(abs(S))-sum(abs(oldS-S));
    temp3=delta*sum(abs(lambda1))-sum(abs(oldlambda1-lambda1));
    temp=[temp1,temp2,temp3]; 
    test=min(temp);
    tests(iter)=test; % record to see convergence
    
    iter=iter+1;
end
    

% Solve R for the optimal control case
% Solve for R forward
R(1)=R0;
for i=2:n
    tk=T(i-1);
    sk=S(i-1);
    rk=R(i-1);
    uk=u(i-1);
    k1=dt*fR(tk,sk,rk,uk);
    k2=dt*fR(tk+k1/2,sk+k1/2,rk+k1/2,uk+k1/2);
    k3=dt*fR(tk+k2/2,sk+k2/2,rk+k2/2,uk+k2/2);
    k4=dt*fR(tk+k3,sk+k3,rk+k3,uk+k3);
    R(i)=rk+1/6*(k1+2*(k2+k3)+k4);
end

%% Plot the results
figure(2)
subplot(1,3,1)
plot(T,u,'LineWidth',4,'color',[0.4940, 0.1840, 0.5560])
title('Optimal control strategy','fontsize',18)
ylabel('Control','fontsize', 16)
xlabel('Time', 'fontsize',16)
set(gca,'FontSize', 18, 'FontWeight', 'bold')

subplot(1,3,2)
plot(T,S,'LineWidth',4,'color',[0, 0.4470, 0.7410])
hold on
plot(T,R,'LineWidth',4,'color',[0.8500, 0.3250, 0.0980])
ylabel('Population size','fontsize', 16)
xlabel('Time','fontsize', 16)
legend('sensitive', 'resistant')
title('Trajectories','fontsize',18)
set(gca,'FontSize', 18, 'FontWeight', 'bold')

subplot(1,3,3)
plot(T,lambda1,'LineWidth',4,'color',[0, 0.4470, 0.7410])
hold on
plot(T,lambda2,'LineWidth',4,'color',[0.8500, 0.3250, 0.0980])
ylabel('Multiplier value','fontsize', 16)
xlabel('Time','fontsize', 16)
legend('sensitive', 'resistant')
title('Multiplier dynamics','fontsize',18)
set(gca,'FontSize', 18, 'FontWeight', 'bold')


cost_opt=dt*trapz(S.*(m0+alfa*u)); % evaluate the cost of the optimal trajectory
u_cum=dt*trapz(u); % the cumulative amount of control used


%% Compare results 
Us=[mean(u),40,60,120,200,300]; % constant doses to compare
%Us=[mean(u),u(1),max(u)];
Ss=zeros(length(Us),n);
Rs=zeros(length(Us),n);
csts=zeros(1,length(Us));
Ss(:,1)=S0;
Rs(:,1)=R0;
leg=strings(1,length(Us)+1);
leg(1)='u_{opt}';

for j=1:length(Us)
    uu=Us(j)*ones(1,n);
    for i=2:n
        tk=T(i-1);
        sk=Ss(j,i-1);
        rk=Rs(j,i-1);
        uk=uu(i-1);
        % solve S 
        k1=dt*fS(tk,sk,rk,uk);
        k2=dt*fS(tk,sk+k1/2,rk+k1/2,uk+k1/2);
        k3=dt*fS(tk,sk+k2/2,rk+k2/2,uk+k2/2);
        k4=dt*fS(tk,sk+k3,rk+k3,uk+k3);
        Ss(j,i)=sk+1/6*(k1+2*k2+2*k3+k4);
        % solve R
        k1=dt*fR(tk,sk,rk,uk);
        k2=dt*fR(tk,sk+k1/2,rk+k1/2,uk+k1/2);
        k3=dt*fR(tk,sk+k2/2,rk+k2/2,uk+k2/2);
        k4=dt*fR(tk,sk+k3,rk+k3,uk+k3);
        Rs(j,i)=rk+1/6*(k1+2*k2+2*k3+k4);
    end
    leg(j+1)=['u=',num2str(Us(j))];
end
leg(2)='mean(u_{opt})';

% subplot(1,3,1)
% plot(T,S,'Linewidth',2,'color','blue');
% hold on
% plot(T,Ss, 'Linewidth',2);
% title('Compare trajectories')
% legend(leg)
% xlabel('Time')
% ylabel('S(t)')
figure(4)
subplot(1,2,2)
plot(T,R,'linewidth',3)
hold on
plot(T,Rs,'linewidth',3)
title('Comparing resistant growth','fontsize',14)
xlim([20, Tf])
xlabel('Time', 'fontsize',12,'fontweight','bold')
ylabel('R(t)','fontsize',12,'fontweight','bold')
legend(leg)


subplot(1,2,1)
plot(T,dt*cumtrapz(S.*(m0+alfa*u)),'linewidth',3)
hold on
for j=1:length(Us)
    csts(j)=dt*trapz(Ss(j,:).*(m0+alfa*Us(j)*ones(1,n)));
    plot(T,dt*cumtrapz(Ss(j,:).*(m0+alfa*Us(j)*ones(1,n))), 'linewidth',3)
end
legend(leg)
title('Comparing cumulative costs','fontsize',14)
xlabel('Time','fontsize',12,'fontweight','bold')
ylabel('Cumulative cost C(t)','fontsize',12,'fontweight','bold')