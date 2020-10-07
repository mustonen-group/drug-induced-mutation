%%% Fig3: Contour plot with MTD vs. alpha
%%% Teemu Kuosmanen

%%% IMPORTANT: Run first defineModel.m 

alphas=logspace(-12,-7,30); % plotting region for parameter alpha
doses=[70:5:95,100:50:1000]; % plotting region for parameter MTD
nn=10000;
US=linspace(0,umax,nn); % discretize U
D=zeros(1,nn); % dose response d(U)
DD=zeros(1,nn); % derivative of dose response d'(U)
for i=1:nn
    D(i)=d(0,US(i)); 
    DD(i)=dd(US(i));
end

% solve the optimal constant dose for each alpha using Eq.5
uu=zeros(1,length(alphas)); % optimal constant doses
for i=1:length(alphas)
    FS=(m0+alphas(i)*US)+alphas(i)./(DD+alphas(i)).*(rs-D-(m0+alphas(i)*US));
    [~,ind]=min(abs(FS));
    uu(i)=US(ind);
end
plot(alphas,uu)

% calculate the cost of the optimal doses
opt_cost=zeros(1,length(uu));
S=ones(1,n);
for k=1:length(uu)
    mu=@(t,u) m0+alphas(k)*u; 
    fS=@(t,S,u) rs*S*(1-S/K)-d(t,u)*S-mu(t,u)*S; 
    for i=2:n
        tk=T(i-1);
        sk=S(i-1);
        uk=uu(k);
        % S
        k1=dt*fS(tk,sk,uk);
        k2=dt*fS(tk+k1/2,sk+k1/2,uk+k1/2);
        k3=dt*fS(tk+k2/2,sk+k2/2,uk+k2/2);
        k4=dt*fS(tk+k3,sk+k3,uk+k3);
        S(i)=sk+1/6*(k1+2*(k2+k3)+k4);
    end
    opt_cost(k)=dt*trapz(S.*(m0+alphas(k)*uu(k)));
end

% calculate the cost for each MTD and alpha
costs=zeros(length(alphas),length(doses));
S=ones(1,n);
for k=1:length(alphas)
    for j=1:length(doses)
        mu=@(t,u) m0+alphas(k)*u; 
        fS=@(t,S,u) rs*S*(1-S/K)-d(t,u)*S-mu(t,u)*S; 
        for i=2:n
            tk=T(i-1);
            sk=S(i-1);
            uk=doses(j);
            % S
            k1=dt*fS(tk,sk,uk);
            k2=dt*fS(tk+k1/2,sk+k1/2,uk+k1/2);
            k3=dt*fS(tk+k2/2,sk+k2/2,uk+k2/2);
            k4=dt*fS(tk+k3,sk+k3,uk+k3);
            S(i)=sk+1/6*(k1+2*(k2+k3)+k4);
        end
        costs(k,j)=dt*trapz(S.*(m0+alphas(k)*doses(j)));
    end
    costs(k,:)=costs(k,:)./opt_cost(k); % normalize with optimal cost
end

% surface plot
[X,Y]=meshgrid(alphas,doses);
Z=costs';
%surf(X,Y,Z)

%% Plot the contours
% suitable red colormap for background plotting
cmred=[linspace(1,0.8,50)',linspace(1,0.2,50)',linspace(1,0.2,50)'];

contourf(X,Y,1+X/m0.*Y,1:10,'linecolor','none')
colormap(cmred)
cb=colorbar;
set(get(cb,'title'),'string','FC to \mu_0');
cb.Ticks=1:10;
cb.TickLabels={'1','2','3','4','5','6','7','8','9','\geq10'};
hold on
xline(alfa,'--','linewidth',3)
[C,h]=contour(X,Y,Z,[1.25,1.5,2,3,4,5],'ShowText','on','color',[0.1 0.2 0.71],'linewidth',3,'labelspacing',880);
clabel(C,h,'FontSize',30,'color',[0.1 0.2 0.71])
plot(alphas,uu,'-.*','color',[0.1 0.2 0.71],'linewidth',3)
text(alphas(10), uu(10), '1','fontsize',30,'color',[0.1 0.2 0.71],'BackgroundColor', 'w','rotation',-45)
text(alphas(17), uu(17), '1','fontsize',30,'color',[0.1 0.2 0.71],'BackgroundColor', 'w','rotation',-30)
text(alphas(11), 105, '1.25','fontsize',25,'color',[0.1 0.2 0.71],'BackgroundColor', 'w','rotation',0)
xlabel('Parameter \alpha')
ylabel('MTD')
set(gca,'xscale','log')
set(gca,'FontSize', 22, 'FontWeight', 'bold')
