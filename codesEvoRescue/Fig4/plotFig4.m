%%% Plot figure 4
%%% Teemu Kuosmanen

%%% Simulate data with stochasticSimulation.m 
%%% to get vectors finalX.mat or plot the ready data below

%% Uncomment if you want to plot the ready data used in article
% load('finalN_NEW.mat');
% load('finalR_NEW.mat');
% load('finalS_NEW.mat');
% Nsim=length(finalN);
% doses=[0,20,40,50,55:5:120,130:10:200,250,300,350]; % run through different constant concentrations
% K=1/m0;

%% 
subplot(1,3,1)
H=histogram(finalN(:,14)/K,'normalization','probability'); % modidy binedges to correct for 0-bin if necessary
title('Bimodal distribution (for u=104.5)')
xlabel('Final population size N(T)/K')
ylabel('Normalized frequency')
set(gca,'FontSize', 18, 'FontWeight', 'bold')
xmean=mean(finalN(finalN(:,14)>0,14)/K);
pcure=sum(finalN(:,14)==0)/Nsim;
line([xmean xmean],[0 0.1],'LineWidth',3,'color','r')
line([0,0.2],[pcure pcure],'LineWidth',3,'color','g')
text(0.21,pcure,'Probability of cure','fontweight','bold','fontsize',14)
text(xmean-0.25,0.11,'Expected rescue fraction','fontweight','bold','fontsize',14)


subplot(1,3,2)
Pcure=sum(finalN==0,1)/Nsim;
err=sqrt(Pcure.*(1-Pcure));
%errorbar(doses,Pcure,err,'-o','linewidth',2)
plot(doses,Pcure,'-o','linewidth',4,'Markersize',10)
ylim([0,1])
xlabel('Drug concentration u','fontsize',15,'fontweight','bold')
ylabel('Probability of cure, N(T)=0','fontsize',15,'fontweight','bold')
xline(104.5,'--','linewidth',2)
title('Probability of cure','fontsize',15,'fontweight','bold')
set(gca,'FontSize', 18, 'FontWeight', 'bold')

% calculate the expected rescue fraction conditioned on non-extinction
RF=zeros(1,length(doses));
for i=1:length(doses)
  RF(i)=mean(finalN(finalN(:,i)>0,i)/K); % calculate mean of non-extinct populations
end
subplot(1,3,3)
plot(doses,RF,'-o','linewidth',4,'Markersize',10)
xlabel('Drug concentration u','fontsize',15,'fontweight','bold')
ylabel('Expected final population size N(T)/K','fontsize',15,'fontweight','bold')
xline(60,'--','linewidth',2)
title('Expected rescue fraction','fontsize',15,'fontweight','bold')
set(gca,'FontSize', 18, 'FontWeight', 'bold')



