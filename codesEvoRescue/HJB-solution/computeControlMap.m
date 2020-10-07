%%% Hamilton-Jacobi-Bellman approach 
%%% Tommi Mononen, Teemu Kuosmanen

%%% IMPORTANT: Run first defineModel.m
%%% Dependence of SkellamVec.m help function

K = 100;        % carrying capacity
range=2;        % K*range must be larger than the largest population size

%alfa=10^-8; rs=0.5; dmax=1; h=40; k=2.3; m0=10^-6; Tf=35;
umax=150; 
zeroI=31; % The number of negative + zero population sizes  

% Transition tensor for all controls 0..umax
Wnnbc=zeros(umax+1,K*range+1,K*range+1); 
    
%% compute transition probabilities
fprintf('Computing transition matrices (u=0..%d) :\n',umax);
for u=0:umax
    death=dmax*(1-1/(1+(u/h)^k));
    %death=1/100*u; % try 'linearized' control leverage
    mutation=m0+alfa*u; 
        
    Wnn=zeros(K*range+1,K*range+zeroI);
 
    % use parallel processing for computing a transition matrix
    parfor S=0:(K*range)    
        Wnn(1+S,:)=skellamVec(-30-S,(K*range)-S,rs*S,death*S+mutation*S+rs*S^2/K);
    end
        
    % merge the probability mass of negative jumps
    tmp=[sum(Wnn(:,1:zeroI),2) Wnn(:,(zeroI+1):end)];
        
    % normalize probabilities
    Wnnbc(u+1,:,:)=tmp./repmat(sum(tmp,2),1,(K*range)+1);
        
    fprintf('%d ',u);
    if mod(u+1,20)==0
        fprintf('\n');
	end;
end;
fprintf('\n\nComputing the control map ...\n');

%% compute the control map with backpropagation
Jhjb=zeros(Tf+1,range*K+1);  % storage space for costs (no endpoint cost)
Chjb=zeros(Tf+1,range*K+1);  % optimal control
    
for t=(Tf+1):-1:2
    J=Jhjb(t,:);

    v=zeros(umax+1,range*K+1);
    costRec=zeros(umax+1,Tf+1,range*K+1);
        
    % for different controls
    for u=0:umax           
        % for different population sizes
        for sz=1:(range*K+1)
            v(u+1,sz)=J(:)'*reshape(Wnnbc(u+1,sz,:),1,range*K+1)';
        end;
        % add control cost of this round
        costRec(u+1,t-1,:)=v(u+1,:)+(0:(range*K))*(m0+alfa*u);%*exp(0.8*rS*(T-t-1)); % discounted version *exp(0.8*rS*(T(end)-t-1));
    end
        
    % find the control that minimizes cost
    [a,b]=min(costRec(:,t-1,:));
    Jhjb(t-1,:)=a; % save costs
    Chjb(t-1,:)=b; % save optimal u:s (index of optimal costs)
end;
Chjb=Chjb(1:Tf,:)-1; % shift and remove time T+1

%% --- plotting the control map ---
figure(1);
imagesc(Chjb'); 
set(gca,'Ydir','normal','linewidth',2,'fontsize',20);
hold on;
plot([1 50],[K K],'--k','linewidth',3); %draw line for carrying capacity
hold off;
xlabel('Time','fontsize',30);
ylabel('Population size','fontsize',30); 
title('Control map (value of u)','fontsize',30);
colorbar;


%% Plot the cost-to-go
% Notice that the surface is constant along time except for the very end
figure(2)
surf(Jhjb)
title('Cost-to-go surface J(t,S)')
ylabel('Time')
xlabel('Population size S')
zlabel('J(t,S)')
set(gca,'FontSize', 18, 'FontWeight', 'bold')

% Calculate and plot the time-derivative of J(t,S)
[~,FT]=gradient(Jhjb);

%%% The derivation of the control-law, Eq. 5, hinges on the fact that this surface is
%%% constant (here zero). 
figure(3)
surf(FT)
title('Time-dependence of cost-to-go')
ylabel('Time')
xlabel('Population size S')
zlabel('\partial_t J(t,S)')
set(gca,'FontSize', 18, 'FontWeight', 'bold')






