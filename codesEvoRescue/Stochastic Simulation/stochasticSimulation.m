%%% Stochastic simulation
%%% Teemu Kuosmanen

%%% IMPORTANT: Run first defineModel.m 

%%% Gillespie-algorithm to simulate the stochastic system [eqs.(8)-(9)].

%%% Run through the desired concentrations and the code records the 
%%% final population sizes.

%%% The number of simulations per dose is set to 100 (Fig4 data 2000)
%%% to allow experiments in reasonable time (~35 min). 

doses=[0,20,40,50,55:5:120,130:10:200,250,300,350]; % run through different constant concentrations
finalR=[]; % collect final No. resistant cells per sim per dose
finalS=[]; % collect final No. sensitive cells per sim per dose
finalN=[]; % collect final pop size (=the rescue fraction) per sim per dose
mutants=cell(1,length(doses)); % collect times of rescue mutants for each dose
exmutants=cell(1,length(doses)); % collect times of rescue mutants for each dose

% Specify birth and death rates for the sensitive and resistant cells
betaS=0.8;
deltaS=0.3;
betaR=0.5;
deltaR=0.1;

% How is K realized? d ≤ theta ≤ b (Anttila et al.)
thetaS=betaS; % constant birth rate, K realizes by increased death
thetaR=betaR;

% Simulation parameters
Nsim=100; % No. of virtual patients simulated with the dose
K=1/m0; % set carrying capacity to match baseline mutation rate
E=0.1; % effective baseline mutation rate

%% Single realization
% Experiment how stochastic trajectories look like
tic
dose=104.5;
% Initialize
t=0; % elapsed time
mut=mu(0,dose); % mutation rate for optimal therapy
S=E*K; % No. sensitive cells at the beginning
R=0; % No. resistant cells
mtimes=[]; % collect mutation times
extimes=[]; % collect stoch. extinction times
pops=[S;R]; 
update=0;
uptime=0.01; %  update interval

while t < Tf
    % Calculate event propensities
    bS=(betaS-(betaS-thetaS)*(S+R)/K)*S; 
    dS=(deltaS+(thetaS-deltaS)*(S+R)/K+d(0,dose))*S; 
    mS=mut*S; 
    bR=(betaR-(betaR-thetaR)*(S+R)/K)*R;
    dR=(deltaR+(thetaR-deltaR)*(S+R)/K)*R; 
    P=bS+dS+mS+bR+dR; % total propensity
    
    % Quit, if all cells are extinct
    if P==0
        break
    end
    
    % Determine the time of the next event  
    tau=-log(rand)/P;
    t=t+tau; % update elapsed time
    
    % Choose events based on ur ~ U(0,1)
    ur=rand; % 
    if ur < bS/P
        S=S+1;
    elseif ur < (bS+dS)/P
        S=S-1;
    elseif ur < (bS+dS+mS)/P
        S=S-1;
        R=R+1;
        % record only newly emerged mutations
        if R==1
            mtimes=[mtimes, t];  % record time of mutation events
        end
    elseif ur < (bS+dS+mS+bR)/P
        R=R+1;
    elseif ur < (bS+dS+mS+bR+dR)/P
        R=R-1;
        if R==0
            extimes=[extimes, t];  % record time of extinction events
        end
    end
    % Update pop sizes at update times
    if t > update + uptime
        pops=[pops(1,:), S; pops(2,:), R];
        update=t;
    end
end
toc
figure(1)
plot(linspace(0,t,length(pops(1,:))),pops(1,:),'linewidth',3)
hold on
plot(linspace(0,t,length(pops(2,:))),pops(2,:),'linewidth',3)
xlabel('Time')
ylabel('Population size')
title('Stochastic trajectories')
set(gca,'yscale','log')
set(gca,'FontSize', 18, 'FontWeight', 'bold')


%%
tic
for u=1:length(doses)
    dose=doses(u);
    finalpops=zeros(2,Nsim); % collect final pop. sizes for each sim. 
    for j=1:Nsim   
        % Initialize
        t=0; % elapsed time
        mut=mu(0,dose); % mutation rate for optimal therapy
        S=E*K; % No. sensitive cells at the beginning
        R=0; % No. resistant cells

        while t < Tf
            % Calculate event propensities
            bS=(betaS-(betaS-thetaS)*(S+R)/K)*S; 
            dS=(deltaS+(thetaS-deltaS)*(S+R)/K+d(0,dose))*S; 
            mS=mut*S; 
            bR=(betaR-(betaR-thetaR)*(S+R)/K)*R;
            dR=(deltaR+(thetaR-deltaR)*(S+R)/K)*R; 
            P=bS+dS+mS+bR+dR; % total propensity
            
            % Determine the time of the next event  
            tau=-log(rand)/P;
            t=t+tau; % update elapsed time
    
            % Choose events based on ur ~ U(0,1)
            ur=rand; % 
            if ur < bS/P
                S=S+1;
            elseif ur < (bS+dS)/P
                S=S-1;
            elseif ur < (bS+dS+mS)/P
                S=S-1;
                R=R+1;
            elseif ur < (bS+dS+mS+bR)/P
                R=R+1;
            elseif ur < (bS+dS+mS+bR+dR)/P
                R=R-1;
            end
        end
        finalpops(:,j)=[S,R];
    end
    finalS=[finalS, finalpops(1,:).']; 
    finalR=[finalR, finalpops(2,:).']; 
    finalN=[finalN, sum(finalpops,1).'];
    disp([num2str(u), ' out of ', num2str(length(doses)), ' doses complete'])
end
toc

% use plotFig4.m to visualize the results
