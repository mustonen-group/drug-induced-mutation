%%% Plotting mutation intensities
%%% Teemu Kuosmanen

%%% IMPORTANT: Run first defineModel.m 

% vector of doses for plotting
doses=[0,40,70,100,200,300];
leg=strings(1,length(doses));
S(1)=0.6; % set lower initial condition to allow growth to K
cols=lines(6); % choose colormap which is also visible in print
% Solve S and R for each dose
for j=1:length(doses)
    dose=doses(j)*ones(1,n);
    leg(j)=['u = ', num2str(doses(j))];
    for i=2:n
        tk=T(i-1);
        sk=S(i-1);
        rk=R(i-1);
        uk=dose(i-1);
        % Solve S
        k1=dt*fS(tk,sk,rk,uk);
        k2=dt*fS(tk+k1/2,sk+k1/2,rk+k1/2,uk+k1/2);
        k3=dt*fS(tk+k2/2,sk+k2/2,rk+k2/2,uk+k2/2);
        k4=dt*fS(tk+k3,sk+k3,rk+k3,uk+k3);
        S(i)=sk+1/6*(k1+2*(k2+k3)+k4);
        % Solve R
        k1=dt*fR(tk,sk,rk,uk);
        k2=dt*fR(tk+k1/2,sk+k1/2,rk+k1/2,uk+k1/2);
        k3=dt*fR(tk+k2/2,sk+k2/2,rk+k2/2,uk+k2/2);
        k4=dt*fR(tk+k3,sk+k3,rk+k3,uk+k3);
        R(i)=rk+1/6*(k1+2*(k2+k3)+k4);
    end
    figure(3) % Fig2
    plot(T,S.*mu(0,dose),'linewidth',4,'color',cols(j,:))
    hold on
    title('Mutation intensities','fontsize',20)
    xlabel('Time')
    ylabel('Probability of rescue mutation')
    xlim([0,20])
    set(gca,'FontSize', 18, 'FontWeight', 'bold')
end
legend(leg)