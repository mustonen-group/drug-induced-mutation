% Help function for computeControlMap.m
% Tommi Mononen

function resvec=skellamVec(kmin,kmax,mu1,mu2)
    maxn=500; % notice that when carrying capacity is bigger, this has to be larger
              % otherwise results can be incorrect ('maxn' could be more automatic)
    
    kn=(-maxn+kmin):(maxn+kmax);
    knvec=poisspdf(kn,mu1);
    
    n=-maxn:maxn;
    nvec=poisspdf(n,mu2);

    resvec=zeros(1,kmax-kmin+1);
    
    for i=1:(kmax-kmin+1)
        resvec(i)=sum(knvec(i:((i-1)+length(n))).*nvec);
    end;
end