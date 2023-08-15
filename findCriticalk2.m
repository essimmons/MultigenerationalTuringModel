function [k2,realEigsAp,k2crit,realEigCrit,w] = findCriticalk2(J,param,Dv)
        
    Du=param.D_u;
    k2 = logspace(-4,-1,1000); % this is k^2
    D = diag([Du;Dv;Du;Dv]);
    eigsAp = zeros(4,length(k2));
    
    for i=1:length(k2)
        eigsAp(:,i)=eigs(J-k2(i)*D);
    end
    
    realEigsAp=sort(real(eigsAp),1);
    [realEigCrit,I]=max(realEigsAp(4,:));
    k2crit=k2(I);
    k_crit=sqrt(k2crit);
    w=2*pi/k_crit;
    
end