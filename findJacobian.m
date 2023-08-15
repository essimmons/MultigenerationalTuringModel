function [J_ss,u1_ss,u2_ss,v1_ss,v2_ss] = findJacobian(param)
    
    [u1_ss,u2_ss,v1_ss,v2_ss] = findSteadyState(param);

    bu1 = param.beta_u1; gu1 = param.gamma_u1;
    bu2 = param.beta_u2; gu2 = param.gamma_u2;
    bv1 = param.beta_v1; gv1 = param.gamma_v1;
    bv2 = param.beta_v2; gv2 = param.gamma_v2;
    ku1u1 = param.ku1u1; ku1u2 = param.ku1u2; ku1v1 = param.ku1v1; ku1v2 = param.ku1v2;
    ku2u1 = param.ku2u1; ku2u2 = param.ku2u2; ku2v1 = param.ku2v1; ku2v2 = param.ku2v2;
    kv1u1 = param.kv1u1; kv1u2 = param.kv1u2; kv1v1 = param.kv1v1; kv1v2 = param.kv1v2;
    kv2u1 = param.kv2u1; kv2u2 = param.kv2u2; kv2v1 = param.kv2v1; kv2v2 = param.kv2v2;
    
    U_u1 = ku1u1*u1_ss+ku1u2*u2_ss;
    U_u2 = ku2u1*u1_ss+ku2u2*u2_ss;
    U_v1 = kv1u1*u1_ss+kv1u2*u2_ss;
    U_v2 = kv2u1*u1_ss+kv2u2*u2_ss;
    V_u1 = ku1v1*v1_ss+ku1v2*v2_ss;
    V_u2 = ku2v1*v1_ss+ku2v2*v2_ss;
    V_v1 = kv1v1*v1_ss+kv1v2*v2_ss;
    V_v2 = kv2v1*v1_ss+kv2v2*v2_ss;
    
    X_u1 = 1+U_u1+V_u1;
    X_u2 = 1+U_u2+V_u2;
    X_v1 = 1+U_v1+V_v1;
    X_v2 = 1+U_v2+V_v2;
    
    f1_u1 = -bu1+(2*gu1*ku1u1*U_u1*(1+V_u1))/(X_u1^3);
    f1_v1 = (-2*gu1*ku1v1*U_u1^2)/(X_u1^3);
    f1_u2 = (2*gu1*ku1u2*U_u1*(1+V_u1))/(X_u1^3);
    f1_v2 = (-2*gu1*ku1v2*U_u1^2)/(X_u1^3);
    
    g1_u1 = (2*gv1*kv1u1*U_v1*(1+V_v1))/(X_v1^3);
    g1_v1 = -bv1+(-2*gv1*kv1v1*U_v1^2)/(X_v1^3);
    g1_u2 = (2*gv1*kv1u2*U_v1*(1+V_v1))/(X_v1^3);
    g1_v2 = (-2*gv1*kv1v2*U_v1^2)/(X_v1^3);
    
    f2_u1 = (2*gu2*ku2u1*U_u2*(1+V_u2))/(X_u2^3);
    f2_v1 = (-2*gu2*ku2v1*U_u2^2)/(X_u2^3);
    f2_u2 = -bu2+(2*gu2*ku2u2*U_u2*(1+V_u2))/(X_u2^3);
    f2_v2 = (-2*gu2*ku2v2*U_u2^2)/(X_u2^3);
    
    g2_u1 = (2*gv2*kv2u1*U_v2*(1+V_v2))/(X_v2^3);
    g2_v1 = (-2*gv2*kv2v1*U_v2^2)/(X_v2^3);
    g2_u2 = (2*gv2*kv2u2*U_v2*(1+V_v2))/(X_v2^3);
    g2_v2 = -bv2+(-2*gv2*kv2v2*U_v2^2)/(X_v2^3);
    
    J_ss = [f1_u1 f1_v1 f1_u2 f1_v2;
           g1_u1 g1_v1 g1_u2 g1_v2;           f2_u1 f2_v1 f2_u2 f2_v2;            
            g2_u1 g2_v1 g2_u2 g2_v2];
end

function [u1_ss,u2_ss,v1_ss,v2_ss] = findSteadyState(param)
    dt=0.1; %step size to be used
    N=1000000; %max steps to run
    %t=(0:dt:(N-1)*dt);
    u1=zeros(1,N); %create empty vector for u
    u1(1)=0.1; %initial value u
    u2=zeros(1,N); %create empty vector for u
    u2(1)=0.1; %initial value u
    v1=zeros(1,N); %create empty vector for v
    v1(1)=0.1; %initial value v
    v2=zeros(1,N); %create empty vector for v
    v2(1)=0.1; %initial value v
    epsilon = 10e-15;
    
    for n=1:N-1
        [f1,f2,g1,g2] = findUdot(u1(n),u2(n),v1(n),v2(n),param);
        u1(n+1)=u1(n)+dt*f1;
        u2(n+1)=u2(n)+dt*f2;
        v1(n+1)=v1(n)+dt*g1;
        v2(n+1)=v2(n)+dt*g2;
        if abs(f1) < abs(epsilon*u1(n)) && abs(f2) < abs(epsilon*u2(n)) && abs(g1) < abs(epsilon*v1(n)) && abs(g2) < abs(epsilon*v2(n))
            break
        end
    end 
    
    u1_ss=u1(n); u2_ss=u2(n); v1_ss=v1(n); v2_ss=v2(n);
    
    %figure(1)
    %plot(t(1:n),u1(1:n),t(1:n),v1(1:n))
    %figure(2)
    %plot(t(1:n),u2(1:n),t(1:n),v2(1:n))
end

function [f1,f2,g1,g2] = findUdot(u1,u2,v1,v2,param)
    au1 = param.alpha_u1; bu1 = param.beta_u1; gu1 = param.gamma_u1;
    au2 = param.alpha_u2; bu2 = param.beta_u2; gu2 = param.gamma_u2;
    av1 = param.alpha_v1; bv1 = param.beta_v1; gv1 = param.gamma_v1;
    av2 = param.alpha_v2; bv2 = param.beta_v2; gv2 = param.gamma_v2;
    ku1u1 = param.ku1u1; ku1u2 = param.ku1u2; ku1v1 = param.ku1v1; ku1v2 = param.ku1v2;
    ku2u1 = param.ku2u1; ku2u2 = param.ku2u2; ku2v1 = param.ku2v1; ku2v2 = param.ku2v2;
    kv1u1 = param.kv1u1; kv1u2 = param.kv1u2; kv1v1 = param.kv1v1; kv1v2 = param.kv1v2;
    kv2u1 = param.kv2u1; kv2u2 = param.kv2u2; kv2v1 = param.kv2v1; kv2v2 = param.kv2v2;
    
    
    f1 = -bu1*u1+au1+gu1*((ku1u1*u1+ku1u2*u2)/(1+ku1u1*u1+ku1u2*u2+ku1v1*v1+ku1v2*v2))^2;
    f2 = -bu2*u2+au2+gu2*((ku2u1*u1+ku2u2*u2)/(1+ku2u1*u1+ku2u2*u2+ku2v1*v1+ku2v2*v2))^2;
    g1 = -bv1*v1+av1+gv1*((kv1u1*u1+kv1u2*u2)/(1+kv1u1*u1+kv1u2*u2+kv1v1*v1+kv1v2*v2))^2;
    g2 = -bv2*v2+av2+gv2*((kv2u1*u1+kv2u2*u2)/(1+kv2u1*u1+kv2u2*u2+kv2v1*v1+kv2v2*v2))^2;
end