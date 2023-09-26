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
            g1_u1 g1_v1 g1_u2 g1_v2;
            f2_u1 f2_v1 f2_u2 f2_v2;            
            g2_u1 g2_v1 g2_u2 g2_v2];
end