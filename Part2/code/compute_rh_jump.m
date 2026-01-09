function [P2, rho2, T2, V2] = compute_rh_jump(y,sp,P1,T1,Vs)

%this function solves the Rankine-Hugniot equation to find the flow properties behind the shock

R_mix = FUN_get_mix_thermo_prop(sp,y,T1,'R_sp');
rho1 = P1 / (R_mix * T1);
h1 = FUN_get_mix_thermo_prop(sp,y,T1,'Cp')*T1;

rho2 = rho1*10; %initial guess
T2=T1; 

tolerance = 10^-10;
R1=100;
R2=100;
R3=100;
n=0;nmax=100;

while (R1>tolerance || R2>tolerance || R3>tolerance)
    V2 = (rho1/rho2)*Vs;
    P2 = P1 + rho1*Vs^2 * (1 - (rho1/rho2));
    h2 = h1 + 0.5*Vs^2 * (1 - (rho1/rho2)^2);

    Cp2 = FUN_get_mix_thermo_prop(sp,y,T2,'Cp');
    % T2 = T2 - ((h2 - h1 - 0.5*Vs^2 *(1 - (rho1/rho2)^2)))/Cp2;
    T2 = h2/Cp2;
    rho2 = P2/(FUN_get_mix_thermo_prop(sp,y,T2,'R_sp')*T2);

    R1 = sqrt(((rho1*Vs - rho2*V2)/(rho1*Vs))^2);
    R2 = sqrt(((P1 + rho1*Vs^2 - (P2 + rho2*V2^2))/(P1 + rho1*Vs^2))^2);
    R3 = sqrt(((h1 + 0.5*Vs^2 - (h2 + 0.5*V2^2))/(h1 + 0.5*Vs^2))^2);
    n=n+1;
    if n>nmax
        error('exceeded max iterations')
    end
end

fprintf('\n-----Final Residual for Rankine-Hugniot calculations-----')
fprintf('\nR1 = %d',R1)
fprintf('\nR2 = %d',R2)
fprintf('\nR3 = %d',R3)
fprintf('\n# of iterations = %d\n',n)
end