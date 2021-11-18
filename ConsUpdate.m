function [P,Q,Pc,K,G,SIG,dL] = ConsUpdate(SIG,C,EPSI,Cr,Cc,Pc,K,G,M,nu)
% FUNCTION TO INTEGRATE STRESSES IN MODIFIED CAM-CLAY (MCC)
% THIS FILE IS CALLED ConsUpdate.m
% NOTE WE USE TENSORS AND VECTORS IN VOIGT NOTATION
% APPROXIMATE INETGRATION OF Pc
% AUTHOR: Enrique M. del Castillo

% PRESCRIBE TOLERANCE FOR N-R ALGO

FTOL=10^(-8);

% CALCULATE TRIAL STRESSES
for i = 1:6
    for j = 1:6
      SIG(i,1) = SIG(i,1) + C(i,j)*EPSI(j,1);
    end
end
    
% CALCULATE P-trial

P_tr = (SIG(1,1)+SIG(2,1)+SIG(3,1))/3;

% CALCULATE DEVIATORIC PART OF SIGMA AND OTHER INVARIANTS

ID=[ 1; 1; 1; 0; 0; 0];
for j = 1:6
    XI(j,1) = SIG(j,1) - P_tr*ID(j,1); % DEVIATORIC TRIAL STRESS TENSOR XI
end

Q_tr= sqrt(3/2)*sqrt(XI(1,1)*XI(1,1) + XI(2,1)*XI(2,1) + ...
    XI(3,1)*XI(3,1)+XI(4,1)*XI(4,1)*2.0+XI(5,1)*XI(5,1)*2 + ...
    XI(6,1)*XI(6,1)*2);

fprintf('P_tr = %f\n', P_tr); 
fprintf('Q_tr = %f\n', Q_tr); 
   
P_k=P_tr; Q_k=Q_tr; Pc_k=Pc; DeltaL_k=0; K_k =K; G_k=G; Pcn=Pc;

F_k=(Q_k)^2/M^2 + P_k*(P_k-Pc_k); %yield criterion using sigma_tr

% Check yielding
if F_k > 0

    % YES --> PLASTIC PHASE: RETURN MAPPING
    n_k = 0;
    while (abs(F_k)>=FTOL && n_k < 10)
        A= (-2*K_k*((DeltaL_k*(2*P_k-Pc_k)+Cc-Cr)^2)*(P_k-Pc_k/2))/...
            (8*K_k*((P_k-Pc_k/2)^2)*DeltaL_k^3+8*(K_k*(Cc-Cr)+P_k/2-Pc_k/4)*...
            (P_k-Pc_k/2)*DeltaL_k^2+2*(Cc-Cr)*(K*(Cc-Cr)+2*P_k-Pc_k-Pcn/2)*...
            DeltaL_k+(Cc-Cr)^2); %partial P/partial Delta Lambda

        B= (-Q_k)/(DeltaL_k+(M^2/(6*G_k))); %partial Q/partial Delta Lambda

        C= ((-2*P_k+Pc_k)*(Cc-Cr)*Pcn)/(8*K_k*((P_k-Pc_k/2)^2)*DeltaL_k^3+...
            8*(K_k*(Cc-Cr)+P_k/2-Pc_k/4)*(P_k-Pc_k/2)*DeltaL_k^2+2*(Cc-Cr)*...
            (K*(Cc-Cr)+2*P_k-Pc_k-Pcn/2)*DeltaL_k+(Cc-Cr)^2); %partial P_c/partial Delta Lambda

        Fp_k=(2*P_k-Pc_k)*A +((2*Q_k)/M^2)*B +(-P_k)*C; %F'(\Delta \Lamda)

        DeltaL_k = DeltaL_k-(F_k/Fp_k);

        % UPDATE P, Pc, and Q

        c=Pcn*(Cc-Cr)*(1+2*K_k*DeltaL_k);
        b=(Cc-Cr)*(1+2*K_k*DeltaL_k)+2*DeltaL_k*P_tr;

        Pc_ka=.5*((b/DeltaL_k)+sqrt((b/DeltaL_k)^2-(4*c)/DeltaL_k));
        Pc_kb=.5*((b/DeltaL_k)-sqrt((b/DeltaL_k)^2-(4*c)/DeltaL_k));

        if (Pc_ka < 0) %ENSURE NEGATIVE ROOT
            Pc_k=Pc_ka;
        else
            Pc_k=Pc_kb;
        end

        P_k=(P_tr+K_k*DeltaL_k*Pc_k)/(1+2*K_k*DeltaL_k);
        Q_k=Q_tr/(1+6*G_k*DeltaL_k/M^2);

        F_k=(Q_k)^2/M^2 + P_k*(P_k-Pc_k);
        fprintf('F = %f\n', F_k);
        
        n_k = n_k + 1;
    end    
    % NO --> ELASTIC STEP: TRIAL VALUES ARE FINAL VALUES
end

P=P_k; % FINAL VALUES AFTER EXITING LOOP
Q=Q_k;
Pc=Pc_k;
dL=DeltaL_k;

K=-P/Cr; % CALCULATE NEW K, G USING n+1 UPDATED P
G=K*(3*(1-2*nu))/(2*(1+nu));

% CALCULATE FULL STRESS RATE TENSOR

% Check to see if there is deviatoric stress. Otherwise division by zero.

quotient=sqrt(XI(1,1)*XI(1,1) + XI(2,1)*XI(2,1) + XI(3,1)*XI(3,1)+...
        XI(4,1)*XI(4,1)*2.0+XI(5,1)*XI(5,1)*2+XI(6,1)*XI(6,1)*2);

if quotient > 0
    nhat=XI/quotient;
else
    nhat = zeros(6,1);
end

SIG=P*ID+sqrt(2/3)*Q*nhat;
    
return




