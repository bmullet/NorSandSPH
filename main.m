% MAIN PROGRAM
%
% PROGRAM TO INTEGRATE STRESSES WITH NORSAND CONSTITUTIVE MODEL
% USING RETURN MAPPING ALGORITHM FOR NON-ASSOCIATIVE FLOW RULE 
%
% MODIFICATION OF CODE ORIGINALLY BY ENRIQUE  Enrique M. del Castillo
%
% AUTHOR: Ben Mullet

clear vars
clearvars
clc
format longEng

% INITIALIZE ARRAYS FOR OUTPUT

p_array=[0.0];
q_array=[0.0];
pi_array=[0.0];
eps_array=[0];
eps_p_acc=[0];
sig_array=[0];
   
% DEFINE CONSTANTS AND PARAMETERS:
% INITIALIZE SIG, K, G, AND C

  for i = 1:6
    SIG(i,1) = 0.0;
  end

nu=0.3; % Poisson's ratio
e=2.0; % void ratio
lambda=.355; % virgin compression index
kappa=.0477; % recompression index
M=1.45; % CSL slope
v=(1+e)/(lambda-kappa); 
Cr=kappa/(1+e);
Cc=Cr+(1/v);

Cc = 0.1183;
Cr = 0.0159;

Pc=-98; % initial preconsolidation pressure
Pc0 = Pc;

K=1.1581e3; % initial K to get K_n+1 going since K depends on pressure which is initially zero.
G=K*(3*(1-2*nu))/(2*(1+nu)); % Shear Modulus

K = 16800000.0;
G = 7753846.154;

E=9*K*G/(3*K+G); % Young's modulus
C=(E/((1+nu)*(1-2*nu)))*[1-nu nu nu 0 0 0;nu 1-nu nu 0 0 0;nu nu 1-nu 0 0 0;
 0 0 0 (1-2*nu)/2 0 0;0 0 0 0 (1-2*nu)/2 0;0 0 0 0 0 (1-2*nu)/2]; % Elastic Tangent Matrix

% PRESCRIBE STRAIN HISTORY AND STRAIN INCREMENTS 

exx=0.0; eyy=0.0;ezz=0.0;gxy=0.099739926;gyz=0;gxz=0; % INPUT HERE THE TOTAL PRESCRIBED STRAINS
EPSI_tot=[exx; eyy; ezz; gxy; gyz; gxz];
numsteps = 500; % NUMBER OF INCREMENTS = numsteps
n=1;

% Try to prescribe a stress path with hydrostatic phase, followed by purely
% deviatoric. Should follow yield surface and then the C.S. line.
EPSI1=EPSI_tot/(numsteps/n);
% EPSI2=[0; 0; 0; 0.1; 0; 0]/((n-1)*numsteps/n);
counter=1;

% INTEGRATE STRESSES AND CONSTITUTIVE UPDATE

EPSI = EPSI1;

for n=1:numsteps  
    disp('Step')
    disp(n)
        
%     if (n <= numsteps)
%         EPSI = EPSI1;
%     else
%         EPSI = EPSI2;
%     end
    
    [P,Q,Pc,K,G,SIG,dL]=ConsUpdate(SIG,C,EPSI,Cr,Cc,Pc,K,G,M,nu); % call constitutive update function

    E=9*K*G/(3*K+G);
    C=(E/((1+nu)*(1-2*nu)))*[1-nu nu nu 0 0 0;nu 1-nu nu 0 0 0;nu nu ...
        1-nu 0 0 0;
     0 0 0 (1-2*nu)/2 0 0;0 0 0 0 (1-2*nu)/2 0;0 0 0 0 0 (1-2*nu)/2];   

    % RECORD OUTPUT

    p_array=[p_array;P];
    q_array=[q_array;Q];
    pc_array=[pc_array;Pc];
    eps_p_acc=[eps_p_acc;eps_p_acc(counter)+dL];
    eps_array=[eps_array;EPSI(4,1)+eps_array(counter)];
    counter=counter+1;
    sig_array=[sig_array;SIG(4,1)];

end
%% PLOTTING

%Figure of stress history in p-q space
figure 
%plot(p_array,q_array,'-r','LineWidth',2)
scatter(p_array,q_array,'or','filled')
set(gca,'FontSize',16,'TickLabelInterpreter','latex')
ylabel('$q$','Interpreter','latex','FontSize',16);
xlabel('$p$','Interpreter','latex','FontSize',16)
title('p-q space','Interpreter','latex','FontSize',16)

% plot initial yield surface
hold on
f=@(x,y) ((y.^2)/M^2) + x.^2 - x*Pc0;
fimplicit(f,'-k','LineWidth',2)
hold on

%plot CSL 
f2=@(x,y) -M.*x-y;
fimplicit(f2,'--k','LineWidth',2)
axis equal

% plot changing yield surfaces
% for n=1:numsteps+1
%     if n%100 == 0
%         pc=pc_array(n);
%         hold on
%         fimplicit(@(x,y) ((y.^2)/M^2) + x.^2- x*pc)
%     end
% end

%Figure of shear stress vs strain
figure
plot(eps_array,sig_array)
ylabel('$\sigma_{xy}$','Interpreter','latex','FontSize',16);
xlabel('$\varepsilon_{xy}$','Interpreter','latex','FontSize',16)
title('Shear Stress vs. Strain','Interpreter','latex','FontSize',16)

%Figure of pc vs accumulated plastic strain
figure
plot(eps_p_acc,pc_array)
ylabel('$p_c$','Interpreter','latex','FontSize',16);
xlabel('$\Delta\bar{\epsilon}^p$','Interpreter','latex','FontSize',16)
title('Preconsolidation Stress vs. Accum. Plastic Strain','Interpreter',...
    'latex','FontSize',16)

