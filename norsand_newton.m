%%%%%%%%%%%%%%%%%%
% NORSAND.m
%
% A simple implementation of the NorSand constitutive model using the
% cutting plane method. For constitutive model verification purposes only.
%
% Author: Ben Mullet
% Date: 12 November 2021

close all;
clc;
clear all;

p_i_tol = 1;                    % p_i loop tolerance
r_tol = 1;                      % Outer loop tolerance

% NorSand model parameters
lambda_tilde = 0.04;            % Compressibility
N = 0.4;                        % Yield function constant
N_bar = 0.2;                    % Plastic potential constant
h = 280;                        % Hardening coefficient
M = 1.45;                       % CSL slope

G = 2000e3;                     % Shear modulus (Pa)
nu = 0.3;                       % Poisson's ratio
K = 2*G*(1+nu)/(3*(1-2*nu));    % Bulk modulus (Pa)

p_i0 = -60e3;                    % Sets up preconsolidation stress ~ -130kPa
vc0 = 1.915;                    % Reference specific volume
v0 = 1.63;                      % Initial specific volume
alpha = 3.5;                    % Multiplier for max plastic dilatancy
beta = (1-N)/(1-N_bar);         % Non-associativity parameter
alpha_bar = alpha/beta;         % Multiplier for pstar calculation

% Struct is used for function params passing 
params.N = N;
params.N_bar = N_bar;
params.h = h;
params.M = M;
params.lambda_tilde = lambda_tilde;
params.G = G;
params.K = K;
params.vc0 = vc0;
params.v0 = v0;
params.alpha_bar = alpha_bar;
params.beta = beta;

% Failure function
eta_func = @(p, pi) (M/N)*(1 - (1-N)*(p./pi).^(N/(1-N)));
F_func = @(p, q, pi) q + eta_func(p, pi).*p;

figure
xlim([-200e3, 0])
ylim([0,100]*1e3)

hold on
fimplicit(@(p, q) F_func(p, q, p_i0))
plot(xlim, [0,0], '--r')

% Set up isotropic elastic stiffness matrix
C = get_stiffness(K,G);

% Helper functions
to_matrix_sigma = @(s) [s(1), s(6), s(5);
                        s(6), s(2), s(4);
                        s(5), s(4), s(3)];

to_voigt_epsilon = @(e) [e(1,1); e(2,2); e(3,3); 2*e(2,3); 2*e(1,3); 2*e(1,2)];



% get_q defined below

%%% TEST ELASTICITY CALCULATIONS %%%
eps_ = [1, 0, 0; 0, 0, 0; 0, 0, 0]; 
sig_ = hookes(eps_, 3, 2);

% should be     [5.6667         0         0
%                0         1.6667         0
%                0              0    1.6667]


eps_ = [0, 1, 0; 1, 0, 0; 0, 0, 0]; 
sig_ = hookes(eps_, 3, 2);

% should be    [0     4     0
%               4     0     0
%               0     0     0]

%%% END TEST %%%%

% Initial hydrostatic compression
eps0 = -0.007692*eye(3); % sets up p = -1e5 Pa
sigma0 = hookes(eps0, K, G);

plot_p_q(sigma0, -1);
drawnow

% compression in x direction
EPS = 1.5e-1*[-1, 0, 0; 0, 0.2, 0; 0, 0, 0.2]; % matches Ronnie's figure
%EPS = 1.5e-1*[0, 1, 0; 1, 0, 0; 0, 0, 0]; % matches Ronnie's figure

NUM_STEPS = 100; 

sigma = sigma0;         % stress
eps_tot = eps0;         % total strain
eps_e = eps_tot;        % total elastic strain
eps_p = zeros(3,3);     % total plastic strain

% p_i_n is the plastic internal variable value at the current time step
p_i_n = p_i0;

for i =1:NUM_STEPS
   delta_eps = EPS/NUM_STEPS;
   
   eps_tot = eps_tot + delta_eps;
   
   eps_e_tr = eps_e + delta_eps;        % elastic trial strain
   nhat = eps_e_tr / norm(eps_e_tr);
  
   delta_sigma_tr = hookes(delta_eps, K, G);
   sigma_tr = sigma + delta_sigma_tr;
   
   % elastic stress predictor
   %sigma_tr = hookes(eps_e_tr, params.K, params.G);
   
   p_tr = get_p(sigma_tr);
   q_tr = get_q(sigma_tr);
   
   p_i_tr = p_i_n;
   
   Fk = F_func(p_tr, q_tr, p_i_tr);
   
   if Fk < 0
       % No failure, pure elastic deformation
       sigma = sigma_tr;
       eps_e = eps_e_tr;
       p_i_n = p_i_tr;
       
   else
       % Iterate for plastic corrector etc as in steps 5-7, Box 2 of
       % Borja and Andrade (2005)
       
       dlambda = 0;
       
       
       v = params.v0*(1 + trace(eps_tot));
       
       eps_e_k = eps_e_tr; % Iterator for elastic strain
       sigma_k = sigma_tr; % Iterator for sigma
       
       [r, x] = get_residual(eps_e_k, eps_e_tr, dlambda, p_i_n, params);
       
       p_i = p_i_n;
      
       while norm(r) > r_tol
           % Step 5: plastic corrector
           % [done at the end of the loop]

           % Step 6: Update plastic internal variable
           p_i = p_i_n;
           r_pi = 100; % Nominal value to force loop at least once
           
           while r_pi > p_i_tol
               
               psi = get_psi(params, eps_tot, p_i);
               pstar = get_pstar(params, sigma_k, eps_tot, p_i);
               
               r_pi = p_i - p_i_n + dlambda*h*(p_i - pstar);
               
               coeff = params.lambda_tilde*params.alpha_bar*(1 - params.N)...
                   /(params.M - params.alpha_bar*psi*params.N);
               
               r_pi_prime = 1 + h*dlambda*(1 - coeff*(pstar/p_i));
               
               p_i = p_i - r_pi/r_pi_prime;
               
           end
           
           % Step 7: Discrete consistency condition
           [r, x] = get_residual(eps_e_k, eps_e_tr, dlambda, p_i, params);
           r_prime = get_jacobian(eps_tot, eps_e_k, p_i, dlambda, params);
           
           delta_x = r_prime \ (-r);
           x = x + delta_x;
           
           eps_e_v_k = x(1);
           eps_e_s_k = x(2);
           dlambda = x(3);
           
           eps_e_k = 1/3*eps_e_v_k*eye(3) + sqrt(3/2)*eps_e_s_k*nhat;
           sigma_k = hookes(eps_e_k, params.K, params.G);
           
       end
       
       p_i_n = p_i;
       sigma = sigma_k;
       eps_e = eps_e_k;
       Fk = F_from_sigma(sigma, p_i, params);
    
   end
   plot_pi(p_i_n, Fk, i)
   
   plot_p_q(sigma, Fk)
    
end

function r_prime = get_jacobian(eps_tot, eps_e, p_i, dlambda, params)

sigma = hookes(eps_e, params.K, params.G);
p = get_p(sigma);
q = get_q(sigma);

pstar = get_pstar(params, sigma, eps_tot, p_i);

psi = get_psi(params, eps_tot, p_i);

M = params.M;
N = params.N;

dFdp = (M/N)*(1 - (p/p_i)^(N/(1-N)));
dFdq = 1;

dFdpi = M*(p/p_i)^(N/(1-N));

coeff = params.lambda_tilde*params.alpha_bar*(1 - params.N)...
                   /(params.M - params.alpha_bar*psi*params.N);
               
c = 1 + params.h*dlambda*(1 - coeff*(pstar/p));

D11 = params.K;
D12 = 0;
D13 = 0;
D21 = 0;
D22 = 3*params.G;
D23 = 0;
D31 = (1/c)*params.h*dlambda*(pstar/p)*D11;
D32 = (1/c)*params.h*dlambda*(pstar/p)*D12;
D33 = -(1/c)*params.h*(p_i - pstar);

D = [D11, D12, D13;
     D21, D22, D23;
     D31, D32, D33;];
     

H1 = -1/(1-params.N)*(params.M/p)*(p/p_i)^(params.N/(1-N));
H2 = 0;
H3 = 1/(1-params.N)*(params.M/p)*(p/p_i)^(params.N/(1-params.N));

H = [H1, H2, H3];

G = H*D;

r11 = 1+dlambda*params.beta*G(1);
r12 = dlambda*params.beta*G(2);
r13 = params.beta*(dFdp + dlambda*G(3));
r21 = 0;
r22 = 1;
r23 = dFdq;
r31 = D11*dFdp + D21*dFdq + D31*dFdpi;
r32 = D12*dFdp + D22*dFdq + D32*dFdpi;
r33 = dFdpi;

r_prime = [r11, r12, r13;
           r21, r22, r23;
           r31, r32, r33];

end

function [r, x] = get_residual(eps_e, eps_e_tr, dlambda, p_i, params)
[eps_e_v, eps_e_s] = split_epsilon(eps_e);
[eps_e_v_tr, eps_e_s_tr] = split_epsilon(eps_e_tr);

M = params.M;
N = params.N;



sigma = hookes(eps_e, params.K, params.G);
p = get_p(sigma);
q = get_q(sigma);

dFdp = (M/N)*(1 - (p/p_i)^(N/(1-N)));
dFdq = 1;

r1 = eps_e_v - eps_e_v_tr + dlambda*params.beta*dFdp;
r2 = eps_e_s - eps_e_s_tr + dlambda*dFdq;
r3 = F_from_sigma(sigma, p_i, params);

r = [r1; r2; r3];
x = [eps_e_v; eps_e_s; dlambda];

end

function F = F_from_sigma(sigma, p_i, params)

q = get_q(sigma);
p = get_p(sigma);

M = params.M;
N = params.N;

eta = (M/N)*(1 - (1-N)*(p./p_i).^(N/(1-N)));

F = q + eta*p;


end


function [eps_e_v, eps_e_s] = split_epsilon(eps_e)
eps_e_v = trace(eps_e);
e_e = eps_e - 1/3*eps_e_v*eye(3);
eps_e_s = sqrt(2/3)*norm(e_e);

end

function plot_pi(pi, Fk, i)
    figure(2)
    hold on;
    
    if Fk > 0
        color = 'r';
    else
        color = 'b';
    end
    plot(i, pi,['o', color], 'MarkerFaceColor', color)
        
end

function fcq = multiply_F_C_Q(params, dFdsig, dQdsig)
cq = hookes(dQdsig, params.K, params.G);
fcq = sum(sum(dFdsig.*cq));

end

function H =  get_H(params, sigma, pi, pstar)
p = get_p(sigma);
H = params.M*params.h*(p/pi)^(1/(1-params.N))*(pi - pstar);
end

function pstar = get_pstar(params, sigma, eps_tot, pi)
psi = get_psi(params, eps_tot, pi);
%assert(psi < 0, 'Psi should be less than zero');

p = get_p(sigma);
pstar = p*(1 - params.alpha_bar*psi*params.N/params.M)^((params.N - 1)/params.N);


end

function psi = get_psi(params, eps_tot, pi)
v = params.v0*(1 + trace(eps_tot));
psi = v - params.vc0 + params.lambda_tilde*log(-pi);

% clc
% 
% disp(v)
% disp(trace(eps_tot))
% disp(params.v0)
% disp(params.vc0)
% disp(log(-pi))
% disp(params.lambda_tilde)
% disp(psi)
% 
% true

end

function dFdsig = get_dFdsig(params, sigma, eta)
p = get_p(sigma);
s = sigma - p*eye(3);
norm_s = sqrt(double_dot(s, s));

dFdsig = -1/3*(params.M - eta)/(1-params.N)*eye(3) + sqrt(3/2)*s/norm_s;

end

function dQdsig = get_dQdsig(params, sigma, eta)
p = get_p(sigma);
s = sigma - p*eye(3);
norm_s = sqrt(double_dot(s, s));

dQdsig = -1/3*(params.M - eta)/(1-params.N_bar)*eye(3) + sqrt(3/2)*s/norm_s;

end

function plot_p_q(sigma, Fk)
    figure(1)
    p = get_p(sigma);
    q = get_q(sigma);
    
    if Fk > 0
        color = 'r';
    else
        color = 'b';
    end
    plot(p,q,['o', color], 'MarkerFaceColor', color)
        
end



function p = get_p(sigma)
p =  1/3*trace(sigma);
end

function q = get_q(sigma)

p = get_p(sigma);
s = sigma - p*eye(3);
J2 = 1/2 *  double_dot(s, s);
q = sqrt(3*J2);

end

function s = double_dot(A, B)

s = sum(sum(A.*B));

end


function C = get_stiffness(K, G)

lambda = K - 2*G/3;
c11 = lambda + 2*G;
c12 = lambda;
c44 = G;

C = [c11, c12, c12,   0,   0,   0;
     c12, c11, c12,   0,   0,   0;
     c12, c12, c11,   0,   0,   0;
     0,     0,   0, c44,   0,   0;
     0,     0,   0,   0, c44,   0;
     0,     0,   0,   0,   0, c44 ];

end

function sigma = hookes(epsilon, K, G)
    assert(K > G, 'Bulk modulus must be greater than shear modulus')
    assert(all(size(epsilon) == [3,3]), 'Strain must be 3x3 matrix')
    
    lambda = K - 2*G/3;
    
    ekk = sum(sum(epsilon.*eye(3)));
    
    sigma = 2*G*epsilon + lambda*ekk*eye(3);
    

end