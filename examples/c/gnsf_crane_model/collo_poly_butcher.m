function [ Pdotvalues, Pvalues, A, b, c  ] = collo_poly_butcher( stages, dt)
% given the number of stages, and dt,
% the function returns a matrice Pdotvalues, in which the (i,j)th entry
% contains the value of \dot{p_i}(t_j), where p_i is the ith
% Lagrange-polynomial to the nodes t_j;
% Pvalues contains the values of the polynomials at the end time.

% [ Pdotvalues, Pvalues ] = collo_poly( 3, .1 )
if nargin<2
    dt = 0.1;
end
switch stages
    case 1
        legendre_zeros = [0; .5];
    case 2
        legendre_zeros = [0; 0.21132; 0.78867];
    case 3
        legendre_zeros = [0 ; 0.11270; 0.50000; 0.88729];
    case 4
        legendre_zeros = [0 ;0.06943; 0.33000; 0.66999; 0.93056];
end
% from NOC script;
% http://www.ams.org/journals/bull/1942-48-10/S0002-9904-1942-07771-8/S0002-9904-1942-07771-8.pdf
% Polynomials:
collo_times = legendre_zeros * dt;  % scale to [0, dt]
K = length(collo_times)-1; % number of stages

for i = 1:K+1
    y = zeros(K+1,1); 
    y(i) = 1;   % y contains the values of the ith lagrange polynom
                % at the collo_times
    % generate coeffs of ith lagrange polynom            
    P(i,:) = polyfit(collo_times, y, K);
    % generate coeffs of derived lagrange polynom
    Pdot (i,:) = polyder(P(i,:));
    Pdotvalues(:,i) = polyval( Pdot(i,:), collo_times);
    
    Pvalues(i) = polyval( P(i,:), dt);
end
if  1 % generate butcher tableau
    A = zeros(K);
    b = zeros(K,1);
    c = legendre_zeros(2:end);
    % TODO: not sure if 0 is a collocation time...
    % compare: impl. integrators: P_{k,i}= product{j=0:K}
    % versus dissertation: L_i = product{j=1:K}
    for i = 1:K
        y = zeros(K,1); 
        y(i) = 1;   % y contains the values of the ith lagrange polynom
                % at the legendre_zeros
        l(i,:) = polyfit(legendre_zeros(2:end), y,K-1);
        L(i,:) = polyint( l(i,:));
        b(i) = diff(polyval( L(i,:),[0,1]));
    end
    for i = 1:K
        for j = 1 : K
            A(i,j) = diff(polyval( L(j,:), [0,c(i)])); 
        end
    end
end

end