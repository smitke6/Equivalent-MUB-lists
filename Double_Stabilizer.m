clc;
clearvars;
addpath('QetLab');

%% Setup
d = 4;
% Construct the matrix F(x)
x = 0;
F = [1 1 1 1;
    1 1 -1 -1;
    1 -1 1i*exp(1i*x) -1i*exp(1i*x);
    1 -1 -1i*exp(1i*x) 1i*exp(1i*x)]/2;
F_dagger = F';

%% Find the intersection of the torus with the conjugated torus
% We wish to solve: F*T*F' belongs to the torus

% Choose two permutations in S^d
Symm_Group = flipud(perms(0:d-1));
num_perms = size(Symm_Group,1);
num_eqns = num_perms * (num_perms-1);

num_nontriv_sols = 0;

for ind_rho = 1:num_perms
    for ind_sig = 1:num_perms
        rho = Symm_Group(ind_rho,:);
        sig = Symm_Group(ind_sig,:);

        % Invert rho
        rho_inv = zeros(1, d);
        rho_inv(rho+1) = 0:d-1;

        % Construct the matrix of constraints
        A = zeros(d,num_eqns);
        counter = 1;
        for j = 0:d-1
            for m = 0:d-1
                if m == sig(j+1)
                    continue;
                end
              
                % Construct the next row of A to be the constraint
                % corresponding to [F'*U*F]_jm = 0
                for k = 1:d
                    A(counter,k) = F_dagger(j+1,rho_inv(k)+1)*F(k,m+1);
                end
                counter = counter+1;
            end
        end

        % Seek a nontrivial solution
        r = rank(A);
        if r == d
            continue;
        end
        disp('The permutations are:');
        disp('rho = ');
        disp(rho);
        disp('sig = ');
        disp(sig);

        disp('The rank of A is:');
        disp(rank(A));
        disp('The nontrivial solution is:');
        x = null(A);
        x = x / x(1);
        U = zeros(d);
        for i = 1:d
            U(i, rho(i)+1) = x(rho(i)+1);
        end
        disp(U);
        disp('##########################################################');
        num_nontriv_sols = num_nontriv_sols+1;
    end
end

disp('# of nontrivial solutions:');
disp(num_nontriv_sols);