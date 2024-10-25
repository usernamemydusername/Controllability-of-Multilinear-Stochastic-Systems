%% 3. Stochasitc Noises
% A1...Ak: n*n, diagonalizable
% B1...Bk: n*m
% E1...Ek: n*l
% An = Ak kron Ak-1 kron ... A1
% Bn = Bk kron Bk-1 kron ... B1
% En = Ek kron Ek-1 kron ... E1

%
% we have three candidate methods
% a) C = [B, AB, ..., A^(n-1)B, E, AE, ..., A^(n-1)E] is full rank or not, compute time for it.
% b) A^p B = (A_k^p B_k) kron ... kron (A_1^p B_1), and 
%    A^p E = (A_k^p E_k) kron ... kron (A_1^p E_1), test whether LHS =
%    RHS, and compute time needed to check C's rank if we use RHS to construct
%    C
% c) S1 = ctrb(A1, B1), ... , Sk = ctrb(Ak, Bk), 
%    O1 = ctrb(A1, E1), ... , Ok = ctrb(Ak, Ek), get time needed to check
%    all Sp's rank and Op's rank --> controllable if either all Sp's or all Op's are full rank(Coro 9)
%% GO: we have m = n-1, l = n-1
clear;
clc;

%n_vals = 2;
%k_vals = 2:10;
n_vals = 2:10;
k_vals = 2;
num_trials = 100; % # trials for each n and k to average time
time_a = zeros(length(n_vals), length(k_vals));
time_b = zeros(length(n_vals), length(k_vals));
time_c = zeros(length(n_vals), length(k_vals));

for ni = 1:length(n_vals)
    for ki = 1:length(k_vals)
        n = n_vals(ni);
        k = k_vals(ki);
        fprintf('Running for n=%d, k=%d\n', n, k);
        
        time_a_trials = zeros(num_trials, 1);
        time_b_trials = zeros(num_trials, 1);
        time_c_trials = zeros(num_trials, 1);
        
        for trial = 1:num_trials
            disp(trial);
            A = cell(1, k);
            B = cell(1, k);
            E = cell(1, k);
            for p = 1:k
                [V, ~] = eig(rand(n));
                D = diag(rand(1, n));   % Ensure distinct eigenvalues
                A{p} = V * D / V;       % diagonalizable matrix A{p}
                B{p} = rand(n, n-1); % n x m, where m = n-1
                E{p} = rand(n, n-1); % n x l, where l = n-1
            end
            
            % a) Compute controllability matrix C directly
            An = A{k};
            Bn = B{k};
            En = E{k};
            for p = k-1:-1:1
                An = kron(An, A{p});
                Bn = kron(Bn, B{p});
                En = kron(En, E{p});
            end
            % C = [B, AB, ..., A^(n-1)B, E, AE, ..., A^(n-1)E]
            tic;
            C = [ctrb(An, Bn), ctrb(An, En)];
            
            % Test if C is full rank
            [~, S_mat, ~] = svd(C);
            S_diag_mat = diag(S_mat);
            last_nonzero_S_diag_mat = find(S_diag_mat >= 1e-15*S_diag_mat(1), 1, 'last');
            S_mat = numel(S_diag_mat(1:last_nonzero_S_diag_mat));
            rank_C = S_mat;
            time_a_trials(trial) = toc;
            
            % b) Use recursive relation A^p B = (A_k^p B_k) kron ... kron (A_1^p B_1), and A^p E = (A_k^p E_k) kron ... kron (A_1^p E_1),
            % sigma_p: kron of A^p B
            % sigma_p2: kron of A^p E -- same as omega_p in the paper
            tic;
            C_recursive_B = [];
            C_recursive_E = [];
            for i = 0:n-1
                sigma_p = 1;
                sigma_p2 = 1;
                for p = k:-1:1
                    sigma_p = kron(sigma_p, A{p}^i * B{p});
                    sigma_p2 = kron(sigma_p2, A{p}^i * E{p});
                end
                C_recursive_B = [C_recursive_B, sigma_p];  % Append to C
                C_recursive_E = [C_recursive_E, sigma_p2];
            end
            C_recursive = [C_recursive_B, C_recursive_E];
            
            % Compute rank from recursive C
            [~, S_mat, ~] = svd(C_recursive);
            S_diag_mat = diag(S_mat);
            last_nonzero_S_diag_mat = find(S_diag_mat >= 1e-15*S_diag_mat(1), 1, 'last');
            S_mat = numel(S_diag_mat(1:last_nonzero_S_diag_mat));
            rank_C_recursive = S_mat;
            time_b_trials(trial) = toc;
            
            % Now check if sigma_p equals An^i * Bn
            for i = 0:n-1
                sigma_p = 1;
                sigma_p2 = 1;
                for p = k:-1:1
                    sigma_p = kron(sigma_p, A{p}^i * B{p});
                    sigma_p2 = kron(sigma_p2, A{p}^i * E{p});
                end
                
                if norm(sigma_p - (An^i * Bn), 'fro') > eps
                    warning('Mismatch: sigma_p is not equal to An^i * Bn for p=%d, i=%d', k, i);
                end
                if norm(sigma_p2 - (An^i * En), 'fro') > eps
                    warning('Mismatch: sigma_p2 is not equal to An^i * En for p=%d, i=%d', k, i);
                end
            end
            
            % c) Using ctrb for each (Ap, Bp) and (Ap, Ep) individually
            tic;
            ranks_Sp = zeros(k, 1);
            ranks_Op = zeros(k, 1);
            for p = 1:k
                Cp = ctrb(A{p}, B{p});
                Op = ctrb(A{p}, E{p});
                [~, S_mat, ~] = svd(Cp);
                [~, S_mat2, ~] = svd(Op);
                S_diag_mat = diag(S_mat);
                last_nonzero_S_diag_mat = find(S_diag_mat >= 1e-15*S_diag_mat(1), 1, 'last');
                S_mat = numel(S_diag_mat(1:last_nonzero_S_diag_mat));
                ranks_Sp(p) = S_mat;
                
                S_diag_mat2 = diag(S_mat2);
                last_nonzero_S_diag_mat2 = find(S_diag_mat2 >= 1e-15*S_diag_mat2(1), 1, 'last');
                S_mat2 = numel(S_diag_mat2(1:last_nonzero_S_diag_mat2));
                ranks_Op(p) = S_mat2;
            end
            time_c_trials(trial) = toc;
        end
        
        time_a(ni, ki) = mean(time_a_trials);
        time_b(ni, ki) = mean(time_b_trials);
        time_c(ni, ki) = mean(time_c_trials);
        
        % Output the time taken
        fprintf('Time for method a: %.4f seconds\n', time_a(ni, ki));
        fprintf('Time for method b: %.4f seconds\n', time_b(ni, ki));
        fprintf('Time for method c: %.4f seconds\n', time_c(ni, ki));
    end
end

%% fix k, increase n
figure(3)
set(gcf, 'Position', [100, 100, 600, 300]); 
plot(n_vals, log(time_a),'-o', 'DisplayName','Method a', 'LineWidth', 2,  'Color', [1, 0, 0]);
hold on;
plot(n_vals, log(time_b),'-s', 'DisplayName','Method b', 'LineWidth', 2,  'Color', [0, 0, 1]);
hold on;
plot(n_vals, log(time_c),'-^', 'DisplayName','Method c', 'LineWidth', 2,  'Color', [0, 0.8, 0]);
hold off;

title('Log Computational Time for k = 2', 'FontSize', 14)
xlabel('Matrix Size n', 'FontWeight', 'bold', 'FontSize', 20)
ylabel('Time (s)', 'FontWeight', 'bold',  'FontSize', 20)
set(gca, 'FontSize', 15)

legend('Location', 'northwest','FontSize', 14); 

%% fix n, increase k
figure(3)
set(gcf, 'Position', [100, 100, 600, 300]); 
plot(k_vals, log(time_a),'-o', 'DisplayName','Method a', 'LineWidth', 2, 'Color', [1, 0, 0]);
hold on;
plot(k_vals, log(time_b),'-s', 'DisplayName','Method b', 'LineWidth', 2, 'Color', [0, 0, 1]);
hold on;
plot(k_vals, log(time_c),'-^', 'DisplayName','Method c', 'LineWidth', 2, 'Color', [0, 0.8, 0]);
hold off;

title('Log Computational Time for n = 2', 'FontSize', 14)
xlabel('k', 'FontWeight', 'bold', 'FontSize', 20)
ylabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 20)
set(gca, 'FontSize', 15);


legend('Location', 'northwest','FontSize', 14); 
%% Plot times as n and k increases
figure;
for ki = 1:length(k_vals)
    plot(n_vals, time_a(:, ki), '-o', 'DisplayName', ['Method a (k=' num2str(k_vals(ki)) ')']);
    hold on;
    plot(n_vals, time_b(:, ki), '-s', 'DisplayName', ['Method b (k=' num2str(k_vals(ki)) ')']);
    plot(n_vals, time_c(:, ki), '-^', 'DisplayName', ['Method c (k=' num2str(k_vals(ki)) ')']);
end
xlabel('n (size of Ap and Bp)');
ylabel('Time (seconds)');
title('Comparison of Computation Time for Methods a, b, and c');
legend;
grid on;

%%
% Number of subplots based on the number of k values
num_k_vals = length(k_vals);
figure;

for ki = 1:num_k_vals
    subplot(1, num_k_vals, ki);
    hold on;

    % method a
    plot(n_vals, time_a(:, ki), '-o', 'DisplayName', 'Method a');
    
    % method b
    plot(n_vals, time_b(:, ki), '-s', 'DisplayName', 'Method b');
    
    % method c
    plot(n_vals, time_c(:, ki), '-^', 'DisplayName', 'Method c');

    xlabel('n (size of A_p and B_p)');
    ylabel('Time (seconds)');
    title(['Comparison of Computation Time for Methods a, b, and c (k = ' num2str(k_vals(ki)) ')']);
    legend('show', 'Location', 'best');
    grid on; 
    hold off;
end

sgtitle('Computation Time Comparison Across Different k Values');

