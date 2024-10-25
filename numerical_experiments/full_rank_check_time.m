%% 2. ctrb
% A1...Ak: n*n, diagonalizable
% B1...Bk: n*m
% An = Ak kron Ak-1 kron ... A1
% Bn = Bk kron Bk-1 kron ... B1
%
% we have three candidate methods
% a) C = [B, AB, ..., A^(n-1)B] is full rank or not, compute time for it.
% b) A^p B = (A_k^p B_k) kron ... kron (A_1^p B_1), test whether LHS =
%    RHS, and compute time needed to check C's rank if we use RHS to construct
%    C
% c) C1 = ctrb(A1, B1), ... , Ck = ctrb(Ak, Bk), get time needed to check
%    all Cp's rank
A1 = rand(2,2);
ctrb(A1, rand(2,3))

%% GO: we have m = n-1
clear;
clc;

%n_vals = 2:10;
%k_vals = 2;
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

        A = cell(1, k);
        B = cell(1, k);
        
        for p = 1:k
            [V, ~] = eig(rand(n));  
            D = diag(rand(1, n));   
            A{p} = V * D / V;       
            B{p} = rand(n, n-1); 
        end
        
        for trial = 1:num_trials
            disp(trial);
            % a) Compute controllability matrix C directly
            An = A{k};
            Bn = B{k};
            for p = k-1:-1:1
                An = kron(An, A{p});
                Bn = kron(Bn, B{p});
            end
            tic;
            C = ctrb(An, Bn);
            [~, S_mat, ~] = svd(C);
            S_diag_mat = diag(S_mat);
            last_nonzero_S_diag_mat = find(S_diag_mat >= 1e-15*S_diag_mat(1), 1, 'last');
            S_mat = numel(S_diag_mat(1:last_nonzero_S_diag_mat));
            rank_C = S_mat;
            time_a(ni, ki) = time_a(ni, ki) + toc;

            % b) Use recursive relation A^p B = (A_k^p B_k) kron ... kron (A_1^p B_1)
            tic;
            C_recursive = [];
            for i = 0:n-1
                sigma_p = 1;
                for p = k:-1:1
                    sigma_p = kron(sigma_p, A{p}^i * B{p});
                end
                C_recursive = [C_recursive, sigma_p];  
            end
            [~, S_mat, ~] = svd(C_recursive);
            S_diag_mat = diag(S_mat);
            last_nonzero_S_diag_mat = find(S_diag_mat >= 1e-15*S_diag_mat(1), 1, 'last');
            S_mat = numel(S_diag_mat(1:last_nonzero_S_diag_mat));
            rank_C_recursive = S_mat;
            time_b(ni, ki) = time_b(ni, ki) + toc;

            % Now check if sigma_p equals An^i * Bn
            for i = 0:n-1
                sigma_p = 1;
                for p = k:-1:1
                    sigma_p = kron(sigma_p, A{p}^i * B{p});
                end
                
                if norm(sigma_p - (An^i * Bn), 'fro') > eps
                    warning('Mismatch: sigma_p is not equal to An^i * Bn for p=%d, i=%d', k, i);
                end
            end
            
            % c) Using ctrb for each (Ap, Bp) individually
            tic;
            ranks_Cp = zeros(k, 1);
            for p = 1:k
                Cp = ctrb(A{p}, B{p});
                [~, S_mat, ~] = svd(Cp);
                S_diag_mat = diag(S_mat);
                last_nonzero_S_diag_mat = find(S_diag_mat >= 1e-15*S_diag_mat(1), 1, 'last');
                S_mat = numel(S_diag_mat(1:last_nonzero_S_diag_mat));
                ranks_Cp(p) = S_mat;
            end
            time_c(ni, ki) = time_c(ni, ki) + toc; 
        end
        
        % Average time for each method
        time_a(ni, ki) = time_a(ni, ki) / num_trials;
        time_b(ni, ki) = time_b(ni, ki) / num_trials;
        time_c(ni, ki) = time_c(ni, ki) / num_trials;

        % Output the time taken
        fprintf('Time for method a: %.4f seconds\n', time_a(ni, ki));
        fprintf('Time for method b: %.4f seconds\n', time_b(ni, ki));
        fprintf('Time for method c: %.4f seconds\n', time_c(ni, ki));
    end
end


%% fix k, increase n
figure(3)
set(gcf, 'Position', [100, 100, 600, 300]); 
plot(n_vals, log(time_a),'-o', 'DisplayName','Method a', 'LineWidth', 2, 'Color', [1, 0, 0]);
hold on;
plot(n_vals, log(time_b),'-s', 'DisplayName','Method b', 'LineWidth', 2, 'Color', [0, 0, 1]);
hold on;
plot(n_vals, log(time_c),'-^', 'DisplayName','Method c', 'LineWidth', 2, 'Color', [0, 0.8, 0]);
hold off;

title('Log Computational Time for k = 2', 'FontSize', 14)
xlabel('Matrix Size n', 'FontWeight', 'bold', 'FontSize', 20)
ylabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 20)
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
set(gca, 'FontSize', 15)

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




