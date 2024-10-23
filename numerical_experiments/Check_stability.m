%% 1. check time of computing eig(A) vs (\sum_p eig(Ap)) + lambda1 kron lambda2 kron ... kron lambdak
%--------- prop 1 vs corollary 2
% n: each Ap is of size n*n
% k: there are k Ap's

n_values = [2,3,4,5,6,7,8,9,10,11,12,13,14,15];
k_values = [3];
tolerance = 1e-14;

% n*k matrix for time: each row shows time as k increases; each col shows time as n increases
time_eig_An_array = zeros(length(n_values), length(k_values)); 
time_eig_kron_array = zeros(length(n_values), length(k_values));

for ni = 1:length(n_values)
    n = n_values(ni);
    for ki = 1:length(k_values)
        k = k_values(ki);
        
        A = cell(1, k);
        for p = 1:k
            A{p} = rand(n, n);
        end
        
        % eigenvalue for An = Ak kron A_{k-1} kron... kron A_1
        An = A{k};
        for p = (k-1):-1:1  
            An = kron(An, A{p});
        end
        
        % Measure the time for eig(An): Proposition 1
        tic;
        eig_An = sort(eig(An));
        time_eig_An = toc;
        time_eig_An_array(ni, ki) = time_eig_An;
        
        % Kronecker product of eigenvalues from individual matrices
        
        % Measure the time for kron lambda1....lambdak: Corollary 2
        tic;
        lambda_kron = eig(A{k});
        for p = (k-1):-1:1
            lambda_kron = sort(kron(lambda_kron, eig(A{p})));
        end
        time_eig_kron = toc;
        eig_kron = lambda_kron;
        time_eig_kron_array(ni, ki) = time_eig_kron;
        
        % Compute the difference to check whether the two sets of
        % eigenvalues are equal
        diff = max(abs((eig_An) - (eig_kron)));
        
        % Display the times
        fprintf('n = %d, k = %d\n', n, k);
        fprintf('Time for eig(An): %.6f seconds\n', time_eig_An);
        fprintf('Time for kron-eig computation: %.6f seconds\n', time_eig_kron);
        
        % Check if the difference is larger than machine precision
        if diff > tolerance
            warning('Difference between eig(An) and kron-eig computation is larger than machine precision: %.12f', diff);
        else
            fprintf('Difference is within machine precision: %.12f\n\n', diff);
        end
    end
end

%% plot the comparison
figure;
for ki = 1:length(k_values)
    subplot(1, length(k_values), ki);
    plot(n_values, time_eig_An_array(:, ki), '-o', 'DisplayName', 'eig(A)');
    hold on;
    plot(n_values, time_eig_kron_array(:, ki), '-x', 'DisplayName', 'kron eig(Ap)');
    hold off;
    title(['Comparison for k = ', num2str(k_values(ki))]);
    xlabel('n (Matrix Size)');
    ylabel('Computation Time (s)');
    legend('Location', 'northwest');
    grid on;
end
%%
figure(2)
plot(k_values(1:8), log(time_eig_An_array(1:8)),  '-o', 'DisplayName', 'eig(A)', 'LineWidth', 1)
hold on;
plot(k_values(1:8), log(time_eig_kron_array(1:8)),'-x', 'DisplayName', 'kron eig(A_p)', 'LineWidth', 1)
hold off;

title('Comparison for matrix size n = 3', 'FontSize', 14)
xlabel('k', 'FontWeight', 'bold')
ylabel('Time (s)', 'FontWeight', 'bold')

legend('Location', 'northwest'); 
%% fix k, increase n
figure(2)
plot(n_values, log(time_eig_An_array),  '-o', 'DisplayName', 'eig(A)', 'LineWidth', 1)
hold on;
plot(n_values, log(time_eig_kron_array),'-x', 'DisplayName', 'kron eig(A_p)', 'LineWidth', 1)
hold off;

title('Comparison for k = 3', 'FontSize', 14)
xlabel('Matrix Size n', 'FontWeight', 'bold')
ylabel('Time (s)', 'FontWeight', 'bold')

legend('Location', 'northwest'); 

%% Storage
time1 = time_eig_An_array; 
time2 = time_eig_kron_array; % kron(sort(kron(A1, A2)),A3)...
