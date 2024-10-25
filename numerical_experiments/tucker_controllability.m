A1 = rand(4,4);
A2 = rand(4,4);
A3 = rand(4,4);
A4 = rand(5,5);
A5 = rand(5,5);
%An = kron(kron(kron(kron(A5,A4),A3),A2),A1);
An = kron(kron(A3,A2), A1);
%sort(kron(eig(A2), eig(A1))) - sort(eig(An))
%max(abs(sort(kron(eig(A2), eig(A1))) - sort(eig(An))))
sort(eig(An));
%sort(kron(kron(eig(A2),eig(A3)),eig(A1))) - sort(eig(An));
sort(kron(sort(kron(eig(A1),eig(A2))),eig(A3))) - sort(eig(An));
sort(kron(kron(eig(A1),eig(A2)),eig(A3))) - sort(eig(An));
%sort(kron(sort(kron(sort(kron(sort(kron(eig(A1),eig(A2))),eig(A3))),eig(A4))),eig(A5))) - sort(eig(An))
%sort(kron(kron(kron(kron(eig(A1),eig(A2)),eig(A3)),eig(A4)),eig(A5))) - sort(eig(An))
%% 2. ctrb
ctrb(A1, rand(2,3))

%% 1. check time of computing eig(A) vs (\sum_p eig(Ap)) + lambda1 kron lambda2 kron ... kron lambdak
%--------- prop 1 vs corollary 2
% n: each Ap is of size n*n
% k: there are k Ap's

n_values = [2,3,4,5,6,7,8,9,10,11,12];
k_values = [3];
%n_values = 2;
%k_values = 2:12;
tolerance = 1e-14;
ntrials = 100;
% n*k matrix for time: each row shows time as k increases; each col shows time as n increases
time_eig_An_array = zeros(length(n_values), length(k_values)); 
time_eig_kron_array = zeros(length(n_values), length(k_values));
for ni = 1:length(n_values)
    n = n_values(ni);
    for ki = 1:length(k_values)
        k = k_values(ki);

        % Initialize time accumulators for each (n, k)
        time_eig_An_total = 0;
        time_eig_kron_total = 0;

        for trial = 1:ntrials
            disp(trial);
            
            A = cell(1, k);
            for p = 1:k
                [V, ~] = eig(rand(n));  
                D = diag(rand(1, n));   
                A{p} = V * D / V;
            end

            % eigenvalue for An = Ak kron A_{k-1} kron... kron A_1
            An = A{k};
            for p = (k-1):-1:1  
                An = kron(An, A{p});
            end
            
            % Measure the time for eig(An): Proposition 1
            tic;
            eig_An = (eig(An));
            time_eig_An = toc;
            time_eig_An_total = time_eig_An_total + time_eig_An;
            
            % Measure the time for kron lambda1....lambdak: Corollary 2
            tic;
            lambda_kron = eig(A{k});
            for p = (k-1):-1:1
                lambda_kron = (kron(lambda_kron, eig(A{p})));
            end
            time_eig_kron = toc;
            time_eig_kron_total = time_eig_kron_total + time_eig_kron;
        end

        % Average the times for the current (n, k)
        time_eig_An_array(ni, ki) = time_eig_An_total / ntrials;
        time_eig_kron_array(ni, ki) = time_eig_kron_total / ntrials;

        % Compute the difference to check whether the two sets of eigenvalues are equal
        diff = max(abs(sort(eig_An) - sort(lambda_kron)));

        % Display the times
        fprintf('n = %d, k = %d\n', n, k);
        fprintf('Time for eig(An): %.6f seconds\n', time_eig_An_array(ni, ki));
        fprintf('Time for kron-eig computation: %.6f seconds\n', time_eig_kron_array(ni, ki));
        
        % Check if the difference is larger than machine precision
        if diff > tolerance
            warning('Difference between eig(An) and kron-eig computation is larger than machine precision: %.12f', diff);
        else
            fprintf('Difference is within machine precision: %.12f\n\n', diff);
        end
    end
end




%% fix n, increase k
figure(2)
set(gcf, 'Position', [100, 100, 600, 300]); 
plot(k_values, log(time_eig_An_array),  '-o', 'DisplayName', 'eig(A)', 'LineWidth', 2, 'Color', [1, 0, 0])
hold on;
plot(k_values, log(time_eig_kron_array),'-x', 'DisplayName', 'kron eig(A_p)', 'LineWidth', 2, 'Color', [0, 0, 1])
hold off;

title('Comparison for matrix size n = 2', 'FontSize', 14)
xlabel('k', 'FontWeight', 'bold', 'FontSize', 20)
ylabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 20)
set(gca, 'FontSize', 13)

legend('Location', 'northwest','FontSize', 14); 
%% fix k, increase n
figure(2)
set(gcf, 'Position', [100, 100, 600, 300]); 
plot(n_values, log(time_eig_An_array),  '-o', 'DisplayName', 'eig(A)', 'LineWidth', 2, 'Color', [1, 0, 0])
hold on;
plot(n_values, log(time_eig_kron_array),'-x', 'DisplayName', 'kron eig(A_p)', 'LineWidth', 2, 'Color', [0, 0, 1])
hold off;

title('Comparison for k = 3', 'FontSize', 14)
xlabel('Matrix Size n', 'FontWeight', 'bold', 'FontSize', 18)
ylabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 18)
set(gca, 'FontSize', 13)

legend('Location', 'northwest','FontSize', 14); 

%% Storage
time1 = time_eig_An_array; 
time2 = time_eig_kron_array; % kron(sort(kron(A1, A2)),A3)...

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