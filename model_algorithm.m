% Model: cooperative dual-UAV model
% Algorithm: based on the BCD and SCA
% Optimize variables and order: (a - p - traj_s - traj_j)

clc;
close all;

cvx_solver Mosek;

%% Params
T = 40;
P_p_ave = 500;
P_t_ave = 0.1;
R_min = 2;
delta_t = 0.5;
N = T / delta_t;
H_s = 100;
H_j = 80;
V_max = 50;
MAX_X = 400;
MAX_Y = 400;
beta_0 = 1e-6;
sigma2 = 1e-14;

K = 3;
E = 1;
Loc_K = [100 240 350;60 320 100];
Loc_E = [220 160]';

% UAV prop param
param_Utip = 120;
param_P0 = 79.8573;
param_Pi = 88.6269;
param_d0 = 0.6;
param_rho = 1.225;
param_A = 0.503;
param_s = 0.05;
param_v0 = 4.03;

% initial solution
q_s = zeros(2,N);
q_s(:,1) = [0,MAX_Y/2]';
q_s(:,end) = [MAX_X,MAX_Y/2]';
for n=2:N-1
    length = MAX_X * ((n-1)/(N-1));
    q_s(1,n) = q_s(1,1) + length;
    q_s(2,n) = q_s(2,1);
end 
q_j = zeros(2,N);
q_j(1,:) = q_s(1,:);
q_j(2,:) = q_s(2,:) + 150;

init_s = q_s;
init_j = q_j;

ps = P_t_ave * ones(1,N);
ps_max = 4 * P_t_ave;

pj = P_t_ave * ones(1,N);
pj_max = 4 * P_t_ave;

% BCD temporary variable
d_sk = zeros(K,N);
d_se = zeros(1,N);
d_je = zeros(1,N);

h_sk = zeros(K, N);
h_se = zeros(E, N);
h_je = zeros(E, N);
R_sec = zeros(K, N);

snr_sk = zeros(K,N);

gamma_s_1 = zeros(K,N);
gamma_s_2 = zeros(1,N);

v_s = zeros(1,N-1);
g_j = zeros(1,N-1);
p_fly_j = zeros(1,N-1);

v_j = zeros(1,N-1);
g_s = zeros(1,N-1);
p_fly_s = zeros(1,N-1);


h_sk_real = zeros(K, N);
h_se_real = zeros(E, N);
h_je_real = zeros(E, N);
R_sec_real = zeros(K, N);
ave_R_real = zeros(K, 1);
Max_min_ASR = 0;

MAX_Iter = 30;
epsilon = 1e-4;
res_opt = zeros(1,MAX_Iter);
res_real = zeros(1,MAX_Iter);
total_time = 0;


for iter = 1:MAX_Iter
    
    tic;
    
    % 1. Timeslot Scheduling Optimization - A
    for k = 1 : K
        d_sk(k,:) = norms(q_s(:,:) - Loc_K(:, k)).^2 + H_s^2;
        h_sk(k, :) = beta_0./d_sk(k,:);
    end

    d_se(1,:) = norms(q_s(:,:) - Loc_E(:,1)).^2 + H_s^2;
    h_se(1, :) = beta_0 ./ d_se(1,:);

    d_je(1,:) = norms(q_j(:,:) - Loc_E(:,1)).^2 + H_j^2;
    h_je(1, :) = beta_0 ./ d_je(1,:);

    R_sec = log(1 + (ones(K,1) * ps) .* h_sk / sigma2) / log(2) ...
            - log(1 + (ps .* h_se) ./ (pj .* h_je + sigma2)) / log(2);

    cvx_begin quiet
    
        variable a_k_var(K, N)

        ave_R1_var = sum(a_k_var .* R_sec, 2) / N;

        maximize min(ave_R1_var)
        
        subject to

            ave_R1_var >= R_min;

            0 <= a_k_var <= 1;

            sum(a_k_var, 1) <= ones(1, N);                     

    cvx_end

    if isequal(cvx_status, 'Solved')    
        [~, max_row] = max(a_k_var, [], 1);
        a_k_var_tmp = zeros(size(a_k_var));
        a_k_var_tmp(sub2ind(size(a_k_var), max_row, 1:size(a_k_var,2))) = 1;
        a_k = a_k_var_tmp;
    else 
        break;
    end    
    
    % 2. Transmit Power Optimization - P
    snr_sk =  h_sk / sigma2;

    cvx_begin quiet
    
        variable p_s_var(1,N)
        variable p_j_var(1,N)

        R_sec2_var = log(1 + ((ones(K,1) * p_s_var) .* snr_sk)) / log(2) ...
                  + ones(K,1) * (log( p_j_var .* h_je + sigma2) / log(2)) ...
                  - ones(K,1) * (log( pj .* h_je + ps .* h_se + sigma2 ) / log(2)) ...
                  - ones(K,1) * (( h_je .* (p_j_var -pj) + h_se .* (p_s_var -ps) ) ./ ( log(2) * ( pj .* h_je + ps .* h_se + sigma2 ) ));

        ave_R2_var = sum(a_k .* R_sec2_var, 2) / N;

        maximize min(ave_R2_var)

        subject to

            ave_R2_var >= R_min;

            0 <= p_s_var <= ps_max;

            0 <= p_j_var <= pj_max;

            sum(p_s_var + p_j_var, 2) / N <= 2 * P_t_ave;

    cvx_end

    if isequal(cvx_status, 'Solved')    
        ps = p_s_var;
        pj = p_j_var;
    else 
        break;
    end 

    % 3. US Trajectory Optimization - Q_s
    gamma_s_1 = ones(K,1) * (ps * (beta_0 / sigma2));

    gamma_s_2 = ps .* (beta_0 ./ (pj .* h_je + sigma2));

    v_s = norms( q_s(:,2:end) - q_s(:,1:end-1) , 2 ) / delta_t;

    g_s = ( ( 1 + 1/4 * (v_s/param_v0).^4 ).^(1/2) - 1/2 * (v_s/param_v0).^2 ).^(1/2);

    v_j = norms( q_j(:,2:end) - q_j(:,1:end-1) , 2 ) / delta_t;

    p_fly_j = param_P0 * ( 1 + 3 * (v_j/param_Utip).^2 )... 
               + 1/2 * param_d0 * param_rho * param_A * param_s * (v_j).^3 ...
               + param_Pi * ( ( ( 1 + 1/4 * (v_j/param_v0).^4 ).^(1/2) - 1/2 * (v_j/param_v0).^2 ).^(1/2) );
            

    cvx_begin quiet
    
        variable q_s_var(2,N)
        variable d_s_1_var(K,N)
        variable d_s_2_var(1,N)
        variable g_s_var(1,N-1)   

        R_e3_var = ( rel_entr_quad(d_s_2_var ./ gamma_s_2, d_s_2_var ./ gamma_s_2 + 1 ) + rel_entr_quad( d_s_2_var ./ gamma_s_2 + 1 , d_s_2_var ./ gamma_s_2 ) ) / log(2);

        R_sec3_var = log(1 + gamma_s_1 ./ d_sk ) / log(2) ...
                 - ( gamma_s_1 .* (d_s_1_var - d_sk) ) ./ ( log(2) * ( (d_sk).^2 + gamma_s_1 .* d_sk ) ) ...
                 - ones(K,1) * R_e3_var;

        ave_R3_var = sum(a_k .* R_sec3_var, 2) / N;

        v_s3_var = norms( q_s_var(:, 2:end)- q_s_var(:,1:end-1) , 2 ) / delta_t;
    
        p_fly_s3_var = param_P0 * ( 1 + 3 * pow_pos(v_s3_var/param_Utip ,2) ) ... 
                   + 1/2 * param_d0 * param_rho * param_A * param_s * pow_pos(v_s3_var,3) ...
                   + param_Pi * g_s_var;

        maximize min(ave_R3_var)

        subject to

            for k = 1 : K
                d_s_1_var(k,:) >= pow_pos( norms(q_s_var - (Loc_K(:, k) * ones(1,N)) ),2) + H_s^2;
            end

            d_s_2_var <= d_se + sum( 2 * ( q_s - Loc_E * ones(1,N) ) .* (q_s_var - q_s) , 1);

            ave_R3_var >= R_min;

            q_s_var(1,1) == q_s(1, 1);
            q_s_var(2,1) == q_s(2, 1);
            q_s_var(1,N) == q_s(1, end);
            q_s_var(2,N) == q_s(2, end);

            v_s3_var <= V_max;           

            sum(p_fly_s3_var+ p_fly_j,2) <= P_p_ave * 2 * (N-1);

            g_s_var >= 0;

            pow_pos(inv_pos(g_s_var),2) <= (g_s).^2 + 2 * g_s .* ( g_s_var - g_s ) ...
                                         - (v_s/param_v0).^2 + ( 2 / (delta_t * param_v0)^2 ) * sum( (q_s(:,2:end) - q_s(:,1:end-1)) .* (q_s_var(:,2:end) - q_s_var(:,1:end-1)) ,1); 

    cvx_end   
       
    if isequal(cvx_status, 'Solved')    
        q_s = q_s_var;
    else 
        break;
    end 

    % 4. UJ Trajectory Optimization - Q_j
    for k = 1 : K
        h_sk(k, :) = beta_0./(norms(q_s - Loc_K(:, k)).^2 + H_s^2);
    end

    h_se(1, :) = beta_0 ./ (norms(q_s - Loc_E).^2 + H_s^2);

    v_j = norms( q_j(:,2:end) - q_j(:,1:end-1) , 2 ) / delta_t;

    g_j = ( ( 1 + 1/4 * (v_j/param_v0).^4 ).^(1/2) - 1/2 * (v_j/param_v0).^2 ).^(1/2);

    v_s = norms( q_s(:,2:end) - q_s(:,1:end-1) , 2 ) / delta_t;

    p_fly_s = param_P0 * ( 1 + 3 * (v_s/param_Utip).^2 )... 
               + 1/2 * param_d0 * param_rho * param_A * param_s * (v_s).^3 ...
               + param_Pi * ( ( ( 1 + 1/4 * (v_s/param_v0).^4 ).^(1/2) - 1/2 * (v_s/param_v0).^2 ).^(1/2) );
  
    cvx_begin quiet
    
        variable q_j_var(2,N) 
        variable d_j_1_var(1,N)
        variable d_j_2_var(1,N)
        variable g_j_var(1,N-1)     

        R_sec4_var = log(1 + ((ones(K,1) * ps) .* h_sk) / sigma2) / log(2) ...
                     - ones(K,1) * ( (rel_entr_quad( d_j_1_var .* (ps .* (h_se/beta_0) + (sigma2/beta_0) ) ./ pj , ( d_j_1_var .* (ps .* (h_se/beta_0) + (sigma2/beta_0) ) ./ pj ) + 1 ) ...
                                      + rel_entr_quad( ( d_j_1_var .* (ps .* (h_se/beta_0) + (sigma2/beta_0) ) ./ pj ) + 1 , d_j_1_var .* (ps .* (h_se/beta_0) + (sigma2/beta_0) ) ./ pj ) ) / log(2) ) ...
                     - ones(K,1) * (log( ps .* h_se + sigma2 ) / log(2)) ...
                     + ones(K,1) * (log( pj .* h_je + sigma2) / log(2)) ...
                     - ones(K,1) * ( ( pj .* (d_j_2_var - d_je) ) ./ (log(2) * ( ( pj .* d_je ) + ( (sigma2/beta_0) .* ((d_je).^2) ) ) ) );

        ave_R4_var = sum(a_k .* R_sec4_var, 2) / N;

        v_j4_var = norms( q_j_var(:, 2:end)- q_j_var(:,1:end-1) , 2 ) / delta_t;

        p_fly_j4_var = param_P0 * ( 1 + 3 * pow_pos(v_j4_var/param_Utip ,2) )... 
                   + 1/2 * param_d0 * param_rho * param_A * param_s * pow_pos(v_j4_var,3) ...
                   + param_Pi * g_j_var;

        maximize min(ave_R4_var)
    
        subject to        

            d_j_1_var <= d_je + sum( 2 * ( q_j - Loc_E * ones(1,N) ) .* (q_j_var - q_j) , 1);

            d_j_2_var >= pow_pos( norms(q_j_var(:,:) - (Loc_E(:,1) * ones(1,N)) ) , 2) + H_j^2;                   

            ave_R4_var >= R_min;

            q_j_var(1,1) == q_j(1, 1);
            q_j_var(2,1) == q_j(2, 1);
            q_j_var(1, N) == q_j(1, end);
            q_j_var(2,N) == q_j(2, end);

            v_j4_var <= V_max;

            sum(p_fly_s+ p_fly_j4_var,2) <= P_p_ave * 2 * (N-1);

            g_j_var >= 0;

            pow_pos(inv_pos(g_j_var),2) <= (g_j).^2 + 2 * g_j .* ( g_j_var - g_j ) ...
                                         - (v_j/param_v0).^2 + ( 2 / (delta_t * param_v0)^2 ) * sum( (q_j(:,2:end) - q_j(:,1:end-1)) .* (q_j_var(:,2:end) - q_j_var(:,1:end-1)) , 1);

    
    cvx_end   

    if isequal(cvx_status, 'Solved')    
        q_j = q_j_var; 
    else 
        break;
    end 

    iter_time = toc;
    fprintf('iteration %d, elapsed time: %.4f seconds\n', iter, iter_time);

    total_time = total_time + iter_time;

    % After one round of optimization is completed, record the results
    res_opt(1,iter) = cvx_optval;

    for k = 1 : K
        h_sk_real(k, :) = beta_0./(norms(q_s - Loc_K(:, k)).^2 + H_s^2);
    end

    h_se_real(1, :) = beta_0 ./ (norms(q_s - Loc_E(:,1)).^2 + H_s^2);

    h_je_real(1, :) = beta_0 ./ (norms(q_j - Loc_E(:,1)).^2 + H_j^2);

    R_sec_real = log(1 + (ones(K,1) * ps) .* h_sk_real / sigma2) / log(2) ...
            - log(1 + (ps .* h_se_real) ./ (pj .* h_je_real + sigma2)) / log(2);

    ave_R_real = sum(a_k .* R_sec_real, 2) / N;

    Max_min_ASR = min(ave_R_real);

    res_real(1,iter) = Max_min_ASR;

    fprintf(' the %d iteration, optimized value = %f, real value = %f \n', iter, cvx_optval, full(Max_min_ASR));

    % Algorithm convergence figure (optimized value)
    figure(1);
    grid on;
    box on;
    plot(0:MAX_Iter,[0,res_opt],'m-o','markersize' ,8 ,'linewidth',1.5);
    xlabel('Number of Iterations');
    ylabel('Max-min ASR (bps/Hz)');
    drawnow;

    % Algorithm convergence figure (real value)
    figure(2);
    grid on;
    box on;
    plot(0:MAX_Iter,[0,res_real],'m-o','markersize' ,8 ,'linewidth',1.5);
    xlabel('Number of Iterations');
    ylabel('Max-min ASR (bps/Hz)');
    drawnow;
    
    if (iter ~= 1) && (abs(res_real(1, iter) - res_real(1, iter-1)) <= epsilon) 
        fprintf('Number of iterations: %d, algorithm converged! \n', iter);
        break;
    end

end

fprintf('Total iterations: %d, total time: %.4f seconds, average time per iteration: %.4f seconds\n', ...
            iter, total_time, total_time/iter);


% cooperative dual-UAV trajectory figure
figure(3);
hold on;
grid on;
box on;
plot(init_s(1,2:end-1), init_s(2,2:end-1), 'k--', 'MarkerSize', 10,'linewidth', 1);
plot(q_s(1,2:end-1), q_s(2,2:end-1), 'm-o', 'MarkerSize', 10,'linewidth', 1.5);
plot(q_j(1,2:end-1), q_j(2,2:end-1), 'b-*', 'MarkerSize', 10,'linewidth', 1.5);
plot(q_s(1,1),q_s(2,1),'ko','MarkerFaceColor','k','MarkerSize', 9);
text(q_s(1,1)+10, q_s(2,1)+20, '$q_s^I$', 'FontSize', 18,'Color', 'k', 'Interpreter', 'latex');
plot(q_s(1,end),q_s(2,end),'ko','MarkerFaceColor','k','MarkerSize', 9);
text(q_s(1,end)-20, q_s(2,end)+20, '$q_s^F$', 'FontSize', 18,'Color', 'k', 'Interpreter', 'latex');
plot(q_s(1,:), q_s(2,:), 'm-', 'linewidth', 1);    
plot(q_j(1,1),q_j(2,1),'ko','MarkerFaceColor','k','MarkerSize', 9);
text(q_j(1,1)+10, q_j(2,1)+20, '$q_j^I$', 'FontSize', 18,'Color', 'k', 'Interpreter', 'latex');
plot(q_j(1,end),q_j(2,end),'ko','MarkerFaceColor','k','MarkerSize', 9);
text(q_j(1,end)-20, q_j(2,end)+20, '$q_j^F$', 'FontSize', 18,'Color', 'k', 'Interpreter', 'latex');
plot(q_j(1,:), q_j(2,:), 'b-', 'linewidth', 1);
plot(Loc_K(1,1), Loc_K(2,1), 'k^', 'MarkerFaceColor','k', 'MarkerSize', 13);
text(Loc_K(1,1)-15, Loc_K(2,1)-22, 'GU 1', 'FontSize', 17,'Color', 'k','Interpreter', 'latex');
plot(Loc_K(1,2), Loc_K(2,2), 'k^', 'MarkerFaceColor','k', 'MarkerSize', 13);
text(Loc_K(1,2)-15, Loc_K(2,2)-22, 'GU 2', 'FontSize', 17,'Color', 'k', 'Interpreter', 'latex');
plot(Loc_K(1,3), Loc_K(2,3), 'k^', 'MarkerFaceColor','k', 'MarkerSize', 13);
text(Loc_K(1,3)-15, Loc_K(2,3)-22, 'GU 3', 'FontSize', 17,'Color', 'k', 'Interpreter', 'latex');
plot(Loc_E(1,1), Loc_E(2,1), 'kp','MarkerFaceColor','k','MarkerSize', 17);
text(Loc_E(1,1)-10, Loc_E(2,1)-22, 'GE', 'FontSize', 17,'Color', 'k', 'Interpreter', 'latex');
plot(init_j(1,2:end-1), init_j(2,2:end-1), 'k--', 'MarkerSize', 10,'linewidth', 1);
xlabel('x (m)');
ylabel('y (m)');
xlim([0, MAX_X]);
ylim([0, MAX_Y]);
xticks(0:50:MAX_X); 
yticks(0:50:MAX_Y);
legend('Initial trajectory','Optimal US trajectory', 'Optimal UJ trajectory', ...
    'Interpreter','latex','FontSize',14,'FontName','Times New Roman',...
'Position', [0.47, 0.13, 0.2, 0.2]);


