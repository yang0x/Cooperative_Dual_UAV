% Model: cooperative dual-UAV system
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


d_sk = zeros(K,N);
d_se = zeros(1,N);
d_je = zeros(1,N);

h_sk = zeros(K, N);
h_se = zeros(E, N);
h_je = zeros(E, N);
R_sec1 = zeros(K, N);

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
res_opt = zeros(1,MAX_Iter);
res_real = zeros(1,MAX_Iter);

epsilon = 1e-4;

for iter = 1:MAX_Iter
    
    % 1. Timeslot Scheduling Optimization - A
    for k = 1 : K
        d_sk(k,:) = norms(q_s(:,:) - Loc_K(:, k)).^2 + H_s^2;
        h_sk(k, :) = beta_0./d_sk(k,:);
    end

    d_se(1,:) = norms(q_s(:,:) - Loc_E(:,1)).^2 + H_s^2;
    h_se(1, :) = beta_0 ./ d_se(1,:);

    d_je(1,:) = norms(q_j(:,:) - Loc_E(:,1)).^2 + H_j^2;
    h_je(1, :) = beta_0 ./ d_je(1,:);

    R_sec1 = log(1 + (ones(K,1) * ps) .* h_sk / sigma2) / log(2) ...
            - log(1 + (ps .* h_se) ./ (pj .* h_je + sigma2)) / log(2);

    cvx_begin quiet
    
        variable a_k1(K, N)

        ave_R1 = sum(a_k1 .* R_sec1, 2) / N;

        maximize min(ave_R1)
        
        subject to

            ave_R1 >= R_min;

            0 <= a_k1 <= 1;

            sum(a_k1, 1) <= ones(1, N);                     

    cvx_end
    
    [~, max_row] = max(a_k1, [], 1);
    a_k1_tmp = zeros(size(a_k1));
    a_k1_tmp(sub2ind(size(a_k1), max_row, 1:size(a_k1,2))) = 1;
    a_k = a_k1_tmp;
    
    % 2. Transmit Power Optimization - P
    snr_sk =  h_sk / sigma2;

    cvx_begin quiet
    
        variable ps2(1,N)
        variable pj2(1,N)

        R_sec2 = log(1 + ((ones(K,1) * ps2) .* snr_sk)) / log(2) ...
                  + ones(K,1) * (log( pj2 .* h_je + sigma2) / log(2)) ...
                  - ones(K,1) * (log( pj .* h_je + ps .* h_se + sigma2 ) / log(2)) ...
                  - ones(K,1) * (( h_je .* (pj2 -pj) + h_se .* (ps2 -ps) ) ./ ( log(2) * ( pj .* h_je + ps .* h_se + sigma2 ) ));

        ave_R2 = sum(a_k .* R_sec2, 2) / N;

        maximize min(ave_R2)

        subject to

            ave_R2 >= R_min;

            0 <= ps2 <= ps_max;

            0 <= pj2 <= pj_max;

            sum(ps2 + pj2, 2) / N <= 2 * P_t_ave;

    cvx_end

    ps = ps2;
    pj = pj2;

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
    
        variable q_s3(2,N)
        variable d_s3_1(K,N)
        variable d_s3_2(1,N)
        variable g_s3(1,N-1)   

        R_e3 = ( rel_entr_quad(d_s3_2 ./ gamma_s_2, d_s3_2 ./ gamma_s_2 + 1 ) + rel_entr_quad( d_s3_2 ./ gamma_s_2 + 1 , d_s3_2 ./ gamma_s_2 ) ) / log(2);

        R_sec3 = log(1 + gamma_s_1 ./ d_sk ) / log(2) ...
                 - ( gamma_s_1 .* (d_s3_1 - d_sk) ) ./ ( log(2) * ( (d_sk).^2 + gamma_s_1 .* d_sk ) ) ...
                 - ones(K,1) * R_e3;

        ave_R3 = sum(a_k .* R_sec3, 2) / N;

        v_s3 = norms( q_s3(:, 2:end)- q_s3(:,1:end-1) , 2 ) / delta_t;
    
        p_fly_s3 = param_P0 * ( 1 + 3 * pow_pos(v_s3/param_Utip ,2) ) ... 
                   + 1/2 * param_d0 * param_rho * param_A * param_s * pow_pos(v_s3,3) ...
                   + param_Pi * g_s3;

        maximize min(ave_R3)

        subject to

            for k = 1 : K
                d_s3_1(k,:) >= pow_pos( norms(q_s3 - (Loc_K(:, k) * ones(1,N)) ),2) + H_s^2;
            end

            d_s3_2 <= d_se + sum( 2 * ( q_s - Loc_E * ones(1,N) ) .* (q_s3 - q_s) , 1);

            ave_R3 >= R_min;

            q_s3(1,1) == q_s(1, 1);
            q_s3(2,1) == q_s(2, 1);
            q_s3(1,N) == q_s(1, end);
            q_s3(2,N) == q_s(2, end);

            v_s3 <= V_max;           

            sum(p_fly_s3+ p_fly_j,2) <= P_p_ave * 2 * (N-1);

            g_s3 >= 0;

            pow_pos(inv_pos(g_s3),2) <= (g_s).^2 + 2 * g_s .* ( g_s3 - g_s ) ...
                                         - (v_s/param_v0).^2 + ( 2 / (delta_t * param_v0)^2 ) * sum( (q_s(:,2:end) - q_s(:,1:end-1)) .* (q_s3(:,2:end) - q_s3(:,1:end-1)) ,1); 

    cvx_end   
       
    q_s(:, :) = q_s3(:, :); 

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
    
        variable q_j4(2,N) 
        variable d_j4_1(1,N)
        variable d_j4_2(1,N)
        variable g_j4(1,N-1)     

        R_sec4 = log(1 + ((ones(K,1) * ps) .* h_sk) / sigma2) / log(2) ...
                     - ones(K,1) * ( (rel_entr_quad( d_j4_1 .* (ps .* (h_se/beta_0) + (sigma2/beta_0) ) ./ pj , ( d_j4_1 .* (ps .* (h_se/beta_0) + (sigma2/beta_0) ) ./ pj ) + 1 ) ...
                                      + rel_entr_quad( ( d_j4_1 .* (ps .* (h_se/beta_0) + (sigma2/beta_0) ) ./ pj ) + 1 , d_j4_1 .* (ps .* (h_se/beta_0) + (sigma2/beta_0) ) ./ pj ) ) / log(2) ) ...
                     - ones(K,1) * (log( ps .* h_se + sigma2 ) / log(2)) ...
                     + ones(K,1) * (log( pj .* h_je + sigma2) / log(2)) ...
                     - ones(K,1) * ( ( pj .* (d_j4_2 - d_je) ) ./ (log(2) * ( ( pj .* d_je ) + ( (sigma2/beta_0) .* ((d_je).^2) ) ) ) );

        ave_R4 = sum(a_k .* R_sec4, 2) / N;

        v_j4 = norms( q_j4(:, 2:end)- q_j4(:,1:end-1) , 2 ) / delta_t;

        p_fly_j4 = param_P0 * ( 1 + 3 * pow_pos(v_j4/param_Utip ,2) )... 
                   + 1/2 * param_d0 * param_rho * param_A * param_s * pow_pos(v_j4,3) ...
                   + param_Pi * g_j4;

        maximize min(ave_R4)
    
        subject to        

            d_j4_1 <= d_je + sum( 2 * ( q_j - Loc_E * ones(1,N) ) .* (q_j4 - q_j) , 1);

            d_j4_2 >= pow_pos( norms(q_j4(:,:) - (Loc_E(:,1) * ones(1,N)) ) , 2) + H_j^2;                   

            ave_R4 >= R_min;

            q_j4(1,1) == q_j(1, 1);
            q_j4(2,1) == q_j(2, 1);
            q_j4(1, N) == q_j(1, end);
            q_j4(2,N) == q_j(2, end);

            v_j4 <= V_max;

            sum(p_fly_s+ p_fly_j4,2) <= P_p_ave * 2 * (N-1);

            g_j4 >= 0;

            pow_pos(inv_pos(g_j4),2) <= (g_j).^2 + 2 * g_j .* ( g_j4 - g_j ) ...
                                         - (v_j/param_v0).^2 + ( 2 / (delta_t * param_v0)^2 ) * sum( (q_j(:,2:end) - q_j(:,1:end-1)) .* (q_j4(:,2:end) - q_j4(:,1:end-1)) , 1);

    
    cvx_end   
       
    q_j(:, :) = q_j4(:, :); 

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
        break;
    end

end


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







