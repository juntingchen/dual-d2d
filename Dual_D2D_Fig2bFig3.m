% This script generates simulation results depicted in Fig. 3 and case B in 
% Fig. 2b in the paper
%
% J. Chen, H. Yin, L. Cottatellucci, and D. Gesbert, "Dual-regularized 
% Feedback and Precoding for D2D Assisted MIMO Systems", IEEE Trans. 
% Wireless Commun., 2017.
%
% Author: Junting CHEN (EURECOM)
% Email: chenju@eurecom.fr
% Date: 14 April 2016
%
% Modification
% - data: 10 May 2016
% - Changes: 
%       - non-equal power allocaiton among users, change sum_rate -->
%         sum_rate_neql
%       - Heterogeneous CSI noise
% 
% -------------------------------------------------------------------------
clear
close all
addpath Sub

Nt = 20;    % number of TX antennas
K = 3;      % number of users
Bf = 10;     % Number of bits to feedback to BS

D2D_quality = [1     1     1
               1     1     0
               1     0     1]; % <- Scaling factor of the D2D CSI quality
                               % 1 - some quality, 0 - no CSI at all
SNR_vec = [20];      % SNR
Nsnr = length(SNR_vec) - 1;

Bd2d = [0:15:150 0];
Squared_alpha = 1 - 2 .^ (- Bd2d / Nt);
% Bd2d(length(Bd2d) - 1) = 180;
% Squared_alpha = [0:0.1:1, 0.84];
sd_vec = sqrt(1 - Squared_alpha.^2);
% sd_vec = [0:0.1:1 0.4];     % \sigma_d
Nsd = length(sd_vec) - 1;

Nmc = 100000;                 % Number of Monte-Carlo simulations for channel 
                            % realizations on each case (SNR, sigma_d)

str_sch = {
    'Precoder FB, perft CSI'
    'Proposed'
    'Precoder feedback'
    'CSI feedback'
    };
Nsch = 4;
schemes_to_show = [2 3 4];

% Universal codebook design
Mbf = 2 ^ Bf;
Cb = zeros(Mbf, Nt);
for i = 1:Mbf
    u = randn(Nt, 1) + 1i * randn(Nt, 1);
    u = u / norm(u);
    Cb(i, :) = u.';
end
Sbf = 2 ^ (- Bf / Nt);  % Regularization parameter from [Dabbagh & Love 2008]
Sbf1 = 2 ^ (- Bf / (Nt - 1));

Ncase = Nsnr + Nsd;
Rsum = zeros(Nsch, Ncase);
Rates = zeros(Nsch, Ncase, K);

for i_case = 1:Ncase
    % Case arrangement: (i) go over different SNRs, (ii) go over different Sigma
    if i_case <= Nsnr
        sd = sd_vec(end);
        SNR_dB = SNR_vec(i_case);
    else
        sd = sd_vec(i_case - Nsnr);
        SNR_dB = SNR_vec(end);
    end
    P = 10 ^ (SNR_dB / 10);
    alpha = K / P;
    
    Sigma = ones(K) - (1 - sd) * D2D_quality;
    Sigma = Sigma .* (ones(K) - eye(K));        % Diagonal elements being zeros
    
    
    rsum = zeros(Nsch, Nmc);
    rates = zeros(K, Nsch, Nmc);
    
    % Parameters for robust-D2D
    Alphas2 = 1 - Sigma .^ 2;
    Varrho = sum(ones(K) - Alphas2, 2) + K / P * ones(K, 1);
    Omega = diag(Varrho .^ 2);
    Phi = sum(Varrho .^2 ) + (sum(Alphas2(:) .^ 2) - K) * (Nt + 1) ...
            + 2 * sum((sum(Alphas2, 2) - ones(K, 1)) .* Varrho);
    Phi = Phi / Nt;
    %
    
    for i = 1:Nmc
        Sigma = ones(K) - (1 - sd) * D2D_quality;
        Sigma = Sigma .* (ones(K) - eye(K));        % Diagonal elements being zeros
    
        H = (randn(Nt, K) + 1i * randn(Nt, K)) / sqrt(2);
        % Hd = sqrt(1 - sd^2) * H + sd * (randn(Nt, K) + 1i * randn(Nt, K)) / sqrt(2);    % Noisy D2D
        Hd = cell(K, 1);
        for k = 1:K
            Hd{k} = H * diag(sqrt(1 - Sigma(k, :) .^2 )) + ...
                    (randn(Nt, K) + 1i * randn(Nt, K))  / sqrt(2) * diag(Sigma(k, :));
            % Hd{k}(:, k) = H(:, k);
        end
        Gain = diag(H' * H);
        
%         % D2D Perfect
%         Q = abs(conj(Cb) * H) .^ 2;
%         Qsum = sum(Q, 2);
%         Qsum_alpha = Qsum + ones(Mbf, 1) * alpha;
%         W1 = zeros(Nt, K);
%         for k = 1:K
%             Qk = Q(:, k) ./ (Qsum_alpha - Q(:, k));
%             [~, imk] = max(Qk);
%             W1(:, k) = Cb(imk, :).';
%         end
%         [rsum(1, i), rates(:, 1, i)] = sum_rate_neql(H, W1, P);
        
        % D2D robust (Proposed dual regularized scheme)
        H2 = zeros(Nt, K);
        for k = 1:K
            Hk = Hd{k};
            Q = abs(conj(Cb) * (Hk * diag(ones(1, K) - Sigma(k, :) .^ 2))) .^ 2;
            Qsum = sum(Q, 2);
            Qsum_alpha = (Qsum - Q(:, k)) ...
                          + ones(Mbf, 1) * (alpha + sum(Sigma(k, :).^2));
            Qk = Q(:, k) ./ Qsum_alpha;
            
            [~, imk] = max(Qk);
            H2(:, k) = Cb(imk, :).';
        end
        Sigma = ones(K) - (1 - sd) * ones(K) * min(D2D_quality(:));
        Sigma = Sigma .* (ones(K) - eye(K)); 
        W4 = ((1 - Sbf1) * H2 * Omega * H2' + (Sbf1 * Phi + K / P) * eye(Nt)) ...
                \ H2;
        [rsum(2, i), rates(:, 2, i)] = sum_rate(H, W4, P);
        
        % D2D non-robust
        W5 = zeros(Nt, K);
        for k = 1:K
            Hk = Hd{k};
            Q = abs(conj(Cb) * Hk) .^ 2;
            Qsum = sum(Q, 2);
            Qsum_alpha = Qsum + ones(Mbf, 1) * (alpha);
            Qk = Q(:, k) ./ (Qsum_alpha - Q(:, k));
            [~, imk] = max(Qk);
            W5(:, k) = Cb(imk, :).';
        end
        [rsum(3, i), rates(:, 3, i)] = sum_rate_neql(H, W5, P);
        
        % BS precoding based on channel feedback
        Qb = abs(conj(Cb) * H * diag(1 ./ sqrt(Gain))) .^ 2;
        [~, I] = max(Qb);
        Hfb0 = Cb(I, :).';   % Feedback only the direction
        Hfb = Hfb0;

        % BS Robust MMSE
        W6 = Hfb / (Hfb' * Hfb + (alpha + K * Sbf) * eye(K));
        [rsum(4, i), rates(:, 4, i)] = sum_rate_neql(H, W6, P);
    end
    
    Rsum(:, i_case) = mean(rsum, 2);
    for k = 1:K
        for i_sch = 1:Nsch
            Rates(i_sch, i_case, k) = mean(rates(k, i_sch, :));
        end
    end
    
end

%% Figures
my_markers = {'d', 's', '^', 'o', '<', '>', 'v', '+', 'p', '.', '*'}.';
my_colors4 = [0 0 1 %'b'
              0    0.4980         0 %'g'
              1 0 0 %'r'
              0    0.7490    0.7490];%'c'
my_colors3 = my_colors4(2:end, :);
Nschemes_to_show = length(schemes_to_show);

if Nsnr > 1
    hf = figure(3);
    set(gca, 'ColorOrder', my_colors4, 'NextPlot', 'replacechildren');
    p_handle = plot(SNR_vec(1:end - 1), Rsum(schemes_to_show, 1:Nsnr), 'LineWidth', 2, 'MarkerSize', 9);
    set(p_handle, {'marker'}, my_markers(1:Nschemes_to_show));
    set(gca, 'FontSize', 14);
    legend(str_sch{schemes_to_show}, 'location', 'Southeast');
    ylabel('Sum rate (bps/Hz)');
    xlabel('Bits Bd for CSI exchange');

end

if Nsd > 1
    hf = figure(4);
    set(gca, 'ColorOrder', my_colors3, 'NextPlot', 'replacechildren');
    % p_handle = plot(sd_vec(1:end - 1), Rsum(schemes_to_show, Nsnr + 1:end), 'LineWidth', 2, 'MarkerSize', 9);
    p_handle = plot(Bd2d(1: end - 1), Rsum(schemes_to_show, Nsnr + 1:end), 'LineWidth', 2, 'MarkerSize', 9);
    set(p_handle, {'marker'}, my_markers(1:Nschemes_to_show));
    set(gca, 'FontSize', 14);
    legend(str_sch{schemes_to_show}, 'location', 'Southeast');
    ylabel('Sum rate (bps/Hz)');
    xlabel('Bits Bd for CSI exchange');
    % title(sprintf('Nt = %d, K = %d, SNR = %d dB', Nt, K, round(SNR_vec(end))));
    ylim([0 15]);
    

    
    % figure(5),
    for k = 1:K
        % subplot(1, K, k),
        hf = figure(4 + k);
        set(gca, 'ColorOrder', my_colors3, 'NextPlot', 'replacechildren');
        p_handle = plot(Bd2d(1:end - 1), Rates(schemes_to_show, Nsnr + 1:end, k), 'LineWidth', 2, 'MarkerSize', 9);
        set(p_handle, {'marker'}, my_markers(1:Nschemes_to_show));
        set(gca, 'FontSize', 14);
        legend(str_sch{schemes_to_show}, 'location', 'Southwest');
        ylabel('Sum rate (bps/Hz)');
        xlabel('Bits Bd for CSI exchange');
        % title(sprintf('User %d', k));
        ylim([0 5]);
        set(gca, 'YTick', 0:5);
        

    end
end
