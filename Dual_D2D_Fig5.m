% This script generates simulation results depicted in Fig. 5 in the paper
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
% - Data: 10 May 2016
% - Changes: 
%       - non-equal power allocaiton among users, change sum_rate -->
%         sum_rate_neql
%       - Heterogeneous CSI noise
%
% - Date: 17 Sept 2016
% - Changes:
%       - extension to correlated channels
%       - two user case
%
% - Date: May 9, 2017
% - Changes: Modified for TWC 1st round revision: As requested by 
% reviewer 1, we model the power azimuth spectrum as Laplacian distributed 
% and the angular spread to be a log-normal random variable, according to 
% TR 25.996
% -------------------------------------------------------------------------
clear
close all
addpath Sub

warning('off');
testFigID = 2;

Nt = 20;    % number of TX antennas
K = 2;      % number of users
Bf = 10;    % Number of bits to feedback to BS

Distance_ms2bs = 60;       % user-to-user distance in meters

mu_as = 0.810;
eps_as = 0.34;

P0d2d_array = [100 1000];

SNR_vec = [];       % SNR
SNR0 = 20;          % Default SNR
Nsnr = length(SNR_vec);

Dm2m_vec = 2:2:25;
Dm2m0 = 0;
Nsd = length(Dm2m_vec);

Nmc = 500000;               % Number of Monte-Carlo simulations for channel (500000)
                            % realizations on each case (SNR, sigma_d)
                            % (200)
Nmc1 = 200;                            

str_sch = {
    'D2D perfect'
    'Proposed'
    'Precoder feedback'
    'CSI feedback'
    };
Nsch = 4;
schemes_to_show = [2 3 4];

% % Universal codebook design
Mbf = 2 ^ Bf;

Nd2d = length(P0d2d_array);
Ncase0 = Nsnr + Nsd;
Ncase = Ncase0 * Nd2d;
Rsum = zeros(Nsch, Ncase);
% Rates = zeros(Nsch, Ncase, K);
% Rates_cell = cell(Ncase, 1);

% ------------------------------------------------------------------------
parfor i_case = 1:Ncase
    
    iPd2d = floor((i_case - 1) / Ncase0) + 1;
    P0d2d = P0d2d_array(iPd2d);
    
    iDd2d = mod(i_case - 1, Ncase0) + 1;
    dm2m = Dm2m_vec(iDd2d); % distance between the two (or more) users
    SNR_dB = SNR0;
    P = 10 ^ (SNR_dB / 10);
    alpha = K / P;  % power scaling factor 'alpha' (backward compatibility)

    % User Position and Channel Covariance matrices
    Pos_ue = [-dm2m / 2, dm2m / 2];
    Rue = cell(K, 1);
    Rms = cell(K, 1);
    Aod_ue_deg = Pos_ue / Distance_ms2bs * 180 / pi;

    rsum = zeros(Nsch, Nmc);
    rates = zeros(K, Nsch, Nmc);
    for imc = 1:Nmc
        
        if mod(imc, Nmc1) == 1
            % Update the probabilistic model for channel covariance matrix =======
            for k = 1:K
                % Shadowing
                Sf_k = 10 ^ (randn * 8 / 10); % 8 dB shadowing in standard deivation
                
                AS_deg = 10 .^ (eps_as * randn + mu_as);
                Rue{k} = genChannCorr_lapcn(Nt, Aod_ue_deg(k), AS_deg) * Sf_k;
%                 Rue{k} = genChannCorr_lapcn(Nt, Aod_ue_deg(k), 15);
                % Rue{k} = eye(Nt);
                Rms{k} = sqrtm(Rue{k});
            end
            % Find the overlapped subspace
            POWER_TRUNCATION_THRESHOLD = 0.99;
            [U, D] = eig(Rue{1});
            [D1, I] = sort(real(diag(D)), 'descend');
            [~, M1] = power_truncate(D1, POWER_TRUNCATION_THRESHOLD);
            U1 = U(:, I(1:M1));
            D1 = D1(1:M1);
            %
            [U, D] = eig(Rue{2});
            [D2, I] = sort(real(diag(D)), 'descend');
            [~, M2] = power_truncate(D2, POWER_TRUNCATION_THRESHOLD);
            U2 = U(:, I(1:M2));
            D2 = D2(1:M2); 
            %
            R12 = U2 * U2' * Rue{1} * (U2 * U2');
            [U, D] = eig(R12);
            [D12, I] = sort(real(diag(D)), 'descend');
            [~, M12] = power_truncate(D12, POWER_TRUNCATION_THRESHOLD);
            U12 = U(:, I(1:M12));
            D12 = D12(1:M12);
            rho12 = mean(D12);
            %
            R21 = U1 * U1' * Rue{2} * (U1 * U1');
            [U, D] = eig(R21);
            [D21, I] = sort(real(diag(D)), 'descend');
            [~, M21] = power_truncate(D21, POWER_TRUNCATION_THRESHOLD);
            U21 = U(:, I(1:M21));
            D21 = D21(1:M21);
            rho21 = mean(D21);

            % Codebook design
            Cb1 = zeros(Mbf, Nt);
            Cb2 = zeros(Mbf, Nt);
            for i = 1:Mbf
                u = (U1 * U1') * Rms{1} * (randn(Nt, 1) + 1i * randn(Nt, 1));
                u = u / norm(u);
                Cb1(i, :) = u.';

                u = (U2 * U2') * Rms{2} * (randn(Nt, 1) + 1i * randn(Nt, 1));
                u = u / norm(u);
                Cb2(i, :) = u.';
            end
            Cb = cell(K, 1);
            Cb{1} = Cb1;
            Cb{2} = Cb2;
            Sbf = 2 ^ (- Bf / Nt);  % Regularization parameter from [Dabbagh & Love 2008]
            Sbf1 = 2 ^ (- Bf / (M1 - 1));
            Sbf2 = 2 ^ (- Bf / (M2 - 1));


            % D2D signaling model
            Bd2d = 20 * log2(1 + P0d2d * dm2m ^ (-2));
            Alpha2 = [1,   1 - geo_mean(D12) / mean(D12) * 2 ^ (- Bd2d / M12)
                      1 - geo_mean(D21) / mean(D21) * 2 ^ (- Bd2d / M21),       1];

            % Parameters for the DR scheme
            Q1 = (1 - Alpha2(1, 2)) * rho21 * (U21 * U21') + K / P * eye(Nt);
            Q2 = (1 - Alpha2(2, 1)) * rho12 * (U12 * U12') + K / P * eye(Nt);
            BigB = [1 - Sbf1,    0
                    0,      1 - Sbf2];
            Phi1 = (Alpha2(1, 2)^2 * trace(R21) * R21 + Q1 * (U1 * U1') * Q1 ...
                    + Alpha2(1, 2) * (R21 * (U1 * U1') * Q1 + Q1 * (U1 * U1') * R21)) / M1;
            Phi2 = (Alpha2(2, 1)^2 * trace(R12) * R12 + Q2 * (U2 * U2') * Q2 ...
                    + Alpha2(2, 1) * (R12 * (U2 * U2') * Q2 + Q2 * (U2 * U2') * R12)) / M2;
            BigPhi = Phi1 * Sbf1 + Phi2 * Sbf2;    

            % For test
            Q11 = ((1 - Alpha2(1, 2)) + K / P) * eye(Nt);
            Q22 = ((1 - Alpha2(2, 1)) + K / P) * eye(Nt);

            % Variables for fast computation
            CQ1 = zeros(Mbf, 1);
            CQ2 = zeros(Mbf, 1);
            for i = 1:Mbf
                CQ1(i) = real(conj(Cb1(i, :)) * Q1 * Cb1(i, :).');
                CQ2(i) = real(conj(Cb2(i, :)) * Q2 * Cb2(i, :).');
            end

            % ====================================================================
        end
    
        % CSI model 
        H = zeros(Nt, K);
        for k = 1:K
            hk = Rms{k} * (randn(Nt, 1) + 1i * randn(Nt, 1)) / 2;
            H(:, k) = hk;
        end
        % CSI Exchange model
        Hd = cell(K, 1);
        Hd{1} = H;
        Hd{2} = H;
        % -- Model h_2^1
        Xi21 = sqrt(rho21) * U21 * (randn(M21, 1) + 1i * randn(M21, 1)) / sqrt(2);
        Hd{1}(:, 2) = sqrt(Alpha2(1, 2)) * (U1 * U1') * H(:, 2) + sqrt(1 - Alpha2(1, 2)) * Xi21;
        % -- Model h_1^2
        Xi12 = sqrt(rho12) * U12 * (randn(M12, 1) + 1i * randn(M12, 1)) / sqrt(2);
        Hd{2}(:, 1) = sqrt(Alpha2(2, 1)) * (U2 * U2') * H(:, 1) + sqrt(1 - Alpha2(2, 1)) * Xi12;
       
        % Gain = diag(H' * H);
        
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


        % D2D robust (Proposed dual regularized scheme) -------------------
        H2 = zeros(Nt, K);  % FEEDBACK matrix
        % -- UE1 --
        A1 = abs(conj(Cb1) * Hd{1}(:, 1)) .^ 2;
        B1 = Alpha2(1, 2) * abs(conj(Cb1) * Hd{1}(:, 2)) .^ 2 + CQ1;
        C1 = A1 ./ B1;
        [~, imk] = max(C1);
        H2(:, 1) = Cb1(imk, :).';
        % -- UE 2 --
        A2 = abs(conj(Cb2) * Hd{2}(:, 2)) .^ 2;
        B2 = Alpha2(2, 1) * abs(conj(Cb2) * Hd{2}(:, 1)) .^ 2 + CQ2;
        C2 = A2 ./ B2;
        [~, imk] = max(C2);
        H2(:, 2) = Cb2(imk, :).';

        BigU = [Q1 * H2(:, 1), Q2 * H2(:, 2)];
        W4 = ((BigU * BigB * BigU') + ((Sbf1 * Alpha2(1, 2)^2  + Sbf2 * Alpha2(2, 1)^2) + K / P) * eye(Nt)) \ H2;
        [rsum(2, imc), rates(:, 2, imc)] = sum_rate(H, W4, P);

        
        % D2D non-robust --------------------------------------------------
        W5 = zeros(Nt, K);
        for k = 1:K
            Hk = Hd{k};
            Q = abs(conj(Cb{k}) * Hk) .^ 2;
            Qsum = sum(Q, 2);
            Qsum_alpha = Qsum + ones(Mbf, 1) * (alpha);
            Qk = Q(:, k) ./ (Qsum_alpha - Q(:, k));
            [~, imk] = max(Qk);
            W5(:, k) = Cb{k}(imk, :).';
        end
        [rsum(3, imc), rates(:, 3, imc)] = sum_rate_neql(H, W5, P);
        
        
        
        % BS precoding based on channel feedback --------------------------
        Hfb = zeros(Nt, 2);
        Qb = abs(conj(Cb1) * H(:, 1)) .^ 2;
        [~, I] = max(Qb);
        Hfb(:, 1) = Cb1(I, :).';   % Feedback only the direction
        %
        Qb = abs(conj(Cb2) * H(:, 2)) .^ 2;
        [~, I] = max(Qb);
        Hfb(:, 2) = Cb2(I, :).';   % Feedback only the direction
        

        % BS Robust MMSE
        W6 = Hfb / (Hfb' * Hfb + (alpha + K * Sbf) * eye(K));
        [rsum(4, imc), rates(:, 4, imc)] = sum_rate_neql(H, W6, P);
    end
    
    Rsum(:, i_case) = mean(rsum, 2);
end

Rsum0 = Rsum(:, 1:Ncase0);
Rsum1 = Rsum(:, Ncase0 + 1: 2 * Ncase0);

%% Figures
my_markers = {'d', 's', '^', 'o', '<', '>', 'v', '+', 'p', '.', '*'}.';
my_colors4 = [0 0 1 %'b'
              0    0.4980         0 %'g'
              1 0 0 %'r'
              0    0.7490    0.7490];%'c'
my_colors3 = my_colors4(2:end, :);
Nschemes_to_show = length(schemes_to_show);

hf = figure(testFigID);

% First set of curves
set(gca, 'ColorOrder', my_colors3, 'NextPlot', 'replacechildren');
p_handle = plot(Dm2m_vec, Rsum0(schemes_to_show, Nsnr + 1:end), 'LineWidth', 2, 'MarkerSize', 9);
set(p_handle, {'marker'}, my_markers(1:Nschemes_to_show));

% Second set of curves
hold on
p_handle = plot(Dm2m_vec, Rsum1(schemes_to_show, Nsnr + 1:end), '--', 'LineWidth', 2, 'MarkerSize', 9);
set(p_handle, {'marker'}, my_markers(1:Nschemes_to_show));
set(gca, 'FontSize', 14);
legend(str_sch{schemes_to_show}, 'location', 'Southeast');
ylabel('Sum rate (bps/Hz)');
xlabel('Inter user distance (m)');
% title(sprintf('Nt = %d, K = %d, SNR = %d dB', Nt, K, round(SNR_vec(end))));
ylim([0 16]);
