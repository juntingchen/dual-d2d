function [R, Rm] = genChannCorr_lapcn(Nt, theta_deg, delta_deg)
%
% Generate the channel spatial correlaiton at the TX side. The pathloss
% is assumed to be 1.
%
% PAS MODEL
%   Laplacian (see. TR25.996 V14, Clause 4.5.4)
%   However, we do not consider antenna pattern as in 
%
% INPUT
%   Nt             The number of transmit antennas
%   theta_deg      The azimuth angle of the signal in degree
%   delta_deg      The angular spread in degree
%
% OUTPUT
%   R              Nt * Nt matrix
%   Rm             Square root of R
%
% Copyright (c) 2014-2016, CHEN Junting, eejtchen@connect.ust.hk

% DEBUG MODE
if nargin < 1
    Nt = 180;
    theta_deg = 0;
    % delta_deg = 20;
    mu_as = 1.18;
    eps_as = 0.210;
    delta_deg = 10 .^ (eps_as * randn + mu_as);
    SHOW_PAS = 1;
else
    SHOW_PAS = 0;
end

theta = theta_deg / 180 * pi;
delta = delta_deg / 180 * pi;

sigma = delta;  % We use the same notation as (4.5-2) [TR25.996 V14]

N0 = 1 / integral(@(x) exp(-sqrt(2) * abs(x-theta)/sigma), theta - pi, theta + pi);

% c = 1/(sqrt(2 * pi) * sigma * erf(pi/(sqrt(8) * sigma)));
% c = 1 / erf(pi / (sqrt(2) * sigma));
c = 1;

Corr = zeros(1, Nt);
for d = 0:Nt - 1
    my_fun = @(x) exp( - sqrt(2) * abs(x - theta) / sigma + 1i * pi * d * sin(x)) * N0;
    Corr(d + 1) = c * integral(my_fun, theta - pi, theta + pi);
end

R = zeros(Nt);
for p = 1:Nt
    for q = p:Nt
%         R(p, q) = real(Corr(q-p+1));
        R(p, q) = Corr(q-p+1);
        if p ~= q
            R(q, p) = conj(R(p, q));
        end
    end
end
     
R = R / trace(R) * Nt;

Rm = sqrtm(R);

if SHOW_PAS
    show_pas(R);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function show_pas(R)
%
% A subfunction to show the power spectrum density in angular domain

n = size(R, 1);
Ut = fft(eye(n));

Ra = Ut' * R * Ut;  % R is the covariance matrix, must real on diagonal
Gain = real(diag(Ra)).';

figure(5473),        % A random number
halfid = round(n / 2);
stem([Gain(halfid + 1:end) Gain(1:halfid)]);
end