function [r, R, Sinr] = sum_rate_neql(H, W0, P)
%
% Evaluate the sum rate
% Channel model: Y = H' * W * Phi + N, where N is the unit
% variance Gaussian noise,and Phi is to make sure to meet the sum power
% constraint. 
%
% NOTE: non equal power allocation is applied.
%
% Copyright (c), CHEN Junting, eejtchen@connect.ust.hk

K = size(H, 2);

W = W0  / norm(W0, 'fro') * sqrt(P);
Y = H' * W;
Y2 = abs(Y) .^ 2;
Y2sum = sum(Y2, 2) + ones(K, 1);

R = zeros(K, 1);
Sinr = zeros(K, 1);
for k = 1:K
    Sinr(k) = Y2(k, k) / (Y2sum(k) - Y2(k, k));
    R(k) = log2(1 + Sinr(k));
end
r = sum(R);