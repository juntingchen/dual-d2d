function [r, R, Sinr, Sig, Int] = sum_rate(H, W0, P)
%
% Evaluate the sum rate
% Channel model: Y = H' * W * Phi + N, where N is the unit
% variance Gaussian noise, and Phi is to make sure that equal power
% allocation is applied with sum power P.
%
% Copyright (c), CHEN Junting, eejtchen@connect.ust.hk

K = size(H, 2);

W = W0 * diag(1 ./ sqrt(diag(W0' * W0))) * sqrt(P / K);
Y = H' * W;
Y2 = abs(Y) .^ 2;
Y2sum = sum(Y2, 2) + ones(K, 1);

R = zeros(K, 1);
Sinr = zeros(K, 1);
Sig = zeros(K, 1);
Int = zeros(K, 1);
for k = 1:K
    Sig(k) = Y2(k, k);
    Int(k) = Y2sum(k) - Y2(k, k) - 1;
    
    Sinr(k) = Y2(k, k) / (Y2sum(k) - Y2(k, k));
    R(k) = log2(1 + Sinr(k));
end
r = sum(R);