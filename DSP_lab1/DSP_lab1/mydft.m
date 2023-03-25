function y = mydft(x)
% Computes the Discrete Fourier Transform (DFT) of a sequence x
% Input:
%     x: input sequence of length N
% Output:
%     y: DFT of x, also of length N

N = length(x);
y = zeros(1, N);

for k = 0:N-1
    for n = 0:N-1
        y(k+1) = y(k+1) + x(n+1) * exp(-1i*2*pi*k*n/N);
    end
end