import mydft.*

% Generate a random input sequence of length 128
x = randn(1, 128);

% Compute DFT using mydft function
y1 = mydft(x);

% Compute DFT using fft function
y2 = fft(x);

% Compare the results using the norm of the difference
error = norm(y1 - y2);

% Display the error
fprintf('Error between mydft and fft: %e\n', error);