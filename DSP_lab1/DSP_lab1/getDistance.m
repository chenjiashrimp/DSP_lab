f_start = 100e6;         % starting frequency of the sweep
f_stop = 250e6;          % stopping frequency of the sweep
sweep_time = 7e-6;       % time to sweep from f_start to f_stop
R = 100;                 % distance between target and radar (in meters)
phi = pi/4;              % initial phase of the signal
K = 0.5;                 % fading factor
Fs = 1.2e12;                % sampling frequency
c = 3e8;
B = f_stop-f_start;
S = B/sweep_time;
fc = 60e9;

% Generate LFMCW signal
t = linspace(0, sweep_time, sweep_time*Fs);
f = linspace(f_start, f_stop, sweep_time*Fs);
chirp_signal = cos(2*pi.*(0.5*S.*t.^2 + fc.*t) + phi);

% Generate radar echo signal
tau = 2*R/c;             % round-trip time
echo_signal= cos(2*pi.*(0.5*S.*(t-tau).^2 + fc.*(t-tau)) + phi);

mixed_signal = chirp_signal.*echo_signal;

Nfft = 2^nextpow2(length(mixed_signal));   % FFT length
spectrum = fftshift(fft(mixed_signal,Nfft));   % FFT of mixed signal
freq = linspace(-Fs/2,Fs/2,length(spectrum));    % frequency vector for plotting

% Find peak position of frequency spectrum
[~,idx] = max(abs(spectrum));   % index of peak magnitude
delta_f = freq(idx);    % frequency

R1 = abs(c*sweep_time*delta_f/(2*B)); %get R
fprintf('R1 = %f, R = %f',R1,R);