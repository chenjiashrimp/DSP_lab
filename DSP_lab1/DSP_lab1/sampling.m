
% 问题2 Define parameters
f_start = 100e6;         % starting frequency of the sweep
f_stop = 250e6;          % stopping frequency of the sweep
sweep_time = 7e-6;       % time to sweep from f_start to f_stop
R = 100;                 % distance between target and radar (in meters)
phi = pi/4;              % initial phase of the signal
K = 0.5;                 % fading factor
Fs = 1e9;                % sampling frequency
c = 3e8;
B = f_stop-f_start;
S = B/sweep_time;
% Generate LFMCW signal
t = linspace(0, sweep_time, sweep_time*Fs);
f = linspace(f_start, f_stop, sweep_time*Fs);
chirp_signal = cos(2*pi.*(0.5*S.*t.^2 + f_start.*t) + phi);

% Generate radar echo signal
tau = 2*R/c;             % round-trip time
echo_signal= cos(2*pi.*(0.5*S.*(t-tau).^2 + f_start.*(t-tau)) + phi);

% Sample signals
num_samples = length(echo_signal);
t_samp = linspace(0, (num_samples-1)/Fs, num_samples);
echo_signal_samp = echo_signal(1:num_samples);
chirp_signal_samp = chirp_signal(1:num_samples);

% Perform DFT transformation
nfft = 2^nextpow2(num_samples);
f_axis = linspace(-Fs/2, Fs/2, nfft);
echo_dft = fftshift(fft(echo_signal_samp, nfft));
chirp_dft = fftshift(fft(chirp_signal_samp, nfft));

% Display time-frequency domain waveform of signal
figure;
subplot(2,1,1); plot(t_samp, chirp_signal_samp);
title('LFMCW Radar Chirp Signal');
xlabel('Time (s)'); ylabel('Amplitude');
subplot(2,1,2); plot(f_axis, abs(chirp_dft));
title('Time-Frequency Domain');
xlabel('Frequency (Hz)'); ylabel('Amplitude');


% Display time-frequency domain waveform of signal
figure;
subplot(2,1,1); plot(t_samp, echo_signal_samp);
title('LFMCW Radar Echo Signal');
xlabel('Time (s)'); ylabel('Amplitude');
subplot(2,1,2); plot(f_axis, abs(echo_dft));
title('Time-Frequency Domain');
xlabel('Frequency (Hz)'); ylabel('Amplitude');


% 问题3 Get the mixed-signal
mixed_signal = chirp_signal.* echo_signal;
mixed_samp = mixed_signal(1:num_samples);
mixed_dft = fftshift(fft(mixed_signal,nfft));

% Display time-frequency domain waveform of mixed-signal
figure;
subplot(2,1,1); plot(t_samp, mixed_samp);
title('LFMCW Mixed Signal');
xlabel('Time (s)'); ylabel('Amplitude');
subplot(2,1,2); plot(f_axis, abs(mixed_dft));
title('Time-Frequency Domain');
xlabel('Frequency (Hz)'); ylabel('Amplitude');

% change sampling frequency to 1024
new_num_samples = 1024;
mixed_dft_1024 = fftshift(fft(mixed_signal,new_num_samples));
new_f_axis = linspace(-Fs/2, Fs/2, new_num_samples);
new_mixed_samp = mixed_signal(1:new_num_samples);
new_t_samp = linspace(0, (new_num_samples-1)/Fs, new_num_samples);

%display the new waveform of the mixed-signal
figure;
subplot(2,1,1); plot(new_t_samp,new_mixed_samp);
title('new LFMCW Mixed Signal');
xlabel('Time (s)'); ylabel('Amplitude');
subplot(2,1,2); plot(new_f_axis, abs(mixed_dft_1024));
title('Time-Frequency Domain');
xlabel('Frequency (Hz)'); ylabel('Amplitude');

