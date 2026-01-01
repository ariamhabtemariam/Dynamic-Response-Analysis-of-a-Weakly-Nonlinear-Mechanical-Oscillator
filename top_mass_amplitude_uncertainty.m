load('10increasing.mat');  % loads time vector 't'
load('10timeincreasing.mat');  % loads signal as 'data'


% Estimate sampling interval
dt = mean(diff(t));    
fs = 1 / dt;
fprintf('Estimated sampling rate: %.2f Hz\n', fs);


% Plot raw data
figure;
plot(t, data);
xlabel('Time (s)');
ylabel('Amplitude');
title('Raw Signal');
grid on;


%% Bandpass Filter
low_cutoff = 1.5e4;       
high_cutoff = 3.5e4;      
filter_order = 2;         


[b, a] = butter(filter_order, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
data_filtered = filtfilt(b, a, data);   % Zero-phase filtering


%  plot filtered signal
figure;
plot(t, data_filtered);
xlabel('Time (s)');
ylabel('Amplitude');
title('Bandpass Filtered Signal (25–35 kHz)');
grid on;


%% Pulse and Echo Detection on Filtered Data
pulserate = 25; % Hz
pulsewidth = 0.00032;            
pulsewidthsteps = round(pulsewidth * fs);   
pulseThreshold = 2;        
echoThreshold = 0.23;       
stepsbetweenpulses = round((1 / pulserate) * fs * 0.5);


c = 343;  % Speed of sound in air (m/s)


indexes = []; echoindexes = []; times = [];


i = 1;
while i < length(data_filtered)
   if abs(data_filtered(i)) > pulseThreshold
       j_start = i + pulsewidthsteps;
       j_end = min(j_start + stepsbetweenpulses, length(data_filtered));
       for j = j_start:j_end
           if abs(data_filtered(j)) > echoThreshold
               indexes(end+1) = i;     
               echoindexes(end+1) = j;
               times(end+1) = (t(j) - t(i));
               i = j_end;                  
               break
           end
       end
   end
   i = i + 1;
end


% Convert time to distance in cm
pulse_amplitudes = ((times * c) / 2) * 100;


% Plotting distance results
figure;
plot(t(indexes), pulse_amplitudes, 'b-');
xlabel('Time (s)');
ylabel('Distance (cm)');
title('Distance vs Time');
grid on;


% Extract the uniformly spaced time vector for interpolation
t_uniform = linspace(t(indexes(1)), t(indexes(end)), 1024);  % You can adjust 1024 as needed


% Interpolate the distance data onto the uniform grid
pulse_distances_interp = interp1(t(indexes), pulse_amplitudes, t_uniform, 'linear', 'extrap');


% Apply FFT
N = length(pulse_distances_interp);
Y = fft(pulse_distances_interp);
amplitude = abs(Y) / (N/2);
f = (0:N-1)*(1/(t_uniform(2)-t_uniform(1))/N);


% Only keep first half for real-valued signals
half_N = floor(N/2);
f_plot = f(1:half_N);
amplitude_plot = amplitude(1:half_N);


% Plot the FFT result
figure;
plot(f_plot, amplitude_plot, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('FFT of Distance vs Time ');
grid on;
xlim([0, max(f_plot)]);
% Ignore the DC component (0 Hz) for peak frequency detection
[~, peak_idx] = max(amplitude_plot(2:end));  % Skip the first element (DC)
peak_idx = peak_idx + 1;  % Adjust index because we skipped the first


% Find corresponding frequency
dominant_frequency = f_plot(peak_idx);


% Display the result
fprintf('Dominant Frequency (excluding DC): %.2f Hz\n', dominant_frequency);


% Extract time and distance vectors from pulse detection
time_signal = t(indexes);
distance_signal = pulse_amplitudes;


% Smooth signal (optional if there's noise)
distance_signal_smooth = movmean(distance_signal, 3);


f_estimate = dominant_frequency;  % Hz, rough guess of signal frequency
period_est = 1 / f_estimate;  % in seconds
samples_per_period = round(period_est / mean(diff(time_signal)));


% Loop through each period and grab the max (peak)
peak_vals = [];
peak_times = [];
num_chunks = floor(length(distance_signal_smooth) / samples_per_period);


for k = 1:num_chunks
   idx_start = (k-1)*samples_per_period + 1;
   idx_end = min(k*samples_per_period, length(distance_signal_smooth));
  
   [val, rel_idx] = max(distance_signal_smooth(idx_start:idx_end));
   peak_vals(end+1) = val;
   peak_times(end+1) = time_signal(idx_start + rel_idx - 1);
end


% Frequency from average time between peaks
periods = diff(peak_times);
dominant_frequency = 1 / mean(periods);


% Amplitude = average distance from centerline to peak
signal_mean = mean(distance_signal_smooth);
dominant_amplitude = mean(abs(peak_vals - signal_mean));


% Output
fprintf('Estimated Dominant Frequency: %.2f Hz\n', dominant_frequency);
fprintf('Estimated Amplitude: %.2f cm\n', dominant_amplitude);


deviations = abs(peak_vals - signal_mean);
amplitude_std = std(deviations);
amplitude_uncertainty = amplitude_std / sqrt(length(peak_vals));


fprintf('Amplitude: %.2f cm ± %.2f cm (standard error)\n', dominant_amplitude, amplitude_uncertainty);


period_std = std(periods);
mean_period = mean(periods);
frequency_uncertainty = period_std / (mean_period^2);


fprintf('Estimated Dominant Frequency: %.2f Hz ± %.2f Hz\n', dominant_frequency, frequency_uncertainty);




