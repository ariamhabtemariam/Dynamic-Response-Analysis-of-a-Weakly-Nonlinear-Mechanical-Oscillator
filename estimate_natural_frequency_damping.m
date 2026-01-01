%% Load and Setup
load('freeoscitime.mat');       % loads time vector 't'
load('freeoscidata.mat');       % loads signal, assumed as 'data'


dt = mean(diff(t));    
fs = 1 / dt;
fprintf('Estimated sampling rate: %.2f Hz\n', fs);


%% Plot Raw Signal
figure;
plot(t, data);
xlabel('Time (s)');
ylabel('Amplitude');
title('Raw Signal');
grid on;


%% Bandpass Filter (20–40 kHz)
low_cutoff = 2e4;       
high_cutoff = 4e4;      
filter_order = 4;


[b, a] = butter(filter_order, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
data_filtered = filtfilt(b, a, data);  
% Plot filtered signal
figure;
plot(t, data_filtered);
xlabel('Time (s)');
ylabel('Amplitude');
title('Bandpass Filtered Signal (20–40 kHz)');
grid on;


%% --- Pulse and Echo Detection ---
pulserate = 25; % Hz
pulsewidth = 0.00032;         
pulsewidthsteps = round(pulsewidth * fs);   
pulseThreshold = 2.5;        
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
distances = ((times * c) / 2) * 100;


%% Plot Distance vs Time
figure;
plot(t(indexes), distances, 'b-o');
xlabel('Time (s)');
ylabel('Distance (cm)');
title('Distance vs Time of Pulse Emission');
grid on;


%% FFT of Filtered Signal
N = length(data_filtered);
data_detrended = detrend(data_filtered);   
window = hann(N);                          
data_windowed = data_detrended .* window;  


Y = fft(data_windowed);                    
amplitude = abs(Y)/(sum(window)/2);        


f = (0:N-1)*(fs/N);                        
half_N = floor(N/2);
f_plot = f(1:half_N);
amplitude_plot = amplitude(1:half_N);


figure;
plot(f_plot, amplitude_plot, 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('FFT Spectrum of Filtered Signal');
xlim([0 fs/2]);
grid on;


%% Full Response Plot
figure;
plot(t(indexes), distances);
xlabel('Time (s)');
ylabel('Distance (cm)');
title('Response (Full View)');
grid on;


%% Estimate Steady-State Value
y_inf = mean(distances(end-100:end));  % Approximate steady-state distance
figure;
plot(t(indexes), distances);
xlabel('Time (s)');
ylabel('Amplitude');
title('Zoomed Response for Peak Selection');
grid on;
%% Manual Peak Selection
disp('Select peaks in the zoomed region...');
[tpks, ypks] = ginput(3);   % Click on visible peaks
y0 = ypks(1);
ln_ratios = log((ypks - y_inf) / (y0 - y_inf));


% Logarithmic decrement plot
figure;
plot(tpks, ln_ratios, 'o-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('ln((y(t) - y∞)/(y0 - y∞))');
title('Logarithmic Decrement Plot');
grid on;


% Linear fit to extract damping ratio * natural frequency
coeffs = polyfit(tpks, ln_ratios, 1);
slope = coeffs(1);
fprintf('Slope = %.4f, so ζω_n = %.4f rad/s\n', slope, -slope);


% Estimate damped frequency from time between peaks
Td = mean(diff(tpks));         % Period of damped oscillation
omega_d = 2 * pi / Td;         % Damped natural frequency (rad/s)
fprintf('Estimated ω_d = %.2f rad/s\n', omega_d);
omega_n = sqrt(omega_d^2 + slope^2);
zeta = -slope / omega_n;


fprintf('Estimated ω_n = %.2f rad/s\n', omega_n);
fprintf('Estimated ζ = %.4f (damping ratio)\n', zeta);
