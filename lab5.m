clear all;
close all;

% Generate sinusoidal signal
f = 1000;          
Fs = 4000;         
t = 1/Fs : 1/Fs : 1;
Am = 1.0;
signal = Am * sin(2 * pi * f * t);

figure(1);
plot(t(1:200), signal(1:200))
set(gca, 'ytick', [-1.0 0 1.0])
title('A segment of original sinusoidal waveform')
xlabel('Time (sec)')
ylabel('Amplitude (volt)')
grid on

% Quantization
maxval = max(signal);
minval = min(signal);
interval = (maxval - minval) / 255;
partition = minval : interval : maxval;
codebook = (minval - interval) : interval : maxval;
[index, quants, distor] = quantiz(signal, partition, codebook);

% Binary conversion
indxtrn = index';
for i = 1:4000
    matrix(i, 1:8) = bitget(uint8(indxtrn(i)), 1:8);
end
matrixtps = matrix';
baseband = reshape(matrixtps, 4000 * 8, 1);
Tb = 1 / 32000;
time = 0:Tb:1;

figure(2);
stairs(time(1:500), baseband(1:500))
title('A segment of baseband signal')
xlabel('Time (sec)')
ylabel('Binary value')
set(gca, 'ytick', [0 1])
axis([0 time(500) 0 1])

% Convolutional Encoding
input_bits = baseband';
trellis = poly2trellis(7, [171 133]);
coded = convenc(input_bits, trellis);

% Interleaving
seed = 4831;
interleaved = randintrlv(coded, seed);

% BPSK Modulation
M = 2;
k = log2(M);
symbols = bi2de(reshape(interleaved, k, length(interleaved) / k).', 'left-msb');
symbols = double(symbols);
modulated = pskmod(symbols, M);

% BPSK Demodulation
demodulated_symbols = pskdemod(modulated, M);
[~, sym_err_rate] = symerr(symbols, demodulated_symbols);

% Bit recovery
retrieved_bits = de2bi(demodulated_symbols, 'left-msb');
retrieved_bits = retrieved_bits(:);

% Deinterleaving
errors = zeros(size(retrieved_bits));
received = bitxor(retrieved_bits, errors);
deinterleaved = randdeintrlv(received, seed);

% Convolutional Decoding
tblen = 3;
decoded = vitdec(deinterleaved, trellis, tblen, 'cont', 'hard');
decoded_trimmed = decoded(tblen+1:end);
decoded_trimmed(end+1:32000) = decoded(1); 
decoded_trimmed = decoded_trimmed(:);

% Bit error comparison
baseband = double(baseband);
[~, biterr1] = biterr(decoded_trimmed, baseband);

% Bit to byte matrix
converted = reshape(decoded_trimmed, 8, 4000);
[~, biterr2] = biterr(converted, matrixtps);
converted = converted';

% Integer and signal reconstruction
int_vals = bi2de(converted);
[~, biterr3] = biterr(int_vals, index');
reconstructed = minval + int_vals .* interval;

% Plot reconstructed vs. original signal
figure(3);
subplot(2,1,1);
plot(time(1:100), signal(1:100));
title('Segment of original audio signal');
xlabel('Time (sec)');
ylabel('Amplitude');
axis([0 time(100) -1 1]);
set(gca, 'ytick', [-1.0 0 1.0]);
grid on

subplot(2,1,2);
plot(time(1:100), reconstructed(1:100));
title('Segment of retrieved audio signal');
xlabel('Time (sec)');
ylabel('Amplitude');
axis([0 time(100) -1 1]);
set(gca, 'ytick', [-1.0 0 1.0]);
grid on
