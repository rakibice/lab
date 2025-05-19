clear all;
close all;

% Generate sinusoidal signal
f = 1000;
Fs = 4000;
t = 1/Fs : 1/Fs : 1;
Am = 1.0;
signal = Am * sin(2*pi*f*t);

figure(1);
plot(t(1:200), signal(1:200));
set(gca, 'ytick', [-1.0 0 1.0]);
title('Segment of Original Sinusoidal Signal');
xlabel('Time (s)');
ylabel('Amplitude (V)');
grid on;

% Quantization
maxval = max(signal);
minval = min(signal);
interval = (maxval - minval) / 255;
partition = minval : interval : maxval;
codebook = (minval - interval) : interval : maxval;
[index, ~, ~] = quantiz(signal, partition, codebook);

% Decimal to binary conversion
indxtrn = index';
for i = 1:4000
    matrix(i, :) = bitget(uint8(indxtrn(i)), 1:8);
end
matrixtps = matrix';
baseband = reshape(matrixtps, [], 1);

Tb = 1 / 32000;
time = 0 : Tb : 1;

figure(2);
stairs(time(1:500), baseband(1:500));
title('Segment of Baseband Signal');
xlabel('Time (s)');
ylabel('Binary Value');
set(gca, 'ytick', [0 1]);
axis([0 time(500) 0 1]);

% FEC Encoding using Convolutional Code
input_data = baseband';
trellis = poly2trellis(7, [171 133]);
encoded = convenc(input_data, trellis);

% Interleaving
seed = 4831;
interleaved = randintrlv(encoded, seed);

% QPSK Modulation
M = 4;
k = log2(M);
symbols = bi2de(reshape(interleaved, k, []).', 'left-msb');
modulated = pskmod(symbols, M);

% QPSK Demodulation
demodulated = pskdemod(modulated, M);
[~, sym_err_ratio] = symerr(symbols, demodulated);

% Binary Conversion
retrieved_bits = de2bi(demodulated, 'left-msb')';
retrieved_bits = reshape(retrieved_bits, [], 1);

% Deinterleaving
interrupted = bitxor(retrieved_bits, zeros(size(retrieved_bits)));
deinterleaved = randdeintrlv(interrupted, seed);

% Viterbi Decoding
tblen = 3;
decoded = vitdec(deinterleaved, trellis, tblen, 'cont', 'hard');
decoded2(1:(length(decoded) - 3)) = decoded(tblen+1:end);
decoded2(length(decoded)) = decoded(1);
decoded2 = decoded2';

% Bit Error Calculation
baseband = double(baseband);
[~, biterr_ratio1] = biterr(decoded2, baseband);

% Reshape and convert back
reshaped = reshape(decoded2, 8, 4000);
matrixtps = double(matrixtps);
[~, biterr_ratio2] = biterr(reshaped, matrixtps);

reshaped = reshaped';
intconv = bi2de(reshaped);
[~, biterr_ratio3] = biterr(intconv, index');

% Reconstruct Signal
reconstructed = minval + intconv .* interval;

figure(3);
subplot(2,1,1);
plot(time(1:100), signal(1:100));
set(gca, 'ytick', [-1.0 0 1.0]);
axis([0 time(100) -1 1]);
title('Segment of Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(time(1:100), reconstructed(1:100));
set(gca, 'ytick', [-1.0 0 1.0]);
axis([0 time(100) -1 1]);
title('Segment of Reconstructed Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
