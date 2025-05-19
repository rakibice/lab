clear all;
close all;

f = 1000;
Fs = 4000;
t = 1/Fs:1/Fs:1;
Am = 1.0;
signal = Am * sin(2 * pi * f * t);

figure(1);
plot(t(1:200), signal(1:200));
set(gca, 'ytick', [-1.0 0 1.0]);
title('A segment of synthetically generated sinusoidal waveform');
grid on;
xlabel('Time (sec)');
ylabel('Amplitude (volt)');

maximumvalue = max(signal);
minimumvalue = min(signal);
interval = (maximumvalue - minimumvalue) / 255;
partition = minimumvalue:interval:maximumvalue;
codebook = (minimumvalue - interval):interval:maximumvalue;
[index, quants, distor] = quantiz(signal, partition, codebook);

indxtrn = index';
for i = 1:4000
    matrix(i, 1:8) = bitget(uint8(indxtrn(i)), 1:8);
end

matrixtps = matrix';
baseband = reshape(matrixtps, 4000 * 8, 1);
Tb = 1 / 32000;
time = 0:Tb:1;

figure(2);
stairs(time(1:500), baseband(1:500));
title('A segment of baseband signal');
xlabel('Time (sec)');
ylabel('Binary value');
set(gca, 'ytick', [0 1]);
axis([0, time(500), 0, 1]);

input_to_Convolutional_encoder = baseband';
trellis = poly2trellis(7, [171 133]);
code = convenc(input_to_Convolutional_encoder, trellis);

st2 = 4831;
data_interleave = randintrlv(code, st2);

M = 4;
k = log2(M);
symbol = bi2de(reshape(data_interleave, k, length(data_interleave)/k).', 'left-msb');
modulated_data = qammod(symbol, M);
demodulated_data = qamdemod(modulated_data, M);

[number, ratio] = symerr(symbol, demodulated_data);

Retrieved_bit = de2bi(demodulated_data, 'left-msb');
Retrieved_bit = Retrieved_bit';
Retrieved_bit = reshape(Retrieved_bit, 64000, 1);

errors = zeros(size(Retrieved_bit));
inter_err = bitxor(Retrieved_bit, errors);
data_deinterleave = randdeintrlv(inter_err, st2);

tblen = 3;
decodx = vitdec(data_deinterleave, trellis, tblen, 'cont', 'hard');
N3 = length(decodx);
decod2(1:N3-3) = decodx(tblen+1:end);
decod2(N3) = decodx(1);
decod2 = decod2';

[number, ratio] = biterr(decod2, baseband);

convert = reshape(decod2, 8, 4000);
[number, ratio] = biterr(convert, matrixtps);

convert = convert';
intconv = bi2de(convert);
[number, ratio] = biterr(intconv, index');
sample_value = minimumvalue + intconv .* interval;

figure(3);
subplot(2,1,1);
plot(time(1:100), signal(1:100));
set(gca, 'ytick', [-1.0 0 1.0]);
axis([0, time(100), -1, 1]);
title('Graph for a segment of recorded Audio signal');
xlabel('Time (sec)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(time(1:100), sample_value(1:100));
axis([0, time(100), -1, 1]);
set(gca, 'ytick', [-1.0 0 1.0]);
title('Graph for a segment of retrieved Audio signal');
xlabel('Time (sec)');
ylabel('Amplitude');
grid on;
