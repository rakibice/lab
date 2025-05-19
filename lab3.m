clear all;
close all;

msg = round(rand(1,1000));
trellis = poly2trellis(3,[6 7]);
user = convenc(msg, trellis);
user(user==0) = -1;
length_user = length(user);

fc = 5000;
eb = 0.5;
bitrate = 1000;
tb = 1/bitrate;
chiprate = 10000;
tc = 1/chiprate;
t = tc:tc:tb*length_user;

basebandsig = [];
for i = 1:length_user
    basebandsig = [basebandsig repmat(user(i), 1, tb/tc)];
end

figure(1)
stairs(t(1:800), basebandsig(1:800))
xlabel('Time (sec)')
ylabel('Binary value')
set(gca,'ytick',[-1 1])
title('Original binary sequence (bipolar NRZ)')

bpskmod = [];
for i = 1:length_user
    for j = tc:tc:tb
        bpskmod = [bpskmod sqrt(2*eb)*user(i)*cos(2*pi*fc*j)];
    end
end

number = length(t);
spectrum = abs(fft(bpskmod));
sampling_frequency = 2 * fc;
sampling_interval = 1 / sampling_frequency;

for i = 1:number
    frequency(i) = (1 / (number * sampling_interval)) * i;
end

figure(2)
plot(frequency, spectrum)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Frequency spectrum of BPSK modulated signal')
grid on

seed = [1 -1 1 -1];
pn = [];
for i = 1:length_user
    for j = 1:10
        pn = [pn seed(4)];
        temp = -1 * (seed(4) == seed(3)) + 1 * (seed(4) ~= seed(3));
        seed = [temp seed(1:3)];
    end
end

pnupsampled = [];
for i = 1:length(pn)
    pnupsampled = [pnupsampled repmat(pn(i), 1, tb/(10*tc))];
end

sigtx = bpskmod .* pnupsampled;

figure(3)
plot(t(1:200), sigtx(1:200))
xlabel('Time (sec)')
ylabel('Amplitude')
title('Transmitted DS-CDMA signal segment')
grid on

chan = ricianchan(1/chiprate, 100, 15);
chan.ResetBeforeFiltering = 0;
fad = abs(filter(chan, ones(size(sigtx))));
fadedsig = fad .* sigtx;

snr_in_dBs = 0:1:10;
for m = 1:length(snr_in_dBs)
    composite_signal = awgn(fadedsig, snr_in_dBs(m), 'measured');
    rx = composite_signal .* pnupsampled;

    demodcar = [];
    for i = 1:length_user
        for j = tc:tc:tb
            demodcar = [demodcar sqrt(2*eb)*cos(2*pi*fc*j)];
        end
    end

    bpskdemod = rx .* demodcar;
    len_dmod = length(bpskdemod);
    sumval = zeros(1, len_dmod / 10);

    for i = 1:length(sumval)
        sumval(i) = sum(bpskdemod((i-1)*10+1:i*10));
    end

    rxbits = sumval > 0;
    tblen = 3;
    delay = tblen;
    decoded = vitdec(rxbits, trellis, tblen, 'cont', 'hard');
    [~, ber(m)] = biterr(decoded(delay+1:end), msg(1:end-delay));
end

figure(4)
plot(snr_in_dBs, ber, '-o')
xlabel('Signal to Noise Ratio (dB)')
ylabel('Bit Error Rate (BER)')
legend('Simulated BER')
title('BER vs. SNR for 1/2 Convolutionally Encoded DS-CDMA (AWGN + Rician)')
grid on
