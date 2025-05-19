clear all;
close all;

% Input bit stream
xbit = [1 0 1 1 0 1 0 0 0 1 1 0];

% Differential Encoding
difencod(1) = ~(1 - xbit(1));
for i = 2:length(xbit)
    difencod(i) = ~(difencod(i-1) - xbit(i));
end

% Differential Decoding
xbit(1) = 1 - ~(difencod(1));
for i = 2:length(xbit)
    xbit(i) = difencod(i-1) - ~(difencod(i));
    if xbit(i) == -1
        xbit(i) = 1;
    end
end

% In-phase and quadrature components (unipolar)
for i = 1:2:(length(difencod)-1)
    inp(i) = difencod(i);
    inp(i+1) = inp(i);
end

for i = 2:2:length(difencod)
    qp(i) = difencod(i);
    qp(i-1) = qp(i);
end

% Bipolar NRZ
for i = 1:length(inp)
    it(i) = 2 * inp(i) - 1;
end
for i = 1:length(qp)
    qt(i) = 2 * qp(i) - 1;
end

% Raised cosine filter
filtorder = 40;
nsamp = 4;
delay = filtorder / (2 * nsamp);
rolloff = 0.5;
rrcfilter = rcosine(1, nsamp, 'fir/normal', rolloff, delay);

figure(1);
impz(rrcfilter, 1);
grid on;
title('Impulse Response of Raised Cosine Filter');

% Transmit signal generation
itx = rcosflt(it, 1, nsamp, 'filter', rrcfilter);
Drate = 64000;
T = 1 / Drate;
Ts = T / nsamp;
time = 0:Ts:(length(itx)-1)*Ts;

figure(2);
plot(time, itx);
xlabel('Time (sec)');
ylabel('Amplitude (V)');
title('Filtered In-Phase Component');
grid on;

qtx = rcosflt(qt, 1, nsamp, 'filter', rrcfilter);
tme = Ts:Ts:(length(itx)-1)*Ts + Ts;

figure(3);
plot(tme, qtx);
xlabel('Time (sec)');
ylabel('Amplitude (V)');
title('Filtered Quadrature Component');
grid on;

fc = 900e6;
dd = 2 * pi * fc * time';
ddd = 2 * pi * fc * tme';

% OQPSK modulation with delay
oqdelay = zeros(1, length(qtx));
oqdelay(nsamp+1:end) = qtx(1:end-nsamp);
mt = cos(dd) .* itx + sin(ddd) .* oqdelay';

figure(4);
plot(time, mt);
xlabel('Time (sec)');
ylabel('Amplitude (V)');
title('OQPSK Modulated Signal');
grid on;

% Add AWGN
snr = 10;
madd = awgn(mt, snr);

figure(5);
plot(time, madd);
xlabel('Time (sec)');
ylabel('Amplitude (V)');
title('OQPSK Signal with AWGN');
grid on;

% Coherent demodulation
cscomp = madd .* cos(dd);
sincomp = madd .* sin(ddd);

lpfin = rcosflt(cscomp, 1, nsamp, 'filter', rrcfilter);
lpfqu = rcosflt(sincomp, 1, nsamp, 'filter', rrcfilter);

tmx = 0:Ts:(length(lpfin)-1)*Ts;
tmy = Ts:Ts:(length(lpfqu)-1)*Ts + Ts;

figure(6);
subplot(2,1,1);
plot(tmx, lpfin);
xlabel('Time (sec)');
ylabel('Amplitude');
title('Demodulated In-Phase Signal');
grid on;

subplot(2,1,2);
plot(tmy, lpfqu);
xlabel('Time (sec)');
ylabel('Amplitude (V)');
title('Demodulated Quadrature Signal');
grid on;

% Sampling and decision
half = filtorder / 2;
itxx = itx(half:nsamp:length(xbit)*nsamp + half - 1);
for i = 1:length(itxx)
    chk1(i) = sign(itxx(i));
end

ityy = qtx(half:nsamp:length(xbit)*nsamp + half - 1);
for i = 1:length(ityy)
    chk2(i) = sign(ityy(i));
end

disp('I channel bit stream checking');
distortion = sum((it - chk1).^2) / length(chk1);
disp(distortion);

disp('Q channel bit stream checking');
distortion = sum((qt - chk2).^2) / length(chk2);
disp(distortion);

% Reconstruct differentially decoded bit stream
for i = 1:2:length(xbit)
    dfd(i) = chk1(i);
end
for i = 2:2:length(xbit)
    dfd(i) = chk2(i);
end
for i = 1:length(xbit)
    dfdecod(i) = (dfd(i) == 1);
end

detected(1) = 1 - ~dfdecod(1);
for i = 2:length(xbit)
    detected(i) = dfdecod(i-1) - ~dfdecod(i);
    if detected(i) == -1
        detected(i) = 1;
    end
end

disp('Distortion between transmitted and received NRZ bit stream');
distortion = sum((xbit - detected).^2) / length(detected);
disp(distortion);

tmx = 0:(1/64000):(1/64000)*(length(xbit)-1);
figure(7);
subplot(2,1,1);
stairs(tmx, xbit);
xlabel('Time (sec)');
ylabel('Binary value');
title('Transmitted Bit Stream');
set(gca, 'ytick', [0 1]);
grid on;

subplot(2,1,2);
stairs(tmx, detected);
xlabel('Time (sec)');
ylabel('Binary value');
title('Received Bit Stream');
set(gca, 'ytick', [0 1]);
grid on;

