
x = randn(2,1024);
fs = 256;

[Nc,Ns] = size(x);

Nw = 8;
Nsw = Ns/Nw;

% GOLD STANDARD (Magnitude squared coherence)
% [Cxy,F] = mscohere(X,Y,WINDOW,NOVERLAP,NFFT,Fs)
[refCxy,f] = mscohere(x(1,:),x(2,:),hamming(Nsw),0.25*Nsw,Nsw,fs);
Nf = numel(f);

% TEST FUNCTIONS
% Spectrogram
cfg.band        = [];                  % Frequency band
cfg.length      = Nsw;       % Segments length
cfg.overlap     = 0.25;              % Segments overlap
cfg.window      = 'hamming';                   % Tapering window name
cfg.preprocess  = 'none';               % Segments data preprocess
cfg.nfft        = Nsw; % FFT computation length
[TF,F,T] = process_sc_spectrogram('Spectrogram',x, cfg);
% Coherency
X = process_sc_connect_spectral('Cross',TF);
tstCxy = process_sc_connect_spectral('Coherency',X);
tstCxy = squeeze(tstCxy);
tstCxy = abs(tstCxy(2,:)).^2; % Magnitude squared coherence


% COMPARISON
N = min(numel(refCxy),numel(tstCxy));
refCxy = refCxy(1:N); 
tstCxy = tstCxy(1:N);
[refCxy(:), tstCxy(:)]