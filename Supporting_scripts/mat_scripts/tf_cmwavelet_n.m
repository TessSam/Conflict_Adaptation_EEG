function WAVELET=tf_cmwavelet_n2(EEG,varargin)
% @input arguments
% 1, EEG data structure
% 2, WAVELET setting structure
% 
% @ANALYSIS STREAM BOOKKEEPING!!-------------------------------------------
% 1, Get signal and convert it to frequency domain
% 2, Prepare wavelet and convert it to frequency domain
% 3, Convolve 1 and 2 (mutiply element by element)
% 4, Bring the product of step 3 into time domain
% 5, Normalize (same regardless of T-F extraction methods)
% 6, Done!

% @WAVELET BOOKKEEPING!!---------------------------------------------------
%time = time of wavelet
%frex= frequency frequencies of each wavelet in log scale
%s = standard deviation of gaussian function 
% wavelet scaling factor = wavelet frequency / frequency standard.d

%  @ cmorwavf function??? // waveinfo('cmor')------------------------------
% PSI(X) = ((pi*FB)^(-0.5))*exp(2*i*pi*FC*X)*exp(-X^2/FB) 
%     Complex Morlet Wavelet
%  
%     Definition: a complex Morlet wavelet is
%         cmor(x) = (pi*Fb)^{-0.5}*exp(2*i*pi*Fc*x)*exp(-(x^2)/Fb)
%     depending on two parameters:
%         Fb is a bandwidth parameter
%         Fc is a wavelet center frequency
% =========================================================================

% Default Wavelet setting--------------------------------------------------
% EEG(wvtime,frex,wvgsd) = process_options(varargin,...
%     'wvtime',           -2.5:1/EEG.srate:2.5,...
%     'frex',             logspace(log10(min(EEG.allFreqs)),log10(max(EEG.allFreqs)),EEG.numFreqs),...
%     'wvgsd',            [3,10]);

wvtime=-2.5:1/EEG.srate:2.5;
frex=logspace(log10(min(EEG.allFreqs)),log10(max(EEG.allFreqs)),EEG.numFreqs);
wvgsd=[3,10];
    
% Summarize parameters
WAVELET=struct;
WAVELET.wvtime=wvtime;
WAVELET.frex=frex;
WAVELET.wvgsd=wvgsd;

% definte FFT parameters
WAVELET.n_data= EEG.pnts*EEG.trials;
WAVELET.n_wavelet=length(WAVELET.wvtime);
WAVELET.n_convolution= WAVELET.n_wavelet+WAVELET.n_data-1;
WAVELET.n_conv_pow2=pow2(nextpow2(WAVELET.n_convolution));
WAVELET.half_of_wavelet_size = (WAVELET.n_wavelet-1)/2;

if length(WAVELET.wvgsd)==1,wvgsd=repmat(wvgsd,1,2);end
% When wvgsd has an interval of values, the number of periodic cyle in wavelet 
% is adjusted depending on the central frequency to compensate 1/F
% property of EEG signal for high freqeuency. 
% EX) % s gets lower for higher freqs -> more cycles in a wavelet
WAVELET.s=logspace(log10(wvgsd(1)),log10(wvgsd(2)),EEG.numFreqs)./(2*pi*frex);


% Default Analysis setting-------------------------------------------------
WAVELET.note='chan X frequencies X time X trials';
WAVELET.eegpower=zeros(length(EEG.chanlocs),EEG.numFreqs,EEG.pnts,EEG.trials);
WAVELET.eegphase=zeros(length(EEG.chanlocs),EEG.numFreqs,EEG.pnts,EEG.trials);
% WAVELET.eegfilt=zeros(length(EEG.chanlocs),EEG.numFreqs,EEG.pnts,EEG.trials);
% save only whats needed!

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Loop thoroough all channels
ProgressBar('Wavlet Convolution for each channel');

for ch=1:length(EEG.chanlocs)
    %LoopVars
    loop.ch=ch;
    
    %Subset data and Extract FFT of data for each channel at a time
    EEG.actchan=EEG.chanlocs(ch).labels;chInd=strcmpi(EEG.actchan,{EEG.chanlocs.labels});
    eegfft=fft(reshape(EEG.eegraw(chInd,:,:),1,EEG.pnts*EEG.trials),WAVELET.n_conv_pow2);
    
    % Final output
    WAVELET=doConvolve(loop,EEG,WAVELET,eegfft);
    ProgressBar(ch/length(EEG.chanlocs));

end

end

% =========================================================================

function [WAVELET]=doConvolve(loop,EEG,WAVELET,eegfft)
% eegpower = (1=channels, 2=frequency, 3=time) = power 
% eegitpc = (1=channels, 2=frequency, 3=time) = phase clustering value
% eegphase
% -------------------------------------------------------------------------

for fi=1:EEG.numFreqs
    %LoopVars
    loop.fi=fi;
   
    %Get parameters
    wvtime=WAVELET.wvtime;
    frex=WAVELET.frex(fi);
    s=WAVELET.s(fi);
    cpw2=WAVELET.n_conv_pow2;
    
    %Create wavelet, then FFT to put into frequency rep!
    sine_wave=exp(1i*2*pi*frex.*wvtime);
    gaus=exp((-wvtime.^2) ./ (2*s^2));
    %wavelet in time domain=sine_wave.*gaus;
    wavelet_f=fft(sine_wave.*gaus,cpw2);%in frequency domain
    wavelet_f=wavelet_f./max(wavelet_f);% to put back to uV units!!
    
    %old wavelet
    %wavelet_old=fft(sqrt(1/(s*sqrt(pi)))*exp(2*1i*pi*frex.*wvtime).* exp((-wvtime.^2)./(2*s^2)),cpw2);
    %wavelet_old=sqrt(1/(s*sqrt(pi)))*exp(2*1i*pi*frex.*wvtime).* exp((-wvtime.^2)./(2*s^2));

    %%Check wavelet
    %wv2test=wavelet_old;%sine_wave.*gaus/wavelet_old
    %figure;plot(wvtime,real(wv2test))
    %figure;plot3(wvtime,real(wv2test),imag(wv2test),'k')
    %xlabel('Time (ms)'), ylabel('real amplitude'), zlabel('imag amplitude')

    %Perform Convolution(in frequency domain),then Inverse FFT!
    eegconv = ifft(wavelet_f.*eegfft);

    %Cut off edges (length of convolution = data + wavelet -1 )
    eegconv = eegconv(1:WAVELET.n_convolution);
    eegconv = eegconv(WAVELET.half_of_wavelet_size+1:end-WAVELET.half_of_wavelet_size);
    eegconv_r=reshape(eegconv,EEG.pnts,EEG.trials);

   %Structure Data depending on the type of binning!
   [WAVELET]=processData(loop,EEG,WAVELET,eegconv_r);
   
end

end


% =========================================================================

function [WAVELET]=processData(loop,EEG,WAVELET,eegconv_r)
% 1, BIND processing type:-------------------------------------------------
% Each bin is meaninglully conditioned, and power and inter-trial phase 
% dispersion are averaged across trials for each bin. Usually, this processing
% is required for EEG data with binnning for ERP analysis. 
% 2, EPOCH processing type:------------------------------------------------
% A bin (i.e., named as "ALL") containes all usable epochs; therefore,
% power and phase value need to be stored seperately for each epoch. This
% requires epoching and binnig without any categorization of conditions 
% (i.e., must contain all epochs of interest). It is best to use when the output 
% need be combined with behavioral data. Also, good for preserving phases!
% =========================================================================

% initialize
ch=loop.ch;fi=loop.fi;

%PHASE! (raw phase angles)+++++++++++++++++++++++++++++++++++++++++++++++++
%Phase angle: Cannot be averaged over trials!
eegphase=zeros(EEG.pnts,EEG.trials);
eegphase(:,:)=angle(eegconv_r);
WAVELET.eegphase(ch,fi,:,:)=eegphase(:,:);

%POWER!(for each epoch)!+++++++++++++++++++++++++++++++++++++++++++++++++++
temppower=(abs(eegconv_r).^2);%(time,trial)
WAVELET.eegpower(ch,fi,:,:)=temppower;%no-baseline

%FILTERED DATA!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% WAVELET.eegfilt(ch,fi,:,:)=real(eegconv_r);

end

