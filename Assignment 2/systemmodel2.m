%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 5XCC0 Assignment 2 - (C) Pieter Harpe %%%
%%% Only for use at TU/e %%%%%%%%%%%%%%%%%%%%
%%% Do not remove copyright %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here you set all block-level specifications.
% ------------------------------------------------
% Amplifier
Amp_BW         = 6000;         % BW in Hz (min ~6kHz)
Amp_IRN        = 2e-6;         % IRN in Vrms (amplifier 'noise budget')
Amp_InputRange = 20e-3;        % Input range in Vpp (slightly >18mV for margin)
Amp_Zin        = 32e6;         % Input impedance in Ohms (matches system spec)
Amp_Gain       = 56;           % Gain of amplifier in V/V

% Filter
Filter_BW         = 6000;      % BW in Hz (same as amp, or slightly above 6kHz)
Filter_IRN        = 112e-6;    % IRN in Vrms (from 'equal budget': ~G * Amp_IRN)
Filter_InputRange = 1.0;       % Input range in Vpp (amplifier’s output can be ~1Vpp)

% ADC
ADC_fsample    = 12000;        % Sample rate in Hz (>= 2*6kHz per Nyquist)
ADC_IRN        = 112e-6;       % IRN in Vrms (same as filter for equal budget)
ADC_InputRange = 1.0;          % Input range in Vpp (given design constraint)
ADC_N          = 12;           % Resolution of the ADC (bits); ensures ~65dB DR
% NOTE: ENOB ~ N-0.5 = 11.5 => ~6.02*11.5 + 1.76 ~ 71dB max

% Here you set the impedance of the ETIs.
% ------------------------------------------------
Zeti_p = 100e3; % Impedance of the positive electrode [ohm]
Zeti_n = 200e3; % Impedance of the negative electrode [ohm]

% Here you decide which test signals to enable.
% ------------------------------------------------
APon    = 1;  % Enable the AP test signal
CMon    = 0;  % Disable the common-mode disturbance signal
DIFFon  = 0;  % Disable the differential disturbance signal
NOISEon = 1;  % Enable the input-referred noise

% Here you define the parameters for the common-mode disturbance signal.
% ------------------------------------------------
CMfreq = 200;          % Frequency of the common-mode  [Hz]
CMamp  = 300e-3;       % Amplitude [Vpp]

% Here you define the parameters for the differential disturbance signal.
% ------------------------------------------------
DIFFfreq = 125;        % Frequency of the differential disturbance [Hz]
DIFFamp  = 12e-3;      % Amplitude [Vpp]

%%% After this line, you don't need to change anything,
%%% but please read the code if you like to understand what is going on.

% Check Nyquist criterion.
if (2 * Filter_BW > ADC_fsample)
    disp ('Violation of Nyquist criterion');
else
    disp ('Nyquist is very happy with your choice of frequencies');
end

% Load exemplary AP recording.
% Data record 1551.mat from: 
%   Benjamin Metcalfe, "Action potentials recorded from the L5 dorsal rootlet of rat 
%   using a multiple electrode array," Mendeley Data, Version 1, June 12, 2020, 
%   CC BY 4.0 License, doi: 10.17632/ybhwtngzmm.1
load('1551.mat');
APtime = rawdata(:,1)';
AP_p   = 0.5 * rawdata(:,2)' / 10000;  % Compensate for ~80dB gain used in the recording
AP_n   = -AP_p;
nop    = length(AP_p);
k      = 1:nop;
Fsample= 1 / (APtime(2) - APtime(1));
TotalTime = nop / Fsample;

% Create a common-mode disturbance interference signal.
CM_p = 0.5 * CMamp * sin(2*pi*CMfreq * k / Fsample);
CM_n = CM_p;

% Create a differential disturbance signal.
DIFF_p = 0.25 * DIFFamp * sin(2*pi*DIFFfreq * k / Fsample);
DIFF_n = -DIFF_p;

% Add signals together, depending on enable/disable settings.
Signal_p = APon*AP_p + CMon*CM_p + DIFFon*DIFF_p;
Signal_n = APon*AP_n + CMon*CM_n + DIFFon*DIFF_n;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODEL OF AMPLIFIER %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate resistive division due to electrode & amplifier impedances.
IN_p = Signal_p * Amp_Zin / (Zeti_p + Amp_Zin);
IN_n = Signal_n * Amp_Zin / (Zeti_n + Amp_Zin);
IN_diff = IN_p - IN_n;

% Clip to amplifier input range.
IN_diff = max(min(IN_diff, 0.5*Amp_InputRange), -0.5*Amp_InputRange);

% Calculate noise vector according to amplifier IRN, 
% scaling for the bandwidth vs sampling rate.
Vnoiseamp = NOISEon * Amp_IRN * randn(1,nop) * sqrt(Fsample/(2*Amp_BW) * 2/pi);

% Amplifier output with gain + added noise.
Vampout1 = Amp_Gain * (IN_diff + Vnoiseamp);

% Limit amplifier bandwidth with a simple one-pole IIR.
alpha1 = 1 - 2*pi*Amp_BW / Fsample;
Vampout = zeros(1,nop);
Vampout(1) = Vampout1(1);
for n = 2:nop
    Vampout(n) = alpha1 * Vampout(n-1) + (1 - alpha1) * Vampout1(n-1);
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% MODEL OF FILTER %%%
%%%%%%%%%%%%%%%%%%%%%%%

% Clip to filter input range.
Vfilterin = max(min(Vampout, 0.5*Filter_InputRange), -0.5*Filter_InputRange);

% Calculate noise vector for filter IRN.
Vnoisefilter = NOISEon * Filter_IRN * randn(1,nop) * sqrt(Fsample/(2*Filter_BW)*2/pi);

% Filter output with noise.
Vfilterout1 = Vfilterin + Vnoisefilter;

% Limit filter bandwidth with a simple one-pole IIR.
alpha2 = 1 - 2*pi*Filter_BW / Fsample;
Vfilterout = zeros(1,nop);
Vfilterout(1) = Vfilterout1(1);
for n = 2:nop
    Vfilterout(n) = alpha2 * Vfilterout(n-1) + (1 - alpha2) * Vfilterout1(n-1);
end

%%%%%%%%%%%%%%%%%%%%
%%% MODEL OF ADC %%%
%%%%%%%%%%%%%%%%%%%%

% Clip to ADC input range.
Vadcin = max(min(Vfilterout, 0.5*ADC_InputRange), -0.5*ADC_InputRange);

% Estimate ADC quantization noise from LSB size.
Vlsb = ADC_InputRange / 2^ADC_N;
Vnoisequantization = 0.289 * Vlsb; % ADC quantization noise (ideal ~ 1/sqrt(12)*LSB).

if (Vnoisequantization > ADC_IRN)
    disp('Error: total ADC noise < ADC quantization noise!');
end

% The remainder of ADC_IRN is random circuit noise on top of quantization.
Vnoiserandom = sqrt(ADC_IRN^2 - Vnoisequantization^2);
Vnoiseadc    = NOISEon * Vnoiserandom * randn(1,nop);

% ADC output: circuit noise + quantization + floor operation.
DCode    = floor((Vadcin + Vnoiseadc) / Vlsb);
MaxCode  = 2^(ADC_N - 1) - 1;
MinCode  = -2^(ADC_N - 1);
DigitalCode = max(min(DCode, MaxCode), MinCode);
k2 = round(1:(Fsample/ADC_fsample):nop);
DigitalCode = DigitalCode(k2);

if (sum(DCode>MaxCode)+sum(DCode<MinCode) > 0)
    disp('Your ADC is saturating');
else
    disp('Your signal fits in the ADC range');
end

plot(k2/Fsample, DigitalCode, 'b'); hold on;
plot(k/Fsample, IN_diff*Amp_Gain/ADC_InputRange * 2^ADC_N, 'g');
hold off;
xlabel('Time [s]');
ylabel('Digital output code');
legend('ADC output','Scaled input signal (reference)');
