%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 5XCC0 Assignment 1 - (C) Pieter Harpe %%%
%%% Only for use at TU/e %%%%%%%%%%%%%%%%%%%%
%%% Do not remove copyright %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Here you set all system-level specifications (to be found based on requirements)
Flow = 0.5; %Start of the bandwidth of the system [Hz]
Fhigh = 5e3; %End of the bandwidth of th-e system [Hz]
InputImpedance = 31.623e9; %Input impedance of the system (Z_in) [ohm] 
InputRefNoise = 10e-6; %Input-referred noise (V_IRN) [Vrms]
InputRange = 20e-3; %Input range of the system (V_inpp) [V]

%Here you set the impedance of the ETIs
Zeti_p = 100e3; %Impedance of the positive electrode [ohm]
Zeti_n = 100e3; %Impedance of the negative electrode [ohm]

%Here you decide which test signals to enable
APon = 1; %Set to 0/1 to disable/enable the AP test signal
CMon = 1; %Set to 0/1 to disable/enable the common-mode disturbance signal
DIFFon = 1; %Set to 0/1 to disable/enable the differential disturbance signal
NOISEon = 1; %Set to 0/1 to disable/enable the input-referred noise

%Here you define the parameters for the common-mode disturbance signal
CMfreq = 200; %Frequency of the common-mode power line signal [Hz]
CMamp = 50e-3; %Amplitude of the common-mode power line signal [Vpp]

%Here you define the parameters for the differential disturbance signal
DIFFfreq = 125; %Frequency of the differential disturbance signal [Hz]
DIFFamp = 50e-3; %Amplitude of the differential disturbance signal [Vpp]

%%% After this line you don't need to change anything,
%%% but please read the code if you like to understand what is going on

%Load exemplary AP recording and turn it into a differential signal
%Data record 1551.mat from: Benjamin Metcalfe, "Action potentials recorded from the L5 dorsal rootlet of rat using a multiple electrode array," Mendeley Data, Version 1, June 12, 2020, CC BY 4.0 License, doi: 10.17632/ybhwtngzmm.1
load ('1551.mat');
APtime = rawdata (:, 1)';
AP_p = 0.5 * rawdata (:, 2)' / 10000; %Compensate for 80dB gain in recorded data
AP_n = -AP_p;
nop = length (AP_p);
k = 1:nop;
Fsample = 1 / (APtime(2) - APtime (1));
TotalTime = nop / Fsample;

%Create a common-mode disturbance interference signal
CM_p = 0.5 * CMamp * sin (2 * pi * CMfreq * k / Fsample);
CM_n = CM_p;

%Create a differential disturbance signal
DIFF_p = 0.25 * DIFFamp * sin (2 * pi * DIFFfreq * k / Fsample);
DIFF_n = -DIFF_p;

%Add signals together, dependent on enable/disable setting
%These are the signals at the source (body) side
Signal_p = APon * AP_p + CMon * CM_p + DIFFon * DIFF_p;
Signal_n = APon * AP_n + CMon * CM_n + DIFFon * DIFF_n;

%Calculate resistive division due to electrode and circuit impedances
%These are the signals at the input of the system
Signal1 = Signal_p * InputImpedance / (Zeti_p + InputImpedance) - Signal_n * InputImpedance / (Zeti_n + InputImpedance);

%Calculate noise vector according to V_IRN, scale up to account for lower BW
Vnoise = NOISEon * InputRefNoise * randn (1, nop) * sqrt (Fsample / 2 / (Fhigh - Flow) * 2 / pi);

Signal2 = Signal1 + Vnoise;

%Use discrete time filter according to flow - fhigh bandwidth
alpha1 = 1 - 2 * pi * Flow / Fsample;
alpha2 = 1 - 2 * pi * Fhigh / Fsample;

Signal3 = zeros (1, nop);
Signal3 (1) = Signal2 (1);
Signal4 = zeros (1, nop);
Signal4 (1) = Signal3 (1);
for n = 2:nop
    Signal3 (n) = alpha1 * Signal3 (n - 1) + Signal2 (n) - Signal2 (n - 1);
    Signal4 (n) = alpha2 * Signal4 (n - 1) + (1 - alpha2) * Signal3 (n - 1);
end

plot (k / Fsample, Signal4, 'b');
hold on;
plot (k / Fsample, Signal1, 'g');
plot ([0 TotalTime], -[1 1] * InputRange / 2, 'r');
plot ([0 TotalTime], [1 1] * InputRange / 2, 'r');
hold off;
xlabel ('Time [s]');
ylabel ('Signal into system [V]');
titletext = sprintf ('Actual signal range: %.3fVpp, Allowed input range: %.3fVpp', max (Signal4) - min (Signal4), InputRange);
legend ('Signal with noise and limited BW', 'Signal without noise and without BW limitation', 'Lower limit input range', 'Upper limit input range');
title (titletext);