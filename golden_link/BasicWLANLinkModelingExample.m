%% Basic WLAN Link Modeling
%
% This example shows how to create a basic WLAN link model using 
% WLAN Toolbox(TM). An IEEE(R) 802.11(TM) [ <#29 1> ] VHT packet
% is created, passed through a TGac channel. The received signal is 
% equalized and decoded in order to recover the transmitted bits.

% Copyright 2015-2021 The MathWorks, Inc.

%% Introduction
% This example shows how a simple transmitter-channel-receiver simulation
% may be created using functions from WLAN Toolbox. A VHT transmit and
% receive link is implemented as shown in the figure below. A VHT packet is
% transmitted through a TGac channel, demodulated and the equalized symbols
% are recovered. The equalized symbols are decoded to recover the
% transmitted bits.
%
% <<../BasicWLANlinkModelingDiagram.png>>

%% Waveform Generation
% An 802.11ac VHT transmission is simulated in this example. The transmit
% parameters for the VHT format of the 802.11(TM) standard are configured
% using a VHT configuration object. The <docid:wlan_ref#buw6eev
% wlanVHTConfig> creates a VHT configuration object. In this example the
% object is configured for a 20 MHz channel bandwidth, MCS 5 and single
% transmit antenna.

% Create a format configuration object for a SISO VHT transmission
cfgVHT = wlanVHTConfig;
cfgVHT.NumTransmitAntennas = 1;    % Transmit antennas
cfgVHT.NumSpaceTimeStreams = 1;    % Space-time streams
cfgVHT.APEPLength = 4096;          % APEP length in bytes
cfgVHT.MCS = 4;                    % Single spatial stream, 16-QAM
cfgVHT.ChannelBandwidth = 'CBW20'; % Transmitted signal bandwidth
Rs = wlanSampleRate(cfgVHT);       % Sampling rate

%%
% A single VHT packet is generated consisting of training, signal and data
% fields:
%
% * Non-HT Short Training Field (L-STF)
% * Non-HT Long Training Field (L-LTF)
% * Non-HT Signal (L-SIG) field
% * VHT Signal A (VHT-SIG-A) field
% * VHT Short Training Field (VHT-STF)
% * VHT Long Training Field (VHT-LTF)
% * VHT Signal B (VHT-SIG-B) field
% * Data field
%
% These fields are generated separately using functions from WLAN Toolbox
% and are concatenated to produce a VHT transmit packet.

%%
% The first field in the PPDU is the L-STF and is used for the start of
% packet detection and automatic gain control (AGC) setting. It is also
% used for initial frequency offset estimation and coarse timing
% synchronization.  The <docid:wlan_ref#buzp8qk-1 wlanLSTF> function
% generates the L-STF field in the time-domain using some of the parameters
% included in configuration object |cfgVHT|.
lstf = wlanLSTF(cfgVHT);  

%%
% The L-LTF is used for fine time synchronization, channel estimation and
% fine frequency offset estimation. The <docid:wlan_ref#buzia15 wlanLLTF>
% function generates the L-LTF in the time-domain.
lltf = wlanLLTF(cfgVHT);  

%%
% The L-SIG field carries packet configuration such as data rate,
% modulation and code rate for non-HT format. The <docid:wlan_ref#buzdcd9
% wlanLSIG> function generates the L-SIG field in the time-domain.
lsig = wlanLSIG(cfgVHT);

%%
% The figure below shows the L-STF, L-LTF and L-SIG fields. These fields
% are common to the VHT, HT-Mixed and non-HT OFDM transmission formats.
nonHTfield = [lstf;lltf;lsig]; % Combine the non-HT preamble fields

%%
%
% <<../BasicWLANlinkModelingNonHTpreamble.png>>
%

%%
% The VHT specific signal and training fields are generated after the
% non-HT preamble fields. The purpose of the VHT-SIG-A field is to provide
% information to allow the receiver to decode the data payload. The
% VHT-SIG-A is composed of two symbols VHT-SIG-A1 and VHT-SIG-A2. The
% <docid:wlan_ref#buz0vbn-1 wlanVHTSIGA> function generates the VHT-SIG-A
% field in the time-domain.
vhtsiga = wlanVHTSIGA(cfgVHT);

%%
% The purpose of the VHT-STF is to improve the gain control estimation in a
% MIMO transmission and help the receiver detect the repeating pattern
% similar to the L-STF field. The <docid:wlan_ref#buzzfyy wlanVHTSTF>
% function generates the VHT-STF field in the time-domain.
vhtstf = wlanVHTSTF(cfgVHT);

%%
% The VHT-LTF provides a mean for the receiver to estimate the channel
% between the transmitter and the receiver. Depending on the number of
% space time streams, it consists of 1,2,4,6 or 8 VHT-LTF symbols. The
% <docid:wlan_ref#buz0n9g-1 wlanVHTLTF> function generates the VHT-LTF in
% the time-domain.
vhtltf = wlanVHTLTF(cfgVHT);

%%
% The VHT-SIG-B field is used to set the data rate and the length of the
% data field payload of the transmitted packet. The
% <docid:wlan_ref#buz01qz-1 wlanVHTSIGB> function generates the VHT-SIG-B
% field in the time-domain.
vhtsigb = wlanVHTSIGB(cfgVHT);

%%
% Construct the preamble with the generated signal and training fields for
% the VHT format.
preamble = [lstf;lltf;lsig;vhtsiga;vhtstf;vhtltf;vhtsigb];

%%
% The <docid:wlan_ref#buz05r4-1 wlanVHTData> function generates the
% time-domain VHT data field. The VHT format configuration |cfgVHT|
% specifies the parameters for generating the data field from the PSDU
% bits. The |cfgVHT.PSDULength| property gives the number of bytes to be
% transmitted in the VHT data field. This property is used to generate the
% random PSDU bits |txPSDU|.

rng('default') % Initialize the random number generator
txPSDU = randi([0 1],cfgVHT.PSDULength*8,1); % Generate PSDU data in bits
data = wlanVHTData(txPSDU,cfgVHT);

% A VHT waveform is constructed by prepending the non-HT and VHT
% preamble fields with data
txWaveform = [preamble;data]; % Transmit VHT PPDU

% numPkts = 2;
% bits = [1;0;0;1];
% scramInit = randi([1 127],numPkts,1);
% txWaveform = wlanWaveformGenerator(bits,cfg,'NumPackets',numPkts,'IdleTime',30e-6,'ScramblerInitialization',scramInit);


%%
% Alternatively the waveform for a given format configuration can also be
% generated using a single function call
% <docid:wlan_ref#buwv1f3 wlanWaveformGenerator> function.
% This function can produce one or more VHT packets. By default OFDM
% windowing is applied to the generated waveform. For more information on
% OFDM windowing, see the reference page for the
% <docid:wlan_ref#buwv1f3 wlanWaveformGenerator>
% function.

%% Channel Impairments
% This section simulates the effects of over-the-air transmission. The
% transmitted signal is impaired by the channel and AWGN. The level of the
% AWGN is given in dBs. In this example the TGac channel model [ <#29 2> ]
% is used with delay profile Model-B. For this delay profile when the
% distance between transmitter and receiver is greater than or equal to 5
% meters, the model is in Non-Line-of-Sight (N-LOS) configuration. This is
% described further in the help for <docid:wlan_ref#buy2ffw-1
% wlanTGacChannel>.

% Parameterize the channel
tgacChannel = wlanTGacChannel;
tgacChannel.DelayProfile = 'Model-B';
tgacChannel.NumTransmitAntennas = cfgVHT.NumTransmitAntennas;
tgacChannel.NumReceiveAntennas = 1;
tgacChannel.LargeScaleFadingEffect = 'None';
tgacChannel.ChannelBandwidth = 'CBW20';
tgacChannel.TransmitReceiveDistance = 5;
tgacChannel.SampleRate = Rs;
tgacChannel.RandomStream = 'mt19937ar with seed';
tgacChannel.Seed = 10;

% Pass signal through the channel. Append zeroes to compensate for channel
% filter delay
txWaveform = [txWaveform;zeros(10,1)];
chanOut = tgacChannel(txWaveform);

snr = 40; % In dBs
rxWaveform = awgn(chanOut,snr,0);

% Display the spectrum of the transmitted and received signals. The
% received signal spectrum is affected by the channel
spectrumScope  = spectrumAnalyzer(SampleRate=Rs, ...            
            AveragingMethod='exponential',ForgettingFactor=0.99, ...
            YLimits=[-30 10],ShowLegend=true, ... 
            ChannelNames={'Transmitted waveform','Received waveform'});
spectrumScope([txWaveform rxWaveform]);

%% Channel Estimation and Equalization
% In this section the time-domain VHT-LTF is extracted from the received
% waveform. The waveform is assumed to be synchronized to the start of the
% packet by taking the channel filter delay into account. The VHT-LTF is
% demodulated and is used to estimate the channel. The received signal is
% then equalized using the channel estimate obtained from the VHT-LTF.

%%
% In this example the received signal is synchronized to the start of the
% packet by compensating for a known channel filter delay. For more
% information on how to automatically detect and synchronize to the
% received signal see the following examples:
%
% * <docid:wlan_ug#example-HTMIMOPacketErrorRateExample
% 802.11n Packet Error Rate Simulation for 2x2 TGn Channel>
% * <docid:wlan_ug#example-VHTMIMOPacketErrorRateExample
% 802.11ac Packet Error Rate Simulation for 8x8 TGac Channel>

% %% Real packet import
% int_i_1 = int_i(905052:914718);
% int_q_1 = int_q(905052:914718);
% 
% rxWaveform = complex(int_i_1,int_q_1);


%% 
% Packet detect and determine coarse packet offset

% rxWaveform = complex(rx_w1(400713:412189,1), rx_w1(400713:412189,2));

% threshold = 0.75;
% [coarsePktOffset,M] = wlanPacketDetect(txWaveform,cfgVHT.ChannelBandwidth,0,1);
% thres_idx = find(round(M,2) > threshold);
% thres_idx = thres_idx(1);
% rxWaveform = rxWaveform(thres_idx+1:end);
% 
% spectrumScope  = spectrumAnalyzer(SampleRate=Rs, ...            
%             AveragingMethod='exponential',ForgettingFactor=0.99, ...
%             YLimits=[-30 10],ShowLegend=true, ... 
%             ChannelNames={'Transmitted waveform','Received waveform'});
% spectrumScope(rxWaveform);

% analyzer = WaveformAnalyzer;
% process(analyzer,rxWaveform,cfgVHT.ChannelBandwidth,20000000);
% Display a summary of the detected packets.

% detectionSummary(analyzer);

% packet search from simulator
% rx_signal_detect = rx_find_packet_edge([zeros(20,1);rxWaveform],threshold);


chInfo = info(tgacChannel); % Get characteristic information
% Channel filter delay, measured in samples 
chDelay  = chInfo.ChannelFilterDelay;
rxWaveform = rxWaveform(chDelay+1:end,:);

% wlan waveform generator
% load('test_signals\wlan_gen_20_MHz_4_MCS_1_packet_without_noise.mat');
load('test_signals\wlan_gen_20_MHz_3_MCS_2_packet_without_noise.mat');
rxWaveform = [waveStruct.waveform];
cfgVHT = waveStruct.config.waveform;
%% Detect packets

threshold = 0.75;
init_value = 1;
[coarsePktOffset,M] = wlanPacketDetect(rxWaveform,cfgVHT.ChannelBandwidth,0,1);

thres_idx_array = zeros(2,1);

for i = 1:2
    thres_idx = find(M(init_value:end) > threshold);
    thres_idx_array(i) = thres_idx(1) + init_value - 1;
    init_value = init_value + 10000;
end

for i = 1:length(thres_idx_array)
    if (i == length(thres_idx_array))
        rxWaveform = rxWaveform(thres_idx_array(i)+1:end);
    else
        rxWaveform = rxWaveform(thres_idx_array(i):thres_idx_array(i+1));
    end

%%
% After synchronization the receiver has to extract the relevant fields
% from the received packet. The <docid:wlan_ref#bu1k_fl-1
% wlanFieldIndices> function is used to return the start and end
% time-domain sample indices of all fields relative to the first sample in
% a packet. These indices are used to extract the required fields for
% further processing.
indField = wlanFieldIndices(cfgVHT);

%%
% An estimate of the noise power after OFDM demodulation is required to
% perform MMSE equalization on the received OFDM symbols. In this example
% the noise power in the VHT fields is estimated using the demodulated
% L-LTF symbols.The L-LTF is extracted from the received waveform and is
% demodulated using the <docid:wlan_ref#buz04s0
% wlanLLTFDemodulate> function.
indLLTF = indField.LLTF(1):indField.LLTF(2);
demodLLTF = wlanLLTFDemodulate(rxWaveform(indLLTF),cfgVHT);
% Estimate noise power in VHT fields
nVar = helperNoiseEstimate(demodLLTF,cfgVHT.ChannelBandwidth,cfgVHT.NumSpaceTimeStreams);

%%
% To extract the VHT-LTF from the received signal the start and end indices
% are used to generate a vector of indices.
indVHTLTF = indField.VHTLTF(1):indField.VHTLTF(2);

%%
% The VHT-LTF is used to estimate the channel between all space-time
% streams and receive antennas. The VHT-LTF is extracted from the received
% waveform and is demodulated using the <docid:wlan_ref#bux3cou-1
% wlanVHTLTFDemodulate> function.
demodVHTLTF = wlanVHTLTFDemodulate(rxWaveform(indVHTLTF,:),cfgVHT);

% Plot equalized symbols
% constellationDiagram = comm.ConstellationDiagram;
% constellationDiagram.ReferenceConstellation = wlanReferenceSymbols(cfgVHT);
% Compare received and reference constellation  
% constellationDiagram(reshape(demodVHTLTF,[],1));      
% constellationDiagram.Title = 'Equalized Data Symbols';
%%
% The channel estimate includes the effect of the applied spatial mapping
% and cyclic shifts at the transmitter for a multi antenna configuration.
% The <docid:wlan_ref#buy3ynv wlanVHTLTFChannelEstimate>
% function returns the estimated channel between all space-time streams and
% receive antennas.
chanEstVHTLTF = wlanVHTLTFChannelEstimate(demodVHTLTF,cfgVHT);

%%
% The transmit signal encounters a deep fade as shown in the channel
% frequency response in the figure below. The effect of channel fades can
% also be seen in the spectrum plot shown previously.
figure
plot(20*log10(abs(chanEstVHTLTF)));
grid on;
title('Estimated Channel Response');
xlabel('Subcarrier index');
ylabel('Power (dB)');

%%
% To extract the data field from the received signal the start and end
% indices for the data field are used to generate a vector of indices.
indData = indField.VHTData(1):indField.VHTData(2);

% Recover the bits and equalized symbols in the VHT Data field using the
% channel estimates from VHT-LTF
[rxPSDU,~,eqSym] = wlanVHTDataRecover(rxWaveform(indData,:),chanEstVHTLTF,nVar,cfgVHT);
        
% Compare transmit and receive PSDU bits       
% numErr = biterr(txPSDU,rxPSDU);     

%% 
% The following plot shows the constellation of the equalized symbols at
% the output of the <docid:wlan_ref#buxdph7 wlanVHTDataRecover>
% function compared against the reference constellation. Increasing the
% channel noise should begin to spread the distinct constellation points.

% Plot equalized symbols
constellationDiagram = comm.ConstellationDiagram;
constellationDiagram.ReferenceConstellation = wlanReferenceSymbols(cfgVHT);
% Compare received and reference constellation  
constellationDiagram(reshape(eqSym,[],1));      
constellationDiagram.Title = 'Equalized Data Symbols';
end
%% Appendix
% This example uses this helper function.
%
% * <matlab:edit('helperNoiseEstimate.m') helperNoiseEstimate.m>

%% Selected Bibliography
% # IEEE Std 802.11(TM)-2020. IEEE Standard for Information Technology
% - Telecommunications and Information Exchange between Systems - Local and
% Metropolitan Area Networks - Specific Requirements - Part 11: Wireless
% LAN Medium Access Control (MAC) and Physical Layer (PHY) Specifications.
% # Breit, G., H. Sampath, S. Vermani, et al. TGac Channel Model Addendum.
% Version 12. IEEE 802.11-09/0308r12, March 2010.
