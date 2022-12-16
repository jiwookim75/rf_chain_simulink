%Running receiver
%Hessam Mohajeri
%12/5/2022
%First run the transmitter and the data is captured in 
%also run the matlab for targets to get the signals for targes if needed
%Then decide as to which data you want to process through the receiver.
%stucture 'out' is the transmitter output
%structure 'target_out is the receiver targets
%This allows to fine tune the transmitter wihtout the receiver and speed
%the simulations.
%once that is run then run the receiver in the file 'rx_ovs_ofdm'
%
tx_setup;
load('tx_os_ofdm_out.mat')
rx_in.signals.values = out.tx_out;
rx_in.time = out.tout;
%sim('rx_ovs_ofdm.slx',1.3e-3); %or run the file itself
%after running the receiver then push the data through the post processing.