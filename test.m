tx_setup;
load('tx_os_ofdm_out.mat')
rx_in.signals.values = out.tx_out;
rx_in.time = out.tout;
S1 = reshape(out.tx_out,Nsc,Nsym);
%
% Rx signal at a reference point of the Rx array (without loss of generality, the first element of the ULA is assumed to be the reference);
% note that a seperate reference signal per target is generated, and the signals from all targets are summed up after adding the azimuth impact  
S1 = zeros(Nsc,Nsym,Num_Targets);
for nn=1:Num_Targets
    % radar channel parameters (delay and Doppler)
    tau_targ = 2.*Range(nn)./c0;
    gam = 2*(Velocity_kmh(nn)/3.6)/c0;
    f_Doppler = -2*(Velocity_kmh(nn)/3.6)/lambda;
    rec_power = 10^(0.05*Relative_power(nn));

    taun = tau_targ/T;
    D_Nc_t = diag_matr_gen(taun,Nsc);
    D_Nc_f = diag_matr_gen(f_Doppler*T/Nsc,Nsc);
    D_Nsym = diag_matr_gen(f_Doppler*Tt*alpha_pri,Nsym);
    S1(:,:,nn) = rec_power * exp(-1i*2*pi*Fc*tau_targ) * D_Nc_f * iF_Nc * conj(D_Nc_t) * tx_signal_ref * D_Nsym;
    if nn == Num_Targets % save the Rx power of the weakest target to set the noise level
        S1_ref = S1(:,:,nn);
    end
end

% add the impact from the angular (azimuth) domain and sum up signals from all targets
S1_ant = zeros(Nsc,Nsym,N_rx);
for nn=1:Num_Targets
    a_st = exp(1i*pi*cos(Azimuth(nn)*pi/180)*(0:N_rx-1).'); % Rx steering vector per target
    for nrx=1:N_rx
        S1_ant(:,:,nrx) = S1_ant(:,:,nrx) + a_st(nrx)*S1(:,:,nn);
    end
end