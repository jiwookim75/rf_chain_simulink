Bringin up the simulator and various files
1/6/2022
Hessam Mohajeri
- rf_system simulink can run sinewaves and check end to end design. It is quick and straight forward.
- tx_os_ofdm use this file to run the data stream into the transmitter and capture the output. You need over 8ms to run this. Long run but if all is okay you can run
 it once and be done. The output is captured in out. and can be stored manually to be used later. the name of the file is tx_ovs_ofdm_out.mat
-tx_setup.m is used to setup the transmitter. Run this file first to get the initial settings up running.
- up_converter.slx and dn_converter.slx are the two subcircuits used for simulations. Make changes to the architecture here to try different designs.