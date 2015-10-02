%Demo for use of function: snr_eff_3dtse.m
%Kerry 2015

[SNReff,SNR,scn_time,sig_loss_T1,mtf_xy,mtf_z,FA] = snr_eff_3dtse([300;200;80],[1.5;1.5;1.6],3.5,80,1010,75);

plotcell(FA); title('FA'); 
plot(mtf_xy); title('MTF xy');
plot(mtf_z);  title('MTF z');


figure; plot([1:lengt(SNReff)]*2, SNReff);