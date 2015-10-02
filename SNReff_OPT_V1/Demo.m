%Demo for use of function: snr_eff_3dtse.m
%Kerry 2015

%Optimal TSE vs Plateau Ratio

Ratio = 0.1:0.05:0.9; 

for n_r=1:length(Ratio)
    disp(['Ratio:' num2str(Ratio(n_r))]);
    [tsefactor, SNReff,SNR,scn_time,sig_loss_T1,mtf_xy,mtf_z,FA] = ...
        snr_eff_3dtse([300;200;80],[1.5;1.5;1.6],3.5,150,Ratio(n_r),1010,75);
    
    Opti_SNReff(n_r)=max(SNReff);
    Opti_TSE(n_r)=tsefactor(find(SNReff==max(SNReff)));  
    
    save (['SNReff_results_ratio_N',num2str(n_r),'.mat'], tsefactor, SNReff,SNR,scn_time,sig_loss_T1,mtf_xy,mtf_z,FA);
    clear tsefactor SNReff SNR scn_time sig_loss_T1 mtf_xy mtf_z FA;
end


