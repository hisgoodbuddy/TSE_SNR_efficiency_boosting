%Demo for use of function: snr_eff_3dtse.m
%Kerry 2015

%Optimal TSE vs Plateau Ratio
TSEes = 3:0.2:8;
Ratio = 0.1:0.05:0.9; 

for n_es=1:length(TSEes)
    disp(['TSEes:' num2str(TSEes(n_es))]);
    for n_r=1:length(Ratio)
        if(n_es==1)&(n_r==1)
        else
            disp(['...Ratio:' num2str(Ratio(n_r))]);
            [tsefactor, SNReff,SNR,scn_time,sig_loss_T1,mtf_xy,mtf_z,FA] = ...
                snr_eff_3dtse([150;60;51],[0.6;0.6;3],TSEes(n_es),120,Ratio(n_r),1000,50);

            Opti_SNReff(n_es,n_r)=max(SNReff);
            Opti_TSE(n_es,n_r)=tsefactor(find(SNReff==max(SNReff)));  

            save (['SNReff_results_es_N',num2str(n_es),'_ratio_N',num2str(n_r),'.mat'], 'tsefactor', 'SNReff','SNR','scn_time','sig_loss_T1','mtf_xy','mtf_z','FA');
            clear tsefactor SNReff SNR scn_time sig_loss_T1 mtf_xy mtf_z FA;
            close all;
        end
    end
end


