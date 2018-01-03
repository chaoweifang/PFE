%Ncut_dir = 'E:\Data\PFE-final-results-20150920\KmSegNcutAffiNcut\ncut\dyn\dis\segimage';
%sc_MPB_dir = 'E:\Data\PFE-final-results-20150920\20150910_km_spectralclustering_mpb_16chs\seg\breg\dynkm\segimage';
%wsc_MPB_dir = 'E:\Data\PFE-final-results-20150920\20150910_km_respectralclustering_mpb_16chs\seg\breg\dynkm\segimage';
%wsc_NCut_dir = 'E:\Data\PFE-final-results-20150920\20150909_km_respectralclustering16chs\seg\breg\dynkm\segimage';


gpb_dir = 'E:\Data\PFE-final-results-20150920\UcmSegGPBSCG\ucm2\test\segimage';
%our_KM_MPB_dir='E:\Data\PFE-final-results-20150920\20150428_KmSegBregGMMAffiMPB\seg\breg\dynkm\segimage';
%our_KM_MPB1_dir='E:\Data\PFE-final-results-20150920\20150922_Km_Seg_BregGMM_AffiMPB\dyn\segimage';
MCG_ucm_dir='E:\Data\MCG-BSDS500-all-ucm\segimage';
our_MCG_ucm_dir='E:\Data\PFE-final-results-20150920\20150920_PFE-MCG-ucm\segimage';
%our_gpb1_dir = 'E:\Data\PFE-final-results-20150920\20150918_hi_seg_breg_step2_gmm_miu1500_r_100_10\ucm2_optweights1\segimage';
our_gpb2_dir = 'E:\Data\PFE-final-results-20150920\UcmSegBregGMMCombined_miu1000_r_100_10\ucm2\segimage';
%MCGvsPFE_dir = 'E:\Data\PFE-final-results-20150920\MCGvsPFE';
img_dir = 'E:\Data\BSDS\BSDS500\data\images\test';
inputfils = dir(fullfile(img_dir,'*.jpg'));

output_dir = 'E:\Data\PFE-final-results-20150920\VS_UCM';
if ~exist(output_dir,'dir'), mkdir(output_dir); end
for i=1:numel(inputfils)
    strimidx = inputfils(i).name(1:end-4);
    
    
    
%      save_path = fullfile(output_dir,[strimidx '_1Ncut.png']);
%     if ~exist(save_path,'file')
%         image = imread(fullfile(Ncut_dir,[strimidx '.bmp']));
%         imwrite(image,save_path);
%         %copyfile(fullfile(Ncut_dir,[strimidx '.png']),save_path);
%     end
%     
%     save_path = fullfile(output_dir,[strimidx '_2sc_MPB.png']);
%     if ~exist(save_path,'file')
%         copyfile(fullfile(sc_MPB_dir,[strimidx '.png']),save_path);
%     end
%     
%     save_path = fullfile(output_dir,[strimidx '_3wsc_MPB.png']);
%     if ~exist(save_path,'file')
%         copyfile(fullfile(wsc_MPB_dir,[strimidx '.png']),save_path);
%     end
%     
%     save_path = fullfile(output_dir,[strimidx '_4wsc_NCut.png']);
%     if ~exist(save_path,'file')
%         copyfile(fullfile(wsc_NCut_dir,[strimidx '.png']),save_path);
%     end
    
    
    save_path = fullfile(output_dir,[strimidx '_3gpb.png']);
    if ~exist(save_path,'file')
        image = imread(fullfile(gpb_dir,[strimidx '.bmp']));
        imwrite(image,save_path);
        %copyfile(fullfile(gpb_dir,[strimidx '.png']),save_path);
    end
    
    save_path = fullfile(output_dir,[strimidx '_4MCG_ucm.png']);
    if ~exist(save_path,'file'),
        copyfile(fullfile(MCG_ucm_dir,[strimidx '.png']),save_path);
    end
    
%      save_path = fullfile(output_dir,[strimidx '_5our_KM_MPB.png']);
%     if ~exist(save_path,'file'), 
%         image = imread(fullfile(our_KM_MPB_dir,[strimidx '.bmp']));
%         imwrite(image,save_path);
%         %copyfile(fullfile(our_KM_MPB_dir,[strimidx '.png']),save_path); 
%     end
%      save_path = fullfile(output_dir,[strimidx '_6our_KM_MPB1.png']);
%     if ~exist(save_path,'file'), 
%         %image = imread(fullfile(our_KM_MPB_dir,[strimidx '.bmp']));
%         %imwrite(image,save_path);
%         copyfile(fullfile(our_KM_MPB1_dir,[strimidx '.png']),save_path); 
%     end
     save_path = fullfile(output_dir,[strimidx '_6our_owt_ucm_la1000_r100_10.png']);
    if ~exist(save_path,'file')
        image = imread(fullfile(our_gpb2_dir,[strimidx '.bmp']));
        imwrite(image,save_path);
        %copyfile(fullfile(gpb_dir,[strimidx '.png']),save_path);
    end
    
    %save_path = fullfile(output_dir,[strimidx '_7our_gpb_la1500_r100_10.png']);
    %if ~exist(save_path,'file'), copyfile(fullfile(our_gpb1_dir,[strimidx '.png']),save_path); end
    
     save_path = fullfile(output_dir,[strimidx '_8our_MCG_ucm.png']);
    if ~exist(save_path,'file'), copyfile(fullfile(our_MCG_ucm_dir,[strimidx '.png']),save_path); end
   
end
