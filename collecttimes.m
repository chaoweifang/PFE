
root = '/home/cwfang/Results/BSDS500/DBswPFEp0.8_journal/ucm/re_multi_fq11_thr0.35_satsPb0.8_or3';
fils=dir(fullfile(root,'*.mat'));
t=zeros(1,4);
%t_ini=0;
%t_pfe=0;
%t_others=0;
for i=1:length(fils)
    aa=load(fullfile(root,fils(i).name));
    t=t+aa.time;
end
t/length(fils)