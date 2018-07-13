
%I = imread('~/MATLAB/BSDS500/data/images/test/100007.jpg');
function im2ucmpnorm_all(database, gt_set,res_dir,param_s,param_ucm,param_ini)

resucm_dir = fullfile(res_dir, 'ucm');
if ~exist(resucm_dir,'dir')
    mkdir(resucm_dir);
end
resvec_dir = fullfile(res_dir,'vec');
if ~exist(resvec_dir,'dir')
    mkdir(resvec_dir);
end

resYor_dir = fullfile(database_root_dir(database),[param_ini.flag 'ini']);
if ~exist(resYor_dir,'dir')
    mkdir(resYor_dir);
end

edg_dir = fullfile(database_root_dir(database),'deepboundaries','mat');
im_ids = database_ids(database,gt_set);

for ii=1:length(im_ids)
    % Read image
    %tic
    strimidx =  im_ids{ii};%im_names(ii).name(1:end-4);
    disp(['Processing image ' strimidx ' idx ' num2str(ii)])
    %image = imread(fullfile(imgpath,[strimidx '.jpg']));
    im = get_image(database, strimidx);
    % Check if the result is already computed
    ucmfile = fullfile(resucm_dir,[strimidx '.mat']);
    vecfile = fullfile(resvec_dir,[strimidx '.mat']);
    Yorfile = fullfile(resYor_dir,[strimidx '.mat']);
    %try
    %    ucm=load(ucmfile);
    %    disp(['ucm of ' strimidx ' exists!']);
    %catch
        %vec_exist_flag=true;
        %try
        %    vecs = load(vecfile);
        %    vecs = vecs.vecs;
        %catch
            vecs = [];
            %vec_exist_flag=false;
        %end
        Yor_exist = true;
        try 
            Yor = load(Yorfile);
            Yor=Yor.Yor;
        catch
            Yor=[];
            Yor_exist=false;
        end
        %edge_hed = readHedData(hedpath,strimidx);
        edg = load(fullfile(edg_dir,[strimidx '.mat']));
        edg = edg.edg;
        if (size(edg,1)~=size(im,1)) || (size(edg,2)~=size(im,2))
            error('wrong input edge data!');
        end
        % Call the actual code[ucm2,vecs]
        [ucm2,~,vecs,Yor,t]= img2ucm(im,vecs,Yor,edg,param_s,param_ucm,param_ini);    
        %time=toc;
        disp(['Finished image ' strimidx ' idx ' num2str(ii) ' time:' num2str(t(1))])
        % Store ucms at each scale separately
        parsaveucm2(ucmfile,ucm2,t);
        %if ~vec_exist_flag
        parsavevecs(vecfile,vecs);
        %end
        if ~Yor_exist
            parsaveYor(Yorfile,Yor);
        end
    %end
end
end

function parsaveucm2(res_file,ucm2,time) %#ok<INUSD>
    save(res_file,'ucm2','time');
end
function parsavevecs(res_file,vecs) %#ok<INUSD>
    save(res_file,'vecs');
end
function parsaveYor(res_file,Yor) %#ok<INUSD>
    save(res_file,'Yor');
end