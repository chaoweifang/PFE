function edg2ucm_all(database, gt_set,res_dir,paras)
if ~exist('database','var')
    database = 'VOCContext';
end
if ~exist('gt_set','var')
    gt_set = 'test';
end

% Adjust your paths in this file

% Results folder
%res_dir = fullfile('results',database,gt_set,'COB');
resucm_dir = res_dir;
if ~exist(resucm_dir,'dir')
    mkdir(resucm_dir);
end
edg_dir = fullfile(database_root_dir(database),'deepboundaries','mat');
thr=paras.thr;
fq=paras.fq;
% Which images to process
im_ids = database_ids(database,gt_set);
for ii=1:length(im_ids)
    %ii=2810;
    % Display evolution
    display(['Processing image ' num2str(ii) ' out of ' num2str(length(im_ids)) ' , name: ' im_ids{ii}]);
    
    % Read image
    %im = db_im(database, im_ids{ii});
    
    % Check if the result is already computed and readable
    % (glusterfs not very reliable)
    resucm_file = fullfile(resucm_dir,[im_ids{ii} '.mat']);
    try 
       ucm2=load(resucm_file);
       disp([resucm_file ' exists'])
    catch
        edg = loadvar(fullfile(edg_dir,[im_ids{ii} '.mat']));
        ucm2 = edg2ucm(edg,thr,fq);
        parsave_ucm(resucm_file,ucm2);
    end
end

end


function parsave_ucm(res_file,ucm2) %#ok<INUSD>
    save(res_file,'ucm2');
end
%function par


