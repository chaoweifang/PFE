addpath('./datasets/')
addpath('./ExtraImages/')
database='MSRC';
ids=database_ids(database,'test');
db_root=database_root_dir(database);
seg_root=fullfile(db_root,'gtimg');
if ~exist(seg_root,'dir')
    mkdir(seg_root);
end
for i=1:numel(ids)
    im=ids{i};
    img=get_image(database,im);
    gt = get_ground_truth(database,im);
    for j=1:numel(gt)
        imgseg = showSegResults(img,double(gt{j}));
        fil = fullfile(seg_root,[im '_' num2str(j,'%02d') '.png']);
        imwrite(imgseg,fil);
    end
end
    