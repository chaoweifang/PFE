function ground_truth = get_ground_truth( database, image_id )
    gt = load(fullfile(database_root_dir(database), 'gt', [image_id '.mat']));
    if strcmp(database,'MSRC') || strcmp(database,'SBD')
        ground_truth{1}=gt.groundTruth{1}.Segmentation;
    elseif strcmp(database,'BSDS500')
        ground_truth={};
        for i=1:numel(gt.groundTruth)
            ground_truth=cat(1,ground_truth,gt.groundTruth{i}.Segmentation);
        end
    elseif strcmp(database,'VOCContext')
        ground_truth{1}=gt.LabelMap;
    else
        error(['Unknown database: ' database]);
    end
end

