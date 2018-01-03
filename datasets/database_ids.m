function ids = database_ids( database, gt_set )
    index_file = fullfile(root_dir,'datasets','gt_sets',database, [gt_set '.txt']);
    fileID = fopen(index_file);
    ids = textscan(fileID, '%s');
    ids = ids{1};
    fclose(fileID);
end

