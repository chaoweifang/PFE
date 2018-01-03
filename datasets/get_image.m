function image = get_image( database, image_id )
    if strcmp(database,'MSRC') 
        image = imread(fullfile(database_root_dir(database), 'images', [image_id '.bmp']));
    elseif strcmp(database,'SBD') || strcmp(database,'SeBD') || ...
            strcmp(database,'BSDS500') || strcmp(database,'VOCContext')
        image = imread(fullfile(database_root_dir(database), 'images', [image_id '.jpg']));
    else
        error(['Unknown database: ' database]);
    end
end

