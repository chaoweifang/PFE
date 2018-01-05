function [impad] = padReflect(im,r)
% function [impad] = padReflect(im,r)
%
% Pad an image with a border of size r, and reflect the image into
% the border.
%
% David R. Martin <dmartin@eecs.berkeley.edu>
% March 2003

impad = zeros(size(im)+2*r);
impad(r+1:end-r,r+1:end-r) = im; % middle
impad(1:r,r+1:end-r) = flipud(im(1:r,:)); % top
impad(end-r+1:end,r+1:end-r) = flipud(im(end-r+1:end,:)); % bottom
impad(r+1:end-r,1:r) = fliplr(im(:,1:r)); % left
impad(r+1:end-r,end-r+1:end) = fliplr(im(:,end-r+1:end)); % right
impad(1:r,1:r) = rot90(im(1:r,1:r),2); % top-left
impad(1:r,end-r+1:end) = rot90(im(1:r,end-r+1:end),2); % top-right
impad(end-r+1:end,1:r) = rot90(im(end-r+1:end,1:r),2); % bottom-left
impad(end-r+1:end,end-r+1:end) = rot90(im(end-r+1:end,end-r+1:end),2); % bottom-right
