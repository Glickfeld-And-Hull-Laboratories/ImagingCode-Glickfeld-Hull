function [bin_img ] = get_frame(image_dest, frame, info,ROI_x, ROI_y, BIN_SIZE)
c_img =  double(imread(image_dest, frame, 'info', info)); %old format Mati wrote
%c_img =  double(readtiff(image_dest));
if(exist('ROI_x' , 'var') && exist('ROI_y' , 'var') &&   ~isempty(ROI_x) && ~isempty(ROI_y))
    c_img = c_img(ROI_x, ROI_y);
end

% --- cut edge
n_sz = size(c_img) - mod( size(c_img), BIN_SIZE);
c_img = c_img(1:n_sz(1), 1:n_sz(2));
% ----- bin data
b_sz =  n_sz/BIN_SIZE;
bin_img=reshape(mean(mean(reshape(c_img,[BIN_SIZE b_sz(1) BIN_SIZE b_sz(2) ]),1),3),[b_sz(1) b_sz(2)]);

