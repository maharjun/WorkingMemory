function [ fig, res ] = DisplayImage(ImageMat, varargin)
%DISPLAYIMAGE 
% 
% fig = DisplayImage(ImageMat)
% fig = DisplayImage(ImageMat, 'Option', OptionValue ...)
% 
% Displays the image

	withCbar = true;
	withAxes = true;
	ActualSize = false;
	Tag = [];
	nParams = floor(length(varargin)/2);
	for i = 1:nParams
		if strcmpi(varargin{2*i-1}, 'axes')
			withAxes = varargin{2*i};
		elseif strcmpi(varargin{2*i-1}, 'cbar')
			withCbar = varargin{2*i};
		elseif strcmpi(varargin{2*i-1}, 'actualsize')
			ActualSize = varargin{2*i};
		elseif strcmpi(varargin{2*i-1}, 'tag')
			Tag = varargin{2*i};
		else
			break;
		end
	end
	
	fig = figure;
	ImgHandle = imshow(ImageMat);
	
	if ~isempty(Tag)
		set(ImgHandle, 'Tag', Tag);
	end
	
	if ndims(ImageMat) == 2 && withCbar
		colorbar;
	end
	h = gca();
	set(h, 'Units', 'pixels');
	if withAxes
		set(h, 'Visible', 'on');
	end
	s = get(h, 'Position');
	s = floor(size(ImageMat, 2)*96/s(3))/96;
	
	% 65060241
	if ActualSize
		set(fig, 'Units', 'pixels');
		set(fig, 'Position', [0, 50, size(ImageMat, 2) + 125, size(ImageMat, 1) + 100]);
		set(h, 'Position', [50 50 size(ImageMat, 2) size(ImageMat, 1)]);
		res = 1;
	else
		set(fig, 'Units', 'pixels');
		set(fig, 'Position', [0, 50, size(ImageMat, 2)/s + 125, size(ImageMat, 1)/s + 100]);
		set(h, 'Position', [52 52 size(ImageMat, 2)/s size(ImageMat, 1)/s]);
		res = s;
	end
end

