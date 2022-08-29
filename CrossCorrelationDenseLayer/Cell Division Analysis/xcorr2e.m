function I = xcorr2e(template, im, shape)
% helper function to modify matlabs xcorr2 function to use padding
  if (nargin == 2) || strcmp(shape,'full')
      I = xcorr2(template, im);
      return
  end
  switch shape
      case 'same'
          pad = floor(size(template)./2);
          center = size(im);
      case 'valid'
          pad = size(template) - 1;
          center = size(im) - pad;
      otherwise
          throw(MException('normxcorr2e:BadInput',...
              'SHAPE must be ''full'', ''same'', or ''valid''.'));
  end
  I = xcorr2(gpuArray(template), gpuArray(im));
  I = I([false(1,pad(1)) true(1,center(1))], ...
        [false(1,pad(2)) true(1,center(2))]);

end

