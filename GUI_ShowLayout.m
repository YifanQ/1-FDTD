function GUI_ShowLayout(pp, plot_h)

xx = pp.dimensionX;
yy = pp.dimensionY;
mat = zeros(xx, yy);

if pp.PMLFlag
    dd = pp.PMLthickness;

    xx = pp.dimensionX;
    yy = pp.dimensionY;
    mat(1:dd, :) = -1;
    mat(xx-dd+1:xx, :) = -1;
    mat(:, 1:dd) = -1;
    mat(:, yy-dd+1:yy) = -1;
end

if pp.PECScattFlag
    xx = pp.PECScattLocationX; lx = pp.PECScattDimensionX;
    yy = pp.PECScattLocationY; ly = pp.PECScattDimensionY;

    mat(xx:xx+lx-1, yy:yy+ly-1) = +1;
end

cla(plot_h,'reset');
plot_h.NextPlot = 'add';
xx = pp.dimensionX;
yy = pp.dimensionY;
imagesc(plot_h, 1:xx, 1:yy, mat, [-1, +1]);
axis(plot_h, 'image');
title('')

xx = pp.sourceLocationX;
yy = pp.sourceLocationY;
plot(plot_h, xx, yy, 'ro');
text(plot_h, xx+2, yy, 'source');

xx = pp.targetPointY;
yy = pp.targetPointX;
plot(plot_h, xx, yy, 'bx');
text(plot_h, xx+2, yy, 'observation point');

% https://www.mathworks.com/help/matlab/ref/matlab.graphics.axis.axes-properties.html#property_d119e54679
plot_h.Title.String = 'FDTD grid layout';
plot_h.XLabel.String = 'x / grid index';
plot_h.YLabel.String = 'y / grid index';

end
