close all; clear; clc;

xdim = 500;
ydim = 500;
pane = ones(xdim, ydim);

BoxLeftSide = zeros(5,1);
BoxTop = zeros(1,6);
BoxBottom = zeros(1,8);
TopStep = ones(9,9);
TopDiag = ones(11,4);

for i = 1:1:9
    if (i > 1)
        TopStep(i,i) = 0;
        TopStep(i-1,i) = 0;
    else
        TopStep(i,i) = 0;
    end
end
BottomStep = TopStep;

for j = 1:1:4
    if (j > 1)
        TopDiag(j:j+3,j) = 0;
    else
        TopDiag(j,j) = 0;
        TopDiag(j+1,j) = 0;
    end
end
BottomDiag = TopDiag;

pane(xdim-400:1:(xdim-400)+4,ydim-400) = BoxLeftSide;
pane((xdim-400)-1,(ydim-400):(ydim-400)+7) = BoxBottom;
pane((xdim-400)+5,(ydim-400):(ydim-400)+5) = BoxTop;
pane((xdim-400)+6:(xdim-400)+14,(ydim-400)+5:(ydim-400)+13) = TopStep;
pane((xdim-400):(xdim-400)+8,(ydim-400)+7:(ydim-400)+15) = pane((xdim-400):(xdim-400)+8,(ydim-400)+7:(ydim-400)+15).*BottomStep;
%pane((xdim-400)+15:(xdim-400)+25,(ydim-400)+13:(ydim-400)+16) = TopDiag;


colormap('jet');
pcolor(pane);
shading interp
colorbar