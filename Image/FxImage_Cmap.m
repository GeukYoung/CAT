function [Cmap] = FxImage_Cmap(option)
    switch option
        case 1
            load Cmap1;
            Cmap = Cmap1;
        case 2
            load Cmap2;
            Cmap = Cmap2;
        case 3
            load Cmap3;
            Cmap = Cmap3;
        case 4
            load Cmap4;
            Cmap = Cmap4;
    end
end