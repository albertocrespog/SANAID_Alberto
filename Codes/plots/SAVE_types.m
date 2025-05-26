function SAVE_types(fname,name,gca,gcf)
saveas(gca,fullfile(fname, name),'jpeg');
saveas(gcf,fullfile(fname, name),'fig');
saveas(gcf,fullfile(fname, name),'pdf');
saveas(gcf,fullfile(fname, name),'epsc');
% saveas(gcf,fullfile(fname, name),'bmp');
saveas(gcf,fullfile(fname, name),'png');
