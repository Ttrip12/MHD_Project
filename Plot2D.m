clear all
close all

axes = readmatrix('xy.csv');
FilesU = dir('DataU');
num_files = length(FilesU);
cd DataU
dataU{:} = zeros(num_files-2);
for i = 3:num_files
   fileU = FilesU(i).name;
   dataU{i-2} = readmatrix(fileU);
end
%num_files = length(FilesU);
cd ..
FilesV = dir('DataV');
cd DataV
dataV{:} = zeros(num_files-2);

for i = 3:num_files
   fileV = FilesV(i).name;
   dataV{i-2} = readmatrix(fileV);
end
cd ..

FilesMag = dir('DataMagnitude');
cd DataMagnitude
dataMag{:} = zeros(num_files-2);

for i = 3:num_files
   fileMag = FilesMag(i).name;
   dataMag{i-2} = readmatrix(fileMag);
end
cd ..
%% Plots
cd Figures
sc = 8;

for i=1:num_files - 2
    f = figure(i);
    plotname = "2D Sine Nu = 0.05 Time = " + 4*(i - 1) + " s";
    fname = 4*(i - 1) + "_2D_Sine_different_n_0_05";
    h = pcolor(axes(:,1),axes(:,2),dataMag{i});
    set(h, 'EdgeColor', 'none');
    grid off
    hold on
    quiver(axes(1:sc:end,1),axes(1:sc:end,2),dataV{i}(1:sc:end,1:sc:end),dataU{i}(1:sc:end,1:sc:end),'Color','k');
    hold off
    %s.EdgeColor = 'none';
    xlim([0 2*pi])
    xticks([0 pi 2*pi]);
    xticklabels({'0','\pi','2\pi'});
    ylim([0 2*pi])
    yticks([0 pi 2*pi]);
    yticklabels({'0','\pi','2\pi'});
    c = colorbar;
    c.Label.String = 'Velocity Magnitude)';
    title(plotname)
    xlabel("X")
    ylabel("Y")
    saveas(f,fname,'jpg')

end

cd ..