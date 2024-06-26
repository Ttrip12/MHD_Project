clear all
close all

axes = readmatrix('xy.csv');
Files = dir('DataMagnitude');
num_files = length(Files);
cd DataMagnitude\
data{:} = zeros(num_files-2);
for i = 3:num_files
   file = Files(i).name;
   data{i-2} = readmatrix(file);
end
cd ..
FilesU = dir('Bxcomponent');
num_files = length(FilesU);
cd Bxcomponent\
dataU{:} = zeros(num_files-2);
for i = 3:num_files
   fileU = FilesU(i).name;
   dataU{i-2} = readmatrix(fileU);
end
%num_files = length(FilesU);
cd ..
FilesV = dir('Bycomponent');
cd Bycomponent\
dataV{:} = zeros(num_files-2);

for i = 3:num_files
   fileV = FilesV(i).name;
   dataV{i-2} = readmatrix(fileV);
end
cd ..
FilesB = dir('B_mag');
cd B_mag\
dataB{:} = zeros(num_files-2);

for i = 3:num_files
   fileB = FilesB(i).name;
   dataB{i-2} = readmatrix(fileB);
end
cd ..

FilesVort = dir('Vorticity');
cd Vorticity\
dataVort{:} = zeros(num_files-2);

for i = 3:num_files
   fileVort = FilesVort(i).name;
   dataVort{i-2} = readmatrix(fileVort);
end
cd ..

FilesPres = dir('Pressure');
cd Pressure\
dataPres{:} = zeros(num_files-2);

for i = 3:num_files
   filePres = FilesPres(i).name;
   dataPres{i-2} = readmatrix(filePres);
end
cd ..

FilesCurr = dir('CurrentDensity');
cd CurrentDensity\
dataCurr{:} = zeros(num_files-2);

for i = 3:num_files
   fileCurr = FilesCurr(i).name;
   dataCurr{i-2} = readmatrix(fileCurr);
end
cd ..


for i = 1:num_files-2
    AVGPRES(i) = mean(dataPres{i},'all','omitnan')

end
%% Plots

%Plot_Data(num_files,dataB)
Plot_Data(num_files,dataCurr,axes,'Current Density','CurrentDensity_Plots')
%%
Plot_Data(num_files,dataVort,axes,'Vorticity','Vorticity_Plots')
%%
Plot_Data(num_files,dataB,axes,'B-field Intensity','Bfield_Plots')
%%
function [] = Plot_Data(num_files,data,axes,colormap,foldername)

cd(foldername)

for i=1:num_files - 2
    f = figure(i);
    plotname = "Ideal MHD Nu = 0.05 Time = " + 0.1*(i - 1) + " s";
    fname = "Time_0_" + 1*(i - 1) + "_2D_Sine_different_n_0_05";
    s = pcolor(axes(:,1),axes(:,2),data{i});
    s.EdgeColor = 'none';
    c = colorbar;
    c.Label.String = 'Vorticity';
    hold on
    [C,h] = contour(axes(:,1),axes(:,2),data{i},30);
    h.EdgeColor = 'k';
   
    %quiver(axes(1:sc:end,1),axes(1:sc:end,2),dataV{i}(1:sc:end,1:sc:end),dataU{i}(1:sc:end,1:sc:end),'Color','k');
    hold off
    %s.EdgeColor = 'none';
    xlim([0 2*pi])
    xticks([0 pi 2*pi]);
    xticklabels({'0','\pi','2\pi'});
    ylim([0 2*pi])
    yticks([0 pi 2*pi]);
    yticklabels({'0','\pi','2\pi'});
    c = colorbar;
    c.Label.String = colormap;
    title(plotname)
    xlabel("X")
    ylabel("Y")
    saveas(f,fname,'jpg')

end 
cd ..
end