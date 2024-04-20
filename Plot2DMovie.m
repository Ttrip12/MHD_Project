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
%FilesU = dir('Bxcomponent');
FilesU = dir('DataU');
num_files = length(FilesU);
%cd Bxcomponent\
cd DataU\
dataU{:} = zeros(num_files-2);
for i = 3:num_files
   fileU = FilesU(i).name;
   dataU{i-2} = readmatrix(fileU);
end
%num_files = length(FilesU);
cd ..
%FilesV = dir('Bycomponent');
FilesV = dir('DataV');
cd DataV\
%cd Bycomponent\
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

FilesPressure = dir('Pressure');
cd Pressure\
dataPressure{:} = zeros(num_files-2);

for i = 3:num_files
   filePressure = FilesPressure(i).name;
   dataPressure{i-2} = readmatrix(filePressure);
end
cd ..

%% Video
v = VideoWriter('Vorticity Resistive with pressure', 'MPEG-4');
v.FrameRate = 30;  % arbitrary
open(v)
sc = 8;
f=figure;
pause(0.2) % let plot wake up
for i=1:num_files - 2
    s = pcolor(axes(:,1),axes(:,2),dataVort{i});
    s.EdgeColor = 'none';
    c = colorbar;
    caxis([0 1.5]);
    c.Label.String = 'B-Field Intensity';
    hold on
    [C,h] = contour(axes(:,1),axes(:,2),dataVort{i},30);
    h.EdgeColor = 'k';
   
    %quiver(axes(1:sc:end,1),axes(1:sc:end,2),dataV{i}(1:sc:end,1:sc:end),dataU{i}(1:sc:end,1:sc:end),'Color','k');
    hold off
   
    xlim([0 2*pi])
    xticks([0 pi 2*pi]);
    xticklabels({'0','\pi','2\pi'});
    ylim([0 2*pi])
    yticks([0 pi 2*pi]);
    yticklabels({'0','\pi','2\pi'});

    im = getframe(f);
    writeVideo(v,im)
end
close(v)
