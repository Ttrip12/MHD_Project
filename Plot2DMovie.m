clear all
close all

axes = readmatrix('xy.csv');
Files = dir('Pressure');
num_files = length(Files);
cd Pressure\
data{:} = zeros(num_files-2);
for i = 3:num_files
   file = Files(i).name;
   data{i-2} = readmatrix(file);
end
cd ..
%% Video
v = VideoWriter('Pressure.mp4', 'MPEG-4');
v.FrameRate = 60;  % arbitrary
open(v)
f=figure;
pause(0.2) % let plot wake up
for i=1:num_files - 2
    s = surf(axes(:,1),axes(:,2),data{i});
    s.EdgeColor = 'none';
    xlim([0 2*pi])
    xticks([0 pi 2*pi]);
    xticklabels({'0','\pi','2\pi'});
    ylim([0 2*pi])
    yticks([0 pi 2*pi]);
    yticklabels({'0','\pi','2\pi'});
    %zlim([0 40]);
    im = getframe(f);
    writeVideo(v,im)
end
close(v)
