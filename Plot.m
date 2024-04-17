clear all
close all


Files = dir('Data');
num_files = length(Files);
cd Data
data{:} = zeros(num_files-2);
for i = 3:num_files
   file = Files(i).name;
   data{i-2} = readmatrix(file);
end
cd ..
%% Video
v = VideoWriter('Spikes.mp4', 'MPEG-4');
%v.LosslessCompression = true;
v.FrameRate = 60;  % arbitrary
open(v)
f=figure;

pause(0.2) % let plot wake up
for i=1:num_files - 2
    plot(data{i}(:,1), data{i}(:,2));
    plot_name = "Frame" + i;
    title(plot_name)
    xlim([0 2*pi])
    xticks([0 pi 2*pi]);
    xticklabels({'0','\pi','2\pi'});
    ylim([-1 5])
    xlabel("X")
    ylabel("Velocity")
    data{i}(:,1);
    im = getframe(f);
    writeVideo(v,im)
end
close(v)

