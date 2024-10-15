clear
close all
clc

BaseVideoName = 'video';

for ii=1:4
    VideoName = [BaseVideoName,'_',num2str(ii),'.mp4'];
    Vid{ii} = VideoReader(VideoName);
    FrameRate(ii) = Vid{ii}.FrameRate;
end

v = VideoWriter('MergedVideo','MPEG-4'); 
v.Quality = 100;
v.FrameRate = min(FrameRate);
open(v)

for ii=1:4
    while hasFrame(Vid{ii}) 
        Video = readFrame(Vid{ii}); 
        writeVideo(v,Video) 
    end
end

close(v)