function [ wavfilename,wav_path_name] = find_pcm( FilePath )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
file = dir(FilePath);

DirNULL = 0;

wavfilename = {};
wav_path_name = {};
for i = 1:length(file)
    if(~file(i).isdir)
        if(length(strfind(file(i).name,'.pcm')))
            wavfilename{i} = [FilePath,'\\',file(i).name];
            wav_path_name{i} = file(i).name;
        else
            DirNULL = DirNULL+1;
        end
    else
        DirNULL = DirNULL+1;
    end        
end

wavfilename = wavfilename(DirNULL+1:end);
wav_path_name = wav_path_name(DirNULL+1:end);

end

