function [ hankelSize ] = chooseHankelSize( hankelSizeOption, T )
%CHOOSEHANKELSIZE Summary of this function goes here
%   Detailed explanation goes here

if hankelSizeOption == 1
    hankelSize = floor(T/2);
end

if hankelSizeOption == 2
    hankelSize = floor(T/3);
end

if hankelSizeOption == 3
    hankelSize = floor(T/4);
end


if hankelSizeOption == 4
    hankelSize = floor(T/100);
end



end

