function s = estimate_noise(I)
    [height,width] = size(I);
    I = double(I);
    M = [1 -2 1; -2 4 -2; 1 -2 1];
    s = sum(abs(conv2(I,M)),'all');
    s = s*sqrt(0.5*pi) ./ (6*(width-2)*(height-2));
end