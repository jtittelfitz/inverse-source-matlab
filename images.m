function I = images( imset,res )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    I = cell(2);

    if strcmp(imset,'phantoms')    
        E1 = [01.00, 0.6900, 0.920, 00.00, 00.0000, 000.0; %"skull"
             -0.50, 0.6624, 0.874, 00.00, -0.0184, 000.0; %"head"
             -0.50, 0.1100, 0.310, 00.22, 00.0000, -18.0; %"left eye"
             -0.50, 0.1600, 0.410, -0.22, 00.0000, 018.0; %"right eye"
             00.15, 0.2100, 0.250, 00.00, -0.5050, 000.0; 
             00.50, 0.0460, 0.046, 00.05, 00.4000, 000.0;
             00.60, 0.0460, 0.046, -0.22, -0.1000, 000.0;	 
             00.30, 0.0400, 0.150, -0.00, 00.3000, 000.0]; 

        E2 = [01.00, 0.6900, 0.920, 00.00, 00.0000, 000.0; %"skull"
             -0.50, 0.6000, 0.800, 00.00, -0.0184, 000.0; %"head"
             -0.50, 0.1100, 0.310, 00.22, 00.0000, -18.0; %"left eye"
             -0.50, 0.1600, 0.410, -0.22, 00.0000, 018.0; %"right eye"
             00.15, 0.2100, 0.250, 00.00, 0.3500, 000.0]; 

        E3 = [01.00, 0.6900, 0.920, 00.00, 00.0000, 000.0; %"skull"
             -0.50, 0.6000, 0.800, 00.00, -0.0184, 000.0; %"head"
             0.75, 0.1100, 0.310, 00.22, 00.0000, -18.0; %"left eye"
             -0.25, 0.1600, 0.410, -0.22, 00.0000, 018.0; %"right eye"
             00.15, 0.2100, 0.250, 00.00, 0.3500, 000.0]; 
         
        I{1} = phantom(E2,res);
        I{2} = phantom(E3,res); 
    elseif strcmp(imset,'basic')
        temp = imread('images/circle.png'); s = size(temp,1);
        I{1} = imresize(temp,res/s);
        temp = imread('images/circle.png'); s = size(temp,1);
        I{2} = imresize(temp,res/s);
    elseif strcmp(imset,'comp')
        temp = imread('images/zebras.png'); s = size(temp,1);
        I{1} = imresize(temp,res/s);
        temp = imread('images/cat.png'); s = size(temp,1);
        I{2} = imresize(temp,res/s);
    elseif strcmp(imset,'gauss')
        temp = imread('images/gaussian.png'); s = size(temp,1);
        I{1} = imresize(temp,res/s);
        temp = imread('images/gauss2.png'); s = size(temp,1);
        I{2} = imresize(temp,res/s);
    elseif strcmp(imset,'bonus')
        temp = imread('images/vader.png'); s = size(temp,1);
        I{1} = imresize(temp,res/s);
        temp = imread('images/cat.png'); s = size(temp,1);
        I{2} = imresize(temp,res/s);
    else
        I{1} = phantom('Modified Shepp-Logan',res);
        I{2} = phantom('Modified Shepp-Logan',res);
    end
end

