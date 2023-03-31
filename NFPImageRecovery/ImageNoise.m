clear all

counter = 0;
for SNR = 5:5:25
counter = counter + 1;

colorimage = imread('Sonic.png');
[rows, columns, numberOfColorBands] = size(colorimage);

% Display the original image.
subplot(2, 3, 1);
imshow(colorimage, []);
fontSize = 15;
title('Original Image', 'FontSize', fontSize, 'Interpreter', 'None');
noise = randn(rows,columns);
for s=1:3   
grayimage = colorimage(:, :, s); % Take green channel.
grayimage = im2double(grayimage);
noise = (norm(grayimage)/10^(SNR/10))*noise/norm(noise);
for i = 1:rows
    for j = 1:columns
        grayimage(i,j) = grayimage(i,j) + noise(i,j);
    end
end
rgb(:,s) = reshape(im2double(grayimage),[],1);
end

Imageestimate = cat(3, reshape(rgb(:,1),rows,columns), reshape(rgb(:,2),rows,columns), reshape(rgb(:,3),rows,columns));
subplot(2, 3, counter+1);
imshow(Imageestimate, []);
fontSize = 15;
title(sprintf('SNR = %d db',SNR), 'FontSize', fontSize, 'Interpreter', 'None');
end
