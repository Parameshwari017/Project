

function SD = standarddeviation(origImg, distImg)

origImg = double(origImg);
distImg = double(distImg);

[M N] = size(origImg);
error = origImg - distImg;

SD = sum(sum(error)) / (M * N);