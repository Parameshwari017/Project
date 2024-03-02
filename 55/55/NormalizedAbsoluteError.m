
function NAE = NormalizedAbsoluteError(origImg, distImg)

origImg = double(origImg);
distImg = double(distImg);

error = origImg - distImg;

NAE = sum(sum(abs(error))) / sum(sum(origImg));