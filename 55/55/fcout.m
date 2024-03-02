function [ output_args ] = Untitled5( input_args )

%% Initializations

persistent r;
r = zeros(10304);

%error checkif (Rectangle(2,1) > size(I,1))
    Rectangle(2,1) = size(I,1);
end
if (Rectangle(4,2) > size(I,2))
    Rectangle(4,2) = size(I,2);
end

if (Rectangle(1,1)==0)
    Rectangle(1,1) = 1;
end
if (Rectangle(2,1)==0)
    Rectangle(2,1) = 1;
end
if (Rectangle(1,2)==0)
    Rectangle(1,2) = 1;
end
if (Rectangle(4,2)==0)
    Rectangle(4,2) = 1;
end

face = I(Rectangle(1,2):Rectangle(4,2),Rectangle(1,1):Rectangle(2,1));
% SAVES THE DETECTED FACE AS A JPG
imwrite(imresize(face, [200 200]),strcat(num2str(number_face),'.jpg'));

imwrite(imresize(face, [112 92]), 'temp.pgm');
a=imread('temp.pgm');
r=reshape(a,size(a,1)*size(a,2),1);
r = uint8(r); % Convert to unsigned 8 bit numbers to save memory.

v = w;
clear w;

N=20;                                
%% Subtracting the mean from v
O=uint8(ones(1,size(v,2)));
m=uint8(mean(v,2));                 % m is the mean of all images.
vzm=v-uint8(single(m)*single(O));   % vzm is v with the mean removed.

 
L=single(vzm)'*single(vzm);
[V,D]=eig(L);
V=single(vzm)*V;
V=V(:,end:-1:end-(N-1));            

 
cv=zeros(size(v,2),N);
for i=1:size(v,2);
    cv(i,:)=single(vzm(:,i))'*V;    % Each row in cv is the signature for one image.
end


%% Recognition
%  Now, we run the algorithm and see if we can correctly recognize the face.
p=r-m;                              % Subtract the mean
s=single(p)'*V;
z=[];
for i=1:size(v,2)
    z=[z,norm(cv(i,:)-s,2)];
end

no_of_min = 3;
value_min = [];
[a,i]=min(z); %find the most closest one

for k=1:no_of_min
    [a,i]=min(z); %find the most closest one
    z(i) = max(z);
    value_min = [value_min i];
end

mu_min = mean(value_min);

fid2 = fopen(strcat(num2str(number_face),'.txt'), 'w');

if ((abs(value_min(1)-mu_min)<10) & (abs(value_min(2)-mu_min)<10) & (abs(value_min(3)-mu_min)<10))
   
    imwrite(imresize((reshape(v(:,value_min(1)),112,92)), [200 200]),strcat('r',num2str(number_face),'.jpg'));
    cd(strcat('s',num2str(floor(value_min(1)/10)-1)));
    
    fid = fopen('ID.txt', 'r');
    str = fread(fid, '*char')';
    fclose(fid);
    cd ..;
    
    fwrite(fid2, str);
    fid3 = fopen('rec.txt', 'w');
    fwrite(fid3, 'false');
    fclose(fid3);

else
    fwrite(fid2, 'UNKNOWN');
    email(strcat(num2str(number_face),'.jpg'));
    fid3 = fopen('rec.txt', 'w');
    fwrite(fid3, 'true');
    fclose(fid3);
end

fclose(fid2);

if number_face == 10
    number_face = 0;
end
number_face = number_face + 1;
save face_counter.mat number_face
clear

 



end

