function varargout = object_tracking(varargin)
 
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @object_tracking_OpeningFcn, ...
    'gui_OutputFcn',  @object_tracking_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
  
function object_tracking_OpeningFcn(hObject, eventdata, handles, varargin)
movegui(hObject,'center')
imaqreset
 
handles.fuente=2;
 
set(handles.inicio,'Enable','off');
set(handles.parar,'Enable','off');
set(hObject,'UserData',0)
set(handles.axes1,'XTickLabel',[],'YTickLabel',[])
% Choose default command line output for object_tracking
handles.output = hObject;
 
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = object_tracking_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- FUNCTION TO GET BACKGROUND
function cap_fondo_Callback(hObject, eventdata, handles)
% Reset imaq device
imaqreset
set(hObject,'UserData',0) %User data 0 (1 stop capture)
% Enable "Start" and "Stop" buttons
set(handles.inicio,'Enable','off');
set(handles.parar,'Enable','off');
% Disable current button
set(hObject,'Enable','off');
% Get default source
sel_fuente=handles.fuente;
switch sel_fuente
    % _________________________________________________________________
    case 1        
         
        seudo
        %
        uiwait
         
        global id es_web_ext
        % Determine format depending on the type of camera to use
        if es_web_ext==0
            formato='YUY2_176x144';
        else
            formato='RGB24_320x240';
        end
        try
           
            vid = videoinput('winvideo',id,formato);
           
            guidata(hObject, handles);            
        catch
            % Message on error
            msgbox('Check the connection of the camera','Camera')
            % Remove axis labels
            set(handles.axes1,'XTick',[ ],'YTick',[ ])
        end
        % Specify how often to acquire frame from video stream
        vid.FrameGrabInterval = 1;
        set(vid,'TriggerRepeat',Inf);
        % Start capture
        % _______Get Background_________
        vid.FramesPerTrigger=50;
        start(vid);
        data = getdata(vid,50);
        if es_web_ext==0
            fondo=double(ycbcr2rgb(data(:,:,:,50)));
        else
            fondo=double(data(:,:,:,50));
        end
         
        imshow(uint8(fondo))
        % Reset video object
        stop(vid);
        clear vid
        imaqreset
    case 2 
        [nombre, ruta]=uigetfile('*.avi','SELECT VIDEO AVI');
        if nombre == 0 %If press cancel button, return
            set(hObject,'Enable','on');
            set(handles.inicio,'Enable','on');
            set(handles.parar,'Enable','on');
            return
        end
        %------------------------------------------
        handles.xyloObj = VideoReader(fullfile(ruta,nombre));
        axes(handles.axes1)
        % Show image
        fondo=double(read(handles.xyloObj,10));
        imagesc(uint8(fondo))
        axis image off
        handles.ruta=ruta;
        handles.nombre=nombre;
end
 
handles.backg=fondo;
guidata(hObject,handles)
 
set(hObject,'Enable','on');
set(handles.inicio,'Enable','on');
set(handles.parar,'Enable','on');

 
function inicio_Callback(hObject, eventdata, handles)
% Disable button
set(handles.inicio,'Enable','off');
% Get threshold
umbral=str2double(get(handles.umbral,'String'));
% Verify whether threshold is numerical value
if umbral<0 || isnan(umbral) || isempty(umbral)
    errordlg('Input numerical value','ERROR')
    return
end
% User data 0 (1 stop capture)
set(handles.parar,'UserData',0)
% Get video source
sel_fuente=handles.fuente;
%Get background
fondo=handles.backg;
% _____________________________
switch sel_fuente
    case 1 % WEB CAM        
        global id es_web_ext
        
        if es_web_ext==0
            formato='YUY2_176x144';
        else
            formato='RGB24_320x240';
        end
        try
            vid = videoinput('winvideo',id,formato);
            guidata(hObject, handles);
        catch
            msgbox('Check the connection of the camera','Camera')
            set(handles.axes1,'XTick',[ ],'YTick',[ ])
        end
        % Specify how often to acquire frame from video stream
        vid.FrameGrabInterval = 2;
        set(vid,'TriggerRepeat',Inf);
        %Start capture
        %  try
        start(vid);
        while 1
            if get(handles.parar,'UserData') % Data from "Stop" button
                break
            end
            % Get image
            if es_web_ext==0
                imagen = ycbcr2rgb(getdata(vid,1));
            else
                imagen = getdata(vid,1);
            end
            % Show image
            image(imagen)
            % Convert image to double
            im_ent = double(imagen);
            axis image off
            % Call "compare" function
            compare(im_ent,fondo,umbral);
            drawnow
        end
        stop(vid);
        delete(vid)
        clear vid
        imaqreset
        % _________________________________________________________________
    case 2          
        nFrames = handles.xyloObj.NumberOfFrames;
        for cnt = 1:2:nFrames
            if get(handles.parar,'UserData')% Stop whether "stop" button is pressed
                break
            end
            % Show image
            the_image=read(handles.xyloObj,cnt);
            image(the_image);
            axis image off
            im_ent=double(the_image);
            % Call "compare" function
            compare(im_ent,fondo,umbral);%
            drawnow;
        end        
end
% Enable "Start"  button
set(handles.inicio,'Enable','on');
% User data 0 (1 stop capture)
set(handles.parar,'UserData',0)
guidata(hObject, handles);

% --- STOP BUTTON ---.
function parar_Callback(hObject, eventdata, handles)
set(hObject,'userdata',1)
guidata(hObject, handles)

% --- SELECTION OF VIDEO SOURCE
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
if hObject==handles.video_op %VIDEO AVI
    handles.fuente=2;
else
    handles.fuente=1;    % WEBCAM
end
guidata(hObject,handles)

function umbral_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function umbral_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vid=VideoReader('12.mp4');
 numFrames = vid.NumberOfFrames;
 n=numFrames;
 for i = 1:n
 frames = read(vid,i);
 imwrite(frames,[  int2str(i), '.jpg']);
 im(i)=image(frames);
 end

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
Background=imread('background.jpg');
CurrentFrame=imread('original.jpg');
subplot(2,2,1);imshow(Background);title('BackGround');
subplot(2,2,2);imshow(CurrentFrame);title('Current Frame');
J = std2(Background);
[Background_hsv]=round(rgb2hsv(Background));
[CurrentFrame_hsv]=round(rgb2hsv(CurrentFrame));
Out = bitxor(Background_hsv,CurrentFrame_hsv);
Out=rgb2gray(Out);
[rows, columns]=size(Out);
for i=1:rows
for j=1:columns
if Out(i,j) >0
BinaryImage(i,j)=1;
else
BinaryImage(i,j)=0;
end
end
end
FilteredImage=medfilt2(BinaryImage,[5 5]); 
[L, num]=bwlabel(FilteredImage);
STATS=regionprops(L,'all');
cc=[];
removed=0; 
for i=1:num
dd=STATS(i).Area;
if (dd < 500)
L(L==i)=0;
removed = removed + 1;
num=num-1;
else
end
end
[L2, num2]=bwlabel(L);
[B,~,~,A] = bwboundaries(L2);
subplot(2,2,3),  imshow(L2);title('BackGround Detected');
subplot(2,2,4),  imshow(L2);title('Blob Detected');
 
origImg = imread('11.jpg');
% distImg = imread('1_out.jpg');

noOfDim = ndims(origImg);
if(noOfDim == 3)
    origImg = rgb2gray(origImg);
end

noOfDim = ndims(L2);
if(noOfDim == 3)
    distImg = rgb2gray(L2);
end

origSiz = size(origImg);
distSiz = size(L2);
sizErr = isequal(origSiz, distSiz);
if(sizErr == 0)
    disp('Error: Original Image & Distorted Image should be of same dimensions');
    return;
end

MSE = MeanSquareError(origImg, L2);
disp('Mean Square Error = ');
disp(MSE);

PSNR = PeakSignaltoNoiseRatio(origImg, L2);
disp('Peak Signal to Noise Ratio = ');
disp(PSNR);

NK = NormalizedCrossCorrelation(origImg, L2);
disp('Entropy = ');
disp(NK);

SD = standarddeviation(origImg, L2);
disp('  standard deviation  = ');
disp(SD);

SC = StructuralContent(origImg, L2);
disp('Structural Content  = ');
disp(SC);

MD = MaximumDifference(origImg, L2);
disp('Maximum Difference = ');
disp(MD);

NAE = NormalizedAbsoluteError(origImg, L2);
disp('Overall Cross Entropy = ');
disp(NAE);
hold on;
for k=1:length(B),
if(~sum(A(k,:)))
boundary = B{k};
plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
for l=find(A(:,k));
boundary = B{k};
plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
end
end
end
 


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
clear;
N=32; 
s=rand(N,1);    %first sequence
h=rand(N,1);    %second sequence


for m=0:N-1
    for n=0:N-1
        dC2e(m+1,n+1)=2*cos(pi*m*(2*n+1)/N);
    end
end
dC2e(end+1,1:N)=0; 

for m=1:N
    for n=0:N-1
        dS2e(m+1,n+1)=2*sin(pi*m*(2*n+1)/N); %Since m is starting from 1 and row index of dS2e is starting with m+1, no need to append row of zeros to dS2e matrix because it is already zeros.                                           
    end
end

for m=0:N
    for n=0:N
        
        if n==0 | n==N
            Zeta_n=1/2;
        else
            Zeta_n=1;
        end    
        
        dC1e(m+1,n+1)=2*Zeta_n*cos(2*pi*m*n/N); %This will calculate the downsampled C1e. if we put 1 instead of 2 it will be C1e
    end
end
for m=1:N-1
    for n=1:N-1
        dS1e(m,n)=2*sin(2*pi*m*n/N); %This will calculate the downsampled S1e. if we put 1 instead of 2 it will be S1e
    end
end

dSc2e=dC2e*s;
dHc2e=dC2e*h;

dSs2e=dS2e*s;
dHs2e=dS2e*h;

T1=dSc2e.*dHc2e - dSs2e.*dHs2e;

T2=dSs2e.*dHc2e + dSc2e.*dHs2e;

T2=T2(2:end-1);

T1(1)=T1(1)*2;
T1(N+1)=T1(N+1)*2;

    
dT1c1e=dC1e*T1; 
  
dT2s1e=dS1e*T2; 


y=1/(8*N)*(dT1c1e(2:end)+[dT2s1e;0]);

ydft= real(ifft(fft(s).*fft(h)));

figure(1);
plot(y); hold on;
plot(ydft,'r*');
legend('Proposed method','threshold');
title('drowing video varation');
open gifure.fig;
