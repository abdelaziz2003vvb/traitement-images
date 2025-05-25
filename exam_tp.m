function varargout = exam_tp(varargin)
% EXAM_TP MATLAB code for exam_tp.fig
%      EXAM_TP, by itself, creates a new EXAM_TP or raises the existing
%      singleton*.
%
%      H = EXAM_TP returns the handle to a new EXAM_TP or the handle to
%      the existing singleton*.
%
%      EXAM_TP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXAM_TP.M with the given input arguments.
%
%      EXAM_TP('Property','Value',...) creates a new EXAM_TP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before exam_tp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to exam_tp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help exam_tp

% Last Modified by GUIDE v2.5 17-Feb-2025 21:34:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @exam_tp_OpeningFcn, ...
                   'gui_OutputFcn',  @exam_tp_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before exam_tp is made visible.
function exam_tp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to exam_tp (see VARARGIN)

% Choose default command line output for exam_tp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes exam_tp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = exam_tp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fichier_Callback(hObject, eventdata, handles)
% hObject    handle to fichier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function operateur_ponctulle_Callback(hObject, eventdata, handles)
% hObject    handle to operateur_ponctulle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function bruit_Callback(hObject, eventdata, handles)
% hObject    handle to bruit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function filtre_pass_bas_Callback(hObject, eventdata, handles)
% hObject    handle to filtre_pass_bas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function filtre_pass_haut_Callback(hObject, eventdata, handles)
% hObject    handle to filtre_pass_haut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function filtrage_frequentiel_Callback(hObject, eventdata, handles)
% hObject    handle to filtrage_frequentiel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function morphologie_Callback(hObject, eventdata, handles)
% hObject    handle to morphologie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function point_interet_Callback(hObject, eventdata, handles)
% hObject    handle to point_interet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function harris_Callback(hObject, eventdata, handles)
% hObject    handle to harris (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.courant_data;
%==========================================================================
if(size(img,3)==3)
    display('l''image est en couleur')
    img=rgb2gray(img);
end
%==========================================================================
axes(handles.axes1);
subimage(img);
lambda=0.04;
sigma=1; seuil=200; r=6; w=5*sigma;
[m,n]=size(img)
imd=double(img);
dx=[-1 0 1
    -2 0 2
    -1 0 1]; % deriv?e horizontale : filtre de Sobel
dy=dx'; % deriv?e verticale : filtre de Sobel

g = fspecial('gaussian',max(1,fix(w)), sigma);
Ix=conv2(imd,dx,'same');
Iy=conv2(imd,dy,'same');
Ix2=conv2(Ix.^2, g, 'same');
Iy2=conv2(Iy.^2, g, 'same');
Ixy=conv2(Ix.*Iy, g,'same');

detM=Ix2.*Iy2-Ixy.^2;
trM=Ix2+Iy2;
R=detM-lambda*trM.^2;
%==========================================================================
R1=(1000/(1+max(max(R))))*R;
%==========================================================================          
%****** Seuillage et extraction des points d'intérêt ********
[u,v]=find(R1<=seuil);
nb=length(u);
for k=1:nb
    R1(u(k),v(k))=0;
end
R11=zeros(m+2*r,n+2*r);
R11(r+1:m+r,r+1:n+r)=R1;
[m1,n1]=size(R11);

for i=r+1:m1-r
    for j=r+1:n1-r
        fenetre=R11(i-r:i+r,j-r:j+r);
        ma=max(max(fenetre));
        if fenetre(r+1,r+1)<ma
            R11(i,j)=0;
        end
    end
end

nv=uint8(img); 
axes(handles.axes2);
subimage(nv);

hold on
R11=R11(r+1:m+r,r+1:n+r);
[x,y]=find(R11);
nb=length(x)
plot(y,x,'.r')

% --------------------------------------------------------------------
function susan_Callback(hObject, eventdata, handles)
% hObject    handle to susan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% =====================chargement et conversion de limage==================
im=handles.courant_data;
if(size(im,3)==3)
    display('l''image est en couleur') 
    ig=rgb2gray(im);
end
axes(handles.axes1);
subimage(im);
% =======================conversion de l'image=============================
image=double(im);
[n,m]=size(image);
% =============================données=====================================
rayon=2;
alpha=50;
r=2;
alpha=alpha/100;
mask=zeros(2*rayon+1);
b=ones(rayon+1);
for i=1:rayon+1
    for j=1:rayon+1
        if (rayon==1)
           if(j>i)
            b(i,j)=0;
           end
         else
             if(j>i+1)
            b(i,j)=0;
         end
        end
    end
end
mask(1:rayon+1,rayon+1:2*rayon+1)=b;
mask(1:rayon+1,1:rayon+1)=rot90(b);
mask0=mask;
mask0=flipdim(mask0,1);
mask=mask0+mask;
mask(rayon+1,:)=mask(rayon+1,:)-1;
% ==========================réponse maximale============================
max_reponse=sum(sum(mask));
% =====================balayage de toute l'image===========================
f=zeros(n,m);
for i=(rayon+1):n-rayon
    for j=(rayon+1):m-rayon
  
          image_courant=image(i-rayon:i+rayon,j-rayon:j+rayon);

    image_courant_mask=image_courant.*mask;

         inteniste_cental= image_courant_mask(rayon+1,rayon+1);
         s=exp(-1*(((image_courant_mask-inteniste_cental)/max_reponse).^6));
       somme=sum(sum(s));
%   si le centre du mask est un 0 il faut soustraire les zeros des filtres
                if (inteniste_cental==0)
                    somme=somme-length((find(mask==0)));
                end       
         f(i,j)=somme;           
     end
end
% =============selection et seuillage des points d'interét=================
ff=f(rayon+1:n-(rayon+1),rayon+1:m-(rayon+1));
minf=min(min(ff));
maxf=max(max(f));
fff=f;
d=2*r+1;
temp1=round(n/d);
if (temp1-n/d)<0.5 &(temp1-n/d)>0
temp1=temp1-1;
end
temp2=round(m/d);
if (temp2-m/d)<0.5 &(temp2-m/d)>0
temp2=temp2-1;
end
fff(n:temp1*d+d,m:temp2*d+d)=0;
for i=(r+1):d:temp1*d+d
for j=(r+1):d:temp2*d+d
window=fff(i-r:i+r,j-r:j+r);
window0=window;
[xx,yy]=find(window0==0);
for k=1:length(xx)
window0(xx(k),yy(k))=max(max(window0));
end
minwindow=min(min(window0));
[y,x]=find(minwindow~=window & window<=minf+alpha*(maxf-minf) & window>0);
[u,v]=find(minwindow==window);
if length(u)>1
for l=2:length(u)
fff(i-r-1+u(l),j-r-1+v(l))=0 ;
end
end
if length(x)~=0
for l=1:length(y)
fff(i-r-1+y(l),j-r-1+x(l))=0 ;
end
end
end
end
seuil=minf+alpha*(maxf-minf);
[u,v]=find(minf<=fff & fff<=seuil );
% ==============affichage des resultats====================================
subplot(1,2,2)
imshow(im)
hold on
plot(v,u,'.r','MarkerSize',10)
nombre_de_point_dinteret=length(v)


% --------------------------------------------------------------------
function modele_electrostatique_Callback(hObject, eventdata, handles)
% hObject    handle to modele_electrostatique (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.courant_data;
if(size(img,3)==3)
    display('l''image est en couleur') 
    img=rgb2gray(img);
end
axes(handles.axes1);
subimage(img);
k=0.04; sigma=1; seuil=100; r=6; w=5*sigma;
[m,n]=size(img); imd=double(img);
dxa=[-sqrt(2)/4 0 sqrt(2)/4 ; -1 0 1 ; -sqrt(2)/4 0 sqrt(2)/4];
% dxa=[sqrt(2)/4 0 -sqrt(2)/4; 1 0 -1; sqrt(2)/4 0 -sqrt(2)/4];
dya=dxa'; % derivée verticale
g=fspecial('gaussian',max(1,fix(5*sigma)),sigma); % gaussien
Ixa=conv2(imd,dxa,'same');
Iya=conv2(imd,dya,'same');
Ixa2 = conv2(Ixa.^2, g, 'same'); 
Iya2 = conv2(Iya.^2, g, 'same');
Ixya = conv2(Ixa.*Iya, g,'same');
R=Ixa2.*Iya2-Ixya.^2-k*(Ixa2+Iya2).^2;
R1=(1000/(max(max(R))))*R; %normalisation
[u,v]=find(R1<=seuil);
nb=length(u);
for k=1:nb
R1(u(k),v(k))=0;
end
R11=zeros(m+2*r,n+2*r);
R11(r+1:m+r,r+1:n+r)=R1;
[m1,n1]=size(R11);
for i=r+1:m1-r
for j=r+1:n1-r
fenetre=R11(i-r:i+r,j-r:j+r);
ma=max(max(fenetre));
if fenetre(r+1,r+1)<ma
R11(i,j)=0;
end
end
end
subplot(1,2,2); imshow(img)
hold on
R11=R11(r+1:m+r,r+1:n+r);
[x,y]=find(R11);
nb=length(x)
plot(y,x,'.r')



% --------------------------------------------------------------------
function erosion_Callback(hObject, eventdata, handles)
% hObject    handle to erosion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;

    if size(image,3) == 3
        image = rgb2gray(image);
    end
    image = double(image);
    radius = 4;

   
    [rows, cols] = size(image);
    [x, y] = meshgrid(-radius:radius, -radius:radius);
    se = (x.^2 + y.^2) <= radius^2;  
    erodedI = ones(rows, cols) * 255; 
    for i = 1:rows
        for j = 1:cols
            min_val = 255;  
            for m = -radius:radius
                for n = -radius:radius
                    if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols
                      
                        if se(m+radius+1, n+radius+1) == 1
                            min_val = min(min_val, image(i+m, j+n));
                        end
                    end
                end
            end
            erodedI(i, j) = min_val;
        end
    end
   
    erodedI = uint8(erodedI);
    axes(handles.axes2);
    imshow(erodedI);

% --------------------------------------------------------------------
function dilatation_Callback(hObject, eventdata, handles)
% hObject    handle to dilatation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
if size(image,3) == 3
        image = rgb2gray(image);
end

    image = double(image);
    radius = 4;
    [rows, cols] = size(image);
    [x, y] = meshgrid(-radius:radius, -radius:radius);
    se = (x.^2 + y.^2) <= radius^2; 

    [rows, cols] = size(image);

   
    dilatedI = zeros(rows, cols);
    for i = 1:rows
        for j = 1:cols
            max_val = 0;
            for m = -radius:radius
                for n = -radius:radius
               
                      % Vérifier si (m, n) se trouve à l'intérieur de l'élément structurant (disque)
                    if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols && se(m+radius+1, n+radius+1) == 1
                        max_val = max(max_val, image(i+m, j+n));
                    end
                end
            end
            dilatedI(i, j) = max_val;
        end
    end
    dilatedI = uint8(dilatedI);
    axes(handles.axes2);
    imshow(dilatedI);


% --------------------------------------------------------------------
function ouverture_Callback(hObject, eventdata, handles)
    % hObject    handle to ouverture (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    
    image = handles.courant_data;
    if size(image,3) == 3
        image = rgb2gray(image);
    end
    image = double(image);
    radius = 4;
    [rows, cols] = size(image);
    [x, y] = meshgrid(-radius:radius, -radius:radius);
    se = (x.^2 + y.^2) <= radius^2;  
   
    erodedI = ones(rows, cols) * 255;
    for i = 1:rows
        for j = 1:cols
            min_val = 255;
            for m = -radius:radius
                for n = -radius:radius
                    if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols
                        if se(m+radius+1, n+radius+1) == 1
                            min_val = min(min_val, image(i+m, j+n));
                        end
                    end
                end
            end
            erodedI(i, j) = min_val;
        end
    end

    % Perform Dilation (max operation) on the eroded image
    dilatedI = zeros(rows, cols);
    for i = 1:rows
        for j = 1:cols
            max_val = 0;  
            for m = -radius:radius
                for n = -radius:radius
                    if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols
                        if se(m+radius+1, n+radius+1) == 1
                            max_val = max(max_val, erodedI(i+m, j+n));
                        end
                    end
                end
            end
            dilatedI(i, j) = max_val;
        end
    end
    dilatedI = uint8(dilatedI);
    axes(handles.axes2);
    imshow(dilatedI);



% --------------------------------------------------------------------
function fermetture_Callback(hObject, eventdata, handles)
    % hObject    handle to fermetture (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

   
    image = handles.courant_data;

    if size(image,3) == 3
        image = rgb2gray(image);
    end
    image = double(image);  

   
    radius = 4;
    [rows, cols] = size(image);
    [x, y] = meshgrid(-radius:radius, -radius:radius);
    se = (x.^2 + y.^2) <= radius^2;  
   
    dilatedI = zeros(rows, cols);
    for i = 1:rows
        for j = 1:cols
            max_val = 0;  
            for m = -radius:radius
                for n = -radius:radius
                    if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols
                        if se(m+radius+1, n+radius+1) == 1
                            max_val = max(max_val, image(i+m, j+n));
                        end
                    end
                end
            end
            dilatedI(i, j) = max_val;
        end
    end

    % Perform Erosion (min operation) on the dilated image
    closedI = ones(rows, cols) * 255;  % Initialize to high value (white)
    for i = 1:rows
        for j = 1:cols
            min_val = 255;  % Start with the maximum value
            for m = -radius:radius
                for n = -radius:radius
                    if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols
                        if se(m+radius+1, n+radius+1) == 1
                            min_val = min(min_val, dilatedI(i+m, j+n));
                        end
                    end
                end
            end
            closedI(i, j) = min_val;
        end
    end
    closedI = uint8(closedI);
    axes(handles.axes2);
    imshow(closedI);



% --------------------------------------------------------------------
function gradient_m_Callback(hObject, eventdata, handles)
% hObject    handle to gradient_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function chapeau_hat_Callback(hObject, eventdata, handles)
% hObject    handle to chapeau_hat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function detection_bruit_Callback(hObject, eventdata, handles)
    % hObject    handle to detection_bruit (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
   
    %overture + fermeture 
    % --- overture (  Erosion followed by Dilation) ---
      % erosion
    image = handles.courant_data
    if size(image,3) == 3
        image = rgb2gray(image);
    end
    image = double(image); 
    radius = 4;
    [rows, cols] = size(image);
    [x, y] = meshgrid(-radius:radius, -radius:radius);
    se = (x.^2 + y.^2) <= radius^2; 
    erodedI = ones(rows, cols) * 255;
    
    
    for i = 1:rows
        for j = 1:cols
            min_val = 255; 
            for m = -radius:radius
                for n = -radius:radius
                    if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols
                        if se(m+radius+1, n+radius+1) == 1
                            min_val = min(min_val, image(i+m, j+n));
                        end
                    end
                end
            end
            erodedI(i, j) = min_val;
        end
    end

    % Dilation (max operation)
    dilatedI = zeros(rows, cols);
    for i = 1:rows
        for j = 1:cols
            max_val = 0;  
            for m = -radius:radius
                for n = -radius:radius
                    if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols
                        if se(m+radius+1, n+radius+1) == 1
                            max_val = max(max_val, erodedI(i+m, j+n));
                        end
                    end
                end
            end
            dilatedI(i, j) = max_val;
        end
    end

    % --- fermeture (Dilation followed by Erosion) ---
    % Dilation (max operation) on the dilated image
    closedI = zeros(rows, cols);
    for i = 1:rows
        for j = 1:cols
            max_val = 0;
            for m = -radius:radius
                for n = -radius:radius
                    if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols
                        if se(m+radius+1, n+radius+1) == 1
                            max_val = max(max_val, dilatedI(i+m, j+n));
                        end
                    end
                end
            end
            closedI(i, j) = max_val;
        end
    end

    % Erosion (min operation) on the dilated image
    finalI = ones(rows, cols) * 255;  
    for i = 1:rows
        for j = 1:cols
            min_val = 255; 
            for m = -radius:radius
                for n = -radius:radius
                    if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols
                        if se(m+radius+1, n+radius+1) == 1
                            min_val = min(min_val, closedI(i+m, j+n));
                        end
                    end
                end
            end
            finalI(i, j) = min_val;
        end
    end
    finalI = uint8(finalI);
    axes(handles.axes2);
    imshow(finalI);


% --------------------------------------------------------------------
function FPB_Callback(hObject, eventdata, handles)
% hObject    handle to FPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = handles.courant_data; 
F=fftshift(fft2(I)); 
% %calcul de la taille de l'image; 
M=size(F,1); 
N=size(F,2); 
P=size(F,3); 
H0=zeros(M,N); 
D0=1; 
M2=round(M/2); 
N2=round(N/2); 
H0(M2-D0:M2+D0,N2-D0:N2+D0)=1; 
for i=1:M 
for j=1:N 
G(i,j)=F(i,j)*H0(i,j); 
end 
end 
g=ifft2(G); 
imshow(abs(g),[0,255]); 


% --------------------------------------------------------------------
function FPBB_Callback(hObject, eventdata, handles)
% hObject    handle to FPBB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = handles.courant_data; 
%I = imread('eight.tif'); 
F=fftshift(fft2(I)); 
%calcul de la taille de l'image; 
M=size(F,1); 
N=size(F,2); 
P=size(F,3); 
H0=zeros(M,N); 
D0=1; 
M2=round(M/2); 
N2=round(N/2); 
H0(M2-D0:M2+D0,N2-D0:N2+D0)=1; 
H1(M2-D0:M2+D0,N2-D0:N2+D0)=1; 
n=3; 
for i=1:M 
for j=1:N 
H1(i,j)=1/(1+(H0(i,j)/D0)^(2*n)); 
G(i,j)=F(i,j)*H1(i,j); 
end 
end 
g=ifft2(G); 
imshow(abs(g),[0,255]);


% --------------------------------------------------------------------
function FPH_Callback(hObject, eventdata, handles)
% hObject    handle to FPH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.courant_data; 
%charge; 
F=fftshift(fft2(I)); 
%calcul de la taille de l'image; 
M=size(F,1); 
N=size(F,2); 
P=size(F,3); 
H1=ones(M,N); 
D0=1; 
M2=round(M/2); 
N2=round(N/2); 
H1(M2-D0:M2+D0,N2-D0:N2+D0)=0; 
for i=1:M 
for j=1:N 
G(i,j)=F(i,j)*H1(i,j); 
end 
end 
g=ifft2(G); 
imshow(255-abs(g),[0,255]); 


% --------------------------------------------------------------------
function FPHB_Callback(hObject, eventdata, handles)
% hObject    handle to FPHB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.courant_data; 
F=fftshift(fft2(I)); 
%calcul de la taille de l'image; 
M=size(F,1); 
N=size(F,2); 
P=size(F,3); 
H1=ones(M,N); 
D0=1; 
M2=round(M/2); 
N2=round(N/2); 
H1(M2-D0:M2+D0,N2-D0:N2+D0)=0; 
n=3; 
for i=1:M 
for j=1:N 
H(i,j)=1/(1+(H1(i,j)/D0)^(2*n)); 
G(i,j)=F(i,j)*H(i,j); 
end 
end 
g=ifft2(G); 
imshow(255-abs(g),[0,255]); 


% --------------------------------------------------------------------
function gradient_Callback(hObject, eventdata, handles)
% hObject    handle to gradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
output=image;
%image=rgb2gray(image);

[m,n] = size(image);
output=zeros(size(image)); 
outputhor=zeros(size(image)); 
outputver=zeros(size(image)); 
maskhor = [0,0,0;-1,0,1;0,0,0]; 
maskver = [0,-1,0;0,0,0;0,1,0];
for i=4:(m-3)
   for j=4:(n-3) 
      for k=1:3         
          for l=1:3
            outputhor(i,j) = outputhor(i,j)+image(i-k,j-l)*maskhor(k,l);            
            outputver(i,j) = outputver(i,j)+image(i-k,j-l)*maskver(k,l);          
          end
      end
    end
end
for i=4:(m-3)
for j=4:(n-3)       
    output(i,j)=sqrt(outputhor(i,j)*outputhor(i,j) + outputver(i,j)*outputver(i,j));
end 
end 
output=uint8(output); 

%b=uint8(b);
axes(handles.axes2);
subimage(output);


% --------------------------------------------------------------------
function sobel_Callback(hObject, eventdata, handles)
% hObject    handle to sobel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
output=image;
%image=rgb2gray(image);
[m,n] = size(image);
output=zeros(size(image)); 
outputhor=zeros(size(image)); 
outputver=zeros(size(image)); 

maskhor = [-1,0,1;-2,0,2;-1,0,1]; 
maskver = [-1,-2,-1;0,0,0;1,2,1];

for i=4:m-3
   for j=4:n-3
      for k=1:3          
          for l=1:3
            outputhor(i,j) = outputhor(i,j)+image(i-k,j-l)*maskhor(k,l);             
            outputver(i,j) = outputver(i,j)+image(i-k,j-l)*maskver(k,l);          
          end
      end
    end
end

for i=4:m-3
   for j=4:n-3 
output(i,j)=sqrt(outputhor(i,j)*outputhor(i,j) + outputver(i,j)*outputver(i,j)); 
   end
end
output=uint8(output); 
axes(handles.axes2);
subimage(output);


% --------------------------------------------------------------------
function Prewitt_Callback(hObject, eventdata, handles)
% hObject    handle to Prewitt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
output=image;
%image=rgb2gray(image);
[m,n] = size(image);
output=zeros(size(image)); 
outputhor=zeros(size(image)); 
outputver=zeros(size(image)); 

maskhor = [-1,0,1;-1,0,1;-1,0,1]; 
maskver = [-1,-1,-1;0,0,0;1,1,1];

for i=4:m-3
   for j=4:n-3 
      for k=1:3          
          for l=1:3
            outputhor(i,j) = outputhor(i,j)+image(i-k,j-l)*maskhor(k,l);             
            outputver(i,j) = outputver(i,j)+image(i-k,j-l)*maskver(k,l);          
          end
      end
    end
end

for i=4:m-3
   for j=4:n-3 
output(i,j)=sqrt(outputhor(i,j)*outputhor(i,j) + outputver(i,j)*outputver(i,j)); 
   end
end
 
output=uint8(output); 
axes(handles.axes2);
subimage(output);


% --------------------------------------------------------------------
function prewit_Callback(hObject, eventdata, handles)
% hObject    handle to prewit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
output=image;
%image=rgb2gray(image);

[m,n] = size(image);
output=zeros(size(image)); 
outputhor=zeros(size(image)); 
outputver=zeros(size(image)); 
maskhor = [-1,0,1;-1,0,1;-1,0,1]; 
maskver = [-1,-1,-1;0,0,0;1,1,1];
for i=4:(m-3)
   for j=4:(n-3) 
      for k=1:3         
          for l=1:3
            outputhor(i,j) = outputhor(i,j)+image(i-k,j-l)*maskhor(k,l);            
            outputver(i,j) = outputver(i,j)+image(i-k,j-l)*maskver(k,l);          
          end
      end
    end
end
for i=4:(m-3)
for j=4:(n-3)       
    output(i,j)=sqrt(outputhor(i,j)*outputhor(i,j) + outputver(i,j)*outputver(i,j));
end 
end 
output=uint8(output); 

%b=uint8(b);
axes(handles.axes2);
subimage(output);


% --------------------------------------------------------------------
function roberts_Callback(hObject, eventdata, handles)
% hObject    handle to roberts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
[n,m]=size(image);
image = double(image);

for x=1:n-1
 for y=1:m-1
  b(x,y)= abs(uint8( double(image(x,y))-double(image(x+1,y+1))))+ abs(uint8( double(image(x,y+1)) - double(image(x+1,y))));
 end
end
        %Seuillage
        [n,m]=size(image);
        for i=1:n-1
         for j=1:m-1
          if b(i,j) < 25
            b(i,j)=0;
          end
         end
        end
           %
b=uint8(b); 
axes(handles.axes2);
subimage(b);



% --------------------------------------------------------------------
function laplacien_Callback(hObject, eventdata, handles)
% hObject    handle to laplacien (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
%image=imnoise(imageO,'salt & pepper', 0.05);
[n,m]=size(image);
image = double(image);
%b=image;
[n m]=size(image);
b=zeros(n,m);

M1=[-1 -1 -1;-1 8 -1;-1 -1 -1];
for i=2:n-1
    for j=2:m-1
        V=image((i-1:i+1),(j-1:j+1));
        S=V.*M1;
        b(i,j)=sum(S(:));
    end
end
b=uint8(b);
axes(handles.axes2);
subimage(b);


% --------------------------------------------------------------------
function canny_Callback(hObject, eventdata, handles)
    % hObject    handle to canny (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Récupérer l'image courante
    image = handles.courant_data;

    if size(image, 3) == 3
        image = rgb2gray(image);
    end
    
 
    image = im2double(image);
    
    % 1. Lissage de l'image avec un filtre gaussien
    sigma = 1.4; 
    gauss_filter = fspecial('gaussian', [5 5], sigma);  % Filtre gaussien de taille 5x5
    smoothed_image = imfilter(image, gauss_filter, 'replicate');
    
   
    % (Sobel)
    Gx = [-1 0 1; -2 0 2; -1 0 1];
    Gy = [-1 -2 -1; 0 0 0; 1 2 1];
    
    % Appliquer les convolutions pour obtenir les gradients
    Ix = imfilter(smoothed_image, Gx, 'replicate');
    Iy = imfilter(smoothed_image, Gy, 'replicate');
    

    magnitude = sqrt(Ix.^2 + Iy.^2);
    
   
    orientation = atan2(Iy, Ix);
    
    % 3. Suppression des non-maxima
    [rows, cols] = size(magnitude);
    nms_image = zeros(rows, cols);
    
    % Parcours de l'image pour appliquer la suppression des non-maxima
    for i = 2:rows-1
        for j = 2:cols-1
            % Identifier l'orientation du gradient
            angle = orientation(i,j);
            if ((angle >= 0 && angle < pi/4) || (angle >= 7*pi/4)) % 0° ou 180°
                neighbor1 = magnitude(i, j+1);
                neighbor2 = magnitude(i, j-1);
            elseif (angle >= pi/4 && angle < 3*pi/4)  % 45°
                neighbor1 = magnitude(i+1, j-1);
                neighbor2 = magnitude(i-1, j+1);
            elseif (angle >= 3*pi/4 && angle < 5*pi/4)  % 90°
                neighbor1 = magnitude(i+1, j);
                neighbor2 = magnitude(i-1, j);
            else  % 135°
                neighbor1 = magnitude(i+1, j+1);
                neighbor2 = magnitude(i-1, j-1);
            end
            
            % Comparer la magnitude du pixel central avec ses voisins
            if (magnitude(i,j) >= neighbor1) && (magnitude(i,j) >= neighbor2)
                nms_image(i,j) = magnitude(i,j); % Conserver le pixel
            else
                nms_image(i,j) = 0; % Éliminer le pixel
            end
        end
    end
    

    % Définir les seuils bas et hauts
    T_low = 0.1;  % Seuil bas
    T_high = 0.3;  % Seuil haut
    
    % Appliquer les seuils
    edge_image = nms_image > T_low;
    
    % Relier les pixels de bords faibles aux pixels de bords forts
    for i = 2:rows-1
        for j = 2:cols-1
            if nms_image(i,j) > T_high
                edge_image(i,j) = 1;  % Bords forts
            elseif nms_image(i,j) > T_low && edge_image(i,j) == 1
                edge_image(i,j) = 1;  % Connecter les bords faibles aux forts
            else
                edge_image(i,j) = 0;  % Pas de bord
            end
        end
    end
    
    axes(handles.axes2);
    imshow(edge_image);
  



% --------------------------------------------------------------------
function kirsch_Callback(hObject, eventdata, handles)
% hObject    handle to kirsch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageIn=handles.courant_data;
% [N M L] = size(imageIn);
% g = double( (1/15)*[5 5 5;-3 0 -3; -3 -3 -3] );
% kirschImage = zeros(N,M,8);
% for j = 1:8
%     theta = (j-1)*45;
%     gDirection = imrotate(g,theta,'crop');
%     kirschImage(:,:,j) = conv2(imageIn,gDirection,'same');
% end
% imageOut = zeros(N,M);
% for n = 1:N
%     for m = 1:M
%         imageOut(n,m) = max(kirschImage(n,m,:));
%     end
% end

    x=double(imageIn);


    g1=[5,5,5; -3,0,-3; -3,-3,-3];
    g2=[5,5,-3; 5,0,-3; -3,-3,-3];
    g3=[5,-3,-3; 5,0,-3; 5,-3,-3];
    g4=[-3,-3,-3; 5,0,-3; 5,5,-3];
    g5=[-3,-3,-3; -3,0,-3; 5,5,5];
    g6=[-3,-3,-3; -3,0,5;-3,5,5];
    g7=[-3,-3,5; -3,0,5;-3,-3,5];
    g8=[-3,5,5; -3,0,5;-3,-3,-3];


    x1=imfilter(x,g1,'replicate');
    x2=imfilter(x,g2,'replicate');
    x3=imfilter(x,g3,'replicate');
    x4=imfilter(x,g4,'replicate');
    x5=imfilter(x,g5,'replicate');
    x6=imfilter(x,g6,'replicate');
    x7=imfilter(x,g7,'replicate');
    x8=imfilter(x,g8,'replicate');

    y1=max(x1,x2);
    y2=max(y1,x3);
    y3=max(y2,x4);
    y4=max(y3,x5);
    y5=max(y4,x6);
    y6=max(y5,x7);
    y7=max(y6,x8);
    y=y7;


axes(handles.axes2);
     subimage(uint8(y));


% --------------------------------------------------------------------
function marrhildreth_Callback(hObject, eventdata, handles)
    % hObject    handle to marrhildreth (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    imageIn = handles.courant_data;
    
   
    if size(imageIn, 3) == 3
        im = rgb2gray(imageIn);
    else
        im = imageIn;
    end
    im = im2double(im);

    % Apply the Marr-Hildreth filter (Laplacian of Gaussian)
    gfilter = [0 0 1 0 0;
               0 1 2 1 0;
               1 2 -16 2 1;
               0 1 2 1 0;
               0 0 1 0 0]; 
    smim = conv2(im, gfilter, 'same'); 

   
    [rr, cc] = size(smim);
    zc = zeros([rr, cc]);

   
    for i = 2:rr-1
        for j = 2:cc-1
            if (smim(i,j) > 0)  
              % Vérifier les voisins pour détecter un changement de signe
                if (smim(i,j+1) < 0 && smim(i,j-1) >= 0) || (smim(i,j+1) >= 0 && smim(i,j-1) < 0)
                    zc(i,j) = 1; 
                elseif (smim(i+1,j) < 0 && smim(i-1,j) >= 0) || (smim(i+1,j) >= 0 && smim(i-1,j) < 0)
                    zc(i,j) = 1;  % Mark zero crossing
                elseif (smim(i+1,j+1) < 0 && smim(i-1,j-1) >= 0) || (smim(i+1,j+1) >= 0 && smim(i-1,j-1) < 0)
                    zc(i,j) = 1; 
                elseif (smim(i-1,j+1) < 0 && smim(i+1,j-1) >= 0) || (smim(i-1,j+1) >= 0 && smim(i+1,j-1) < 0)
                    zc(i,j) = 1;  
                end
            end
        end
    end

    
    otpt = im2uint8(zc);

    
    otptth = otpt > 105;
   
    axes(handles.axes2);
    imshow(otptth); % Show the thresholded output







% --------------------------------------------------------------------
function non_lineaire_Callback(hObject, eventdata, handles)
% hObject    handle to non_lineaire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function lineaire_Callback(hObject, eventdata, handles)
% hObject    handle to lineaire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gaussien_Callback(hObject, eventdata, handles)
% hObject    handle to gaussien (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageO=handles.courant_data;
Imgb =imnoise(imageO,'gaussian',0.01);

axes(handles.axes2)
imshow(Imgb);


% --------------------------------------------------------------------
function poivre_et_sel_Callback(hObject, eventdata, handles)
% hObject    handle to poivre_et_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageO=handles.courant_data;
Imgb =imnoise(imageO,'Salt & Pepper', 0.02);

axes(handles.axes2)
imshow(Imgb);



function poisson_Callback(hObject, eventdata, handles)
% hObject    handle to poisson (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageO = handles.courant_data;
Imgb = imnoise(imageO, 'poisson');
axes(handles.axes2)
imshow(Imgb);

% --------------------------------------------------------------------
function localvar_Callback(hObject, eventdata, handles)
% hObject    handle to localvar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the original image
imageO = handles.courant_data;

% Convert image to double precision if needed
if ~isa(imageO, 'double')
    imageO = double(imageO);
end

% Create intensity-dependent noise
% Normalize image to [0,1] range for noise calculation
normalizedImg = imageO/255;
V = 0.01 * normalizedImg; % Variance increases with intensity

% Add local variance noise
Imgb = imnoise(uint8(imageO), 'localvar', V);

% Display result
axes(handles.axes2);
imshow(Imgb);

% Update handles structure
handles.noisy_image = Imgb;
guidata(hObject, handles);

% --------------------------------------------------------------------
function speckle_Callback(hObject, eventdata, handles)
% hObject    handle to speckle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageO = handles.courant_data;
Imgb = imnoise(imageO, 'speckle', 0.04);
axes(handles.axes2)
imshow(Imgb);


% --------------------------------------------------------------------
function contraste_Callback(hObject, eventdata, handles)
% hObject    handle to contraste (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
%[n,m]=size(image);
image = double(image);
%output=image;
 A = min(min(image(:)));
 B = max(max(image(:)));
P=255/(B-A);
L=-P*A;

%ima=imread('cameraman.tif');
[l c]=size(image);
v=image;
for i=1:l
    for j=1:c
     fpixel = image(i,j)*P +L; 
    % on v?rifie que la valeur obtenue est bien dans [0..255]
    if( fpixel>255 )
      fpixel = 255;
    else if( fpixel<0 )
      fpixel = 0;
        end 
    end
    
   v(i,j) = fpixel;
    end
end  
v=uint8(v); 
axes(handles.axes2);
subimage(v);

% --------------------------------------------------------------------
function histogramme_Callback(hObject, eventdata, handles)
% hObject    handle to histogramme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.courant_data;
%I=rgb2gray(img);
        d = length(size(img));
        if d==3
            I = rgb2gray(img);
        elseif d==2
            I = img
        end
axes(handles.axes1);
subimage(I);


[nl nc]=size(I);
v=double(I);
vec=[1:256];
l=0;
for k=0:255 
    for i=1:nl
        for j=1:nc
            if v(i,j)==k 
               l=l+1;
            end
        end
    end
    vec(k+1)=l;
    l=0;
end
axes(handles.axes2);plot(vec);
% --------------------------------------------------------------------
function decalage_additif_Callback(hObject, eventdata, handles)
% hObject    handle to decalage_additif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image=handles.courant_data;
%[n,m]=size(image);
image = double(image);
%output=image;
L=100;

%ima=imread('cameraman.tif');
[l c]=size(image);
image = double(image);
v=image;
for i=1:l
    for j=1:c
     fpixel = image(i,j)+L; 
    % l'encadrage on verifie que la valeur obtenue est bien dans [0..255]
    if( fpixel>255 )
      fpixel = 255;
    else if( fpixel<0 )
      fpixel = 0;
        end 
    end
    
   v(i,j) = fpixel;
    end
end  
v=uint8(v); 
axes(handles.axes2);
subimage(v);


% --------------------------------------------------------------------
function mise_a_echelle_Callback(hObject, eventdata, handles)
% hObject    handle to mise_a_echelle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
%[n,m]=size(image);
image = double(image);
%output=image;

%P=0.5 %si en augmante la valeur de P l'image sera ?clair?
P=1.5;

%ima=imread('cameraman.tif');
[l c]=size(image);
image = double(image);
v=image;
for i=1:l
    for j=1:c
     fpixel = image(i,j)*P; 
    % l'encadrage on v?rifie que la valeur obtenue est bien dans [0..255]
    if( fpixel>255 )
      fpixel = 255;
    else if( fpixel<0 )
      fpixel = 0;
        end 
    end
    
   v(i,j) = fpixel;
    end
end  
v=uint8(v); 
axes(handles.axes2);
subimage(v);


% --------------------------------------------------------------------
function inversion_Callback(hObject, eventdata, handles)
% hObject    handle to inversion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
%[n,m]=size(image);
image = double(image);
[l c]=size(image);
image = double(image);
v=image;
for i=1:l
   for j=1:c
     v(i,j)=-double(image(i,j))+255;
    end
 end 

v=uint8(v); 
axes(handles.axes2);
subimage(v);


% --------------------------------------------------------------------
function sueillage_Callback(hObject, eventdata, handles)
    % hObject    handle to sueillage (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

     % egalisation d'histogramme
     
    image = handles.courant_data;

    % Convertir l'image en niveaux de gris
    if size(image, 3) == 3
        image = rgb2gray(image);
    end
  
    image = double(image);
    
    [rows, cols] = size(image);
    total_pixels = rows * cols;

    % Calculer l'histogramme de l'image
    hist_counts = zeros(1, 256);
    for i = 1:rows
        for j = 1:cols
            intensity = image(i, j) + 1; 
            hist_counts(intensity) = hist_counts(intensity) + 1;
        end
    end

    %  Calculer l'histogramme cumulé
    hist_cdf = cumsum(hist_counts) / total_pixels;

    % Appliquer la transformation de l'égalisation
    equalized_image = zeros(rows, cols);
    for i = 1:rows
        for j = 1:cols
            intensity = image(i, j) + 1;
            equalized_image(i, j) = hist_cdf(intensity) * 255; % Normalisation sur [0,255]
        end
    end

   
    equalized_image = uint8(equalized_image);

    % Afficher l'image égalisée
    axes(handles.axes2);
    imshow(equalized_image);
  
    
    

% --------------------------------------------------------------------
function ouvrir_Callback(hObject, eventdata, handles)
% hObject    handle to ouvrir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile('*.*');
handles.ima = imread(sprintf('%s',path,file));
axes(handles.axes1)
handles.courant_data = handles.ima;
subimage(handles.courant_data);

axes(handles.axes2)
subimage(handles.courant_data);

handles.output = hObject;
guidata(hObject, handles);



% --------------------------------------------------------------------
function enrigester_Callback(hObject, eventdata, handles)
% hObject    handle to enrigester (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
[file,path] = uiputfile('*.png','Enregistrer Votre Image ...');
imwrite(image, sprintf('%s',path,file),'png');


% --------------------------------------------------------------------
function quitter_Callback(hObject, eventdata, handles)
% hObject    handle to quitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)


% --------------------------------------------------------------------
function moyenneur_Callback(hObject, eventdata, handles)
% hObject    handle to moyenneur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function pyramidal_Callback(hObject, eventdata, handles)
% hObject    handle to pyramidal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;

H=(1/81)*[1 2 3 2 1 ; 2 4 6 4 2 ; 3 6 9 6 3 ; 2 4 6 4 2 ; 1 2 3 2 1];

for x = 3 : n-2
    for y = 3 : m-2
          f=image(x-2:x+2,y-2:y+2);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end

b=uint8(b);
axes(handles.axes2);
subimage(b);


% --------------------------------------------------------------------
function conique_Callback(hObject, eventdata, handles)
% hObject    handle to conique (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[ n,m]=size(image);
image = double(image);
b=image;

H=(1/25)*[0 0 1 0 0 ; 0 2 2 2 0 ; 1 2 5 2 1 ; 0 2 2 2 0 ; 0 0 1 0 0];

for x = 3 : n-2
    for y = 3 : m-2
       f=image(x-2:x+2,y-2:y+2);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end

b=uint8(b);
axes(handles.axes2);
subimage(b);

handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function median_Callback(hObject, eventdata, handles)
% hObject    handle to median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
image=double(image);
[n,m]=size(image);
img=image;

for i=2:n-1
    for j=2:m-1
       fenetre=image(i-1:i+1,j-1:j+1);
       v=[fenetre(1,:) fenetre(2,:) fenetre(3,:)];
       sort(v);
       a=median(v);
       img(i,j)=a;
    end
end

b=uint8(img);
handles.ima_traite = b;
axes(handles.axes2);
subimage(b);

handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function forme_blanc_Callback(hObject, eventdata, handles)
% hObject    handle to forme_blanc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image = handles.courant_data;
%pour les formes blanc en utilise la forme suivant
%image=double(image)-double( imopen(image,se));


if size(image, 3) == 3
    image = rgb2gray(image);
end
image = double(image);
radius = 4;
[rows, cols] = size(image);
[x, y] = meshgrid(-radius:radius, -radius:radius);
se = (x.^2 + y.^2) <= radius^2;
erodedI = ones(rows, cols) * 255;  

% **ÉROSION**
for i = 1:rows
    for j = 1:cols
        min_val = 255;  
        for m = -radius:radius
            for n = -radius:radius
                if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols && se(m+radius+1, n+radius+1)
                    min_val = min(min_val, image(i+m, j+n));
                end
            end
        end
        erodedI(i, j) = min_val;
    end
end

dilatedI = zeros(rows, cols);

% **DILATATION**
for i = 1:rows
    for j = 1:cols
        max_val = 0;
        for m = -radius:radius
            for n = -radius:radius
                if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols && se(m+radius+1, n+radius+1)
                    max_val = max(max_val, erodedI(i+m, j+n));
                end
            end
        end
        dilatedI(i, j) = max_val;
    end
end

finalI = image - dilatedI;
nv = uint8(finalI);
axes(handles.axes2);
imshow(nv);


% --------------------------------------------------------------------
function forme_noir_Callback(hObject, eventdata, handles)
% hObject    handle to forme_noir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image = handles.courant_data;

%pour les formes noir  en utilise la forme suivant
%image=double(imclose(image,se))-double(image);

if size(image, 3) == 3
    image = rgb2gray(image);
end
image = double(image);

radius = 4;


[rows, cols] = size(image);


[x, y] = meshgrid(-radius:radius, -radius:radius);
se = (x.^2 + y.^2) <= radius^2; 


dilatedI = zeros(rows, cols);

% **DILATATION**
for i = 1:rows
    for j = 1:cols
        max_val = 0; 
        for m = -radius:radius
            for n = -radius:radius
                if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols && se(m+radius+1, n+radius+1)
                    max_val = max(max_val, image(i+m, j+n));
                end
            end
        end
        dilatedI(i, j) = max_val;
    end
end

% Initialisation de l'image érodée
erodedI = ones(rows, cols) * 255;

% **ÉROSION**
for i = 1:rows
    for j = 1:cols
        min_val = 255; 
        for m = -radius:radius
            for n = -radius:radius
                if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols && se(m+radius+1, n+radius+1)
                    min_val = min(min_val, dilatedI(i+m, j+n));
                end
            end
        end
        erodedI(i, j) = min_val;
    end
end

finalI = erodedI - image;
nv = uint8(finalI);
axes(handles.axes2);
imshow(nv);



% --------------------------------------------------------------------
function gradient_interne_Callback(hObject, eventdata, handles)
% hObject    handle to gradient_interne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%gradient interne =double(image)-double( imerode(image,se));
image = handles.courant_data;
if size(image, 3) == 3
    image = rgb2gray(image);
end
image = double(image);


radius = 4;
[rows, cols] = size(image);
[x, y] = meshgrid(-radius:radius, -radius:radius);
se = (x.^2 + y.^2) <= radius^2;  

erodedI = ones(rows, cols) * 255; 

% **ÉROSION**
for i = 1:rows
    for j = 1:cols
        min_val = 255; 
        for m = -radius:radius
            for n = -radius:radius
                if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols && se(m+radius+1, n+radius+1)
                    min_val = min(min_val, image(i+m, j+n));
                end
            end
        end
        erodedI(i, j) = min_val;
    end
end

% **gradient interne** : Différence entre l’image originale et son érosion
gradient_interne = image - erodedI;


nv = uint8(gradient_interne);
axes(handles.axes2);
imshow(nv);



% --------------------------------------------------------------------
function gradient_externe_Callback(hObject, eventdata, handles)
% hObject    handle to gradient_externe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image = handles.courant_data;
%gradient externe =double(imdilate(image,se))-double(image);
if size(image, 3) == 3
    image = rgb2gray(image);
end
image = double(image);
radius = 4;
[rows, cols] = size(image);
[x, y] = meshgrid(-radius:radius, -radius:radius);
se = (x.^2 + y.^2) <= radius^2;  


dilatedI = zeros(rows, cols); 
% **DILATATION **
for i = 1:rows
    for j = 1:cols
        max_val = 0; 
        for m = -radius:radius
            for n = -radius:radius
                if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols && se(m+radius+1, n+radius+1)
                    max_val = max(max_val, image(i+m, j+n));
                end
            end
        end
        dilatedI(i, j) = max_val;
    end
end

% **gradient externe** : Différence entre la dilatation et l’image originale
gradient_externe = dilatedI - image;


nv = uint8(gradient_externe);
axes(handles.axes2);
imshow(nv);


% --------------------------------------------------------------------
function gradient_morphologique_Callback(hObject, eventdata, handles)
% hObject    handle to gradient_morphologique (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image = handles.courant_data;
%gradient morphologique =double(imdilate(image,se))-double(imerode(image,se));
if size(image, 3) == 3
    image = rgb2gray(image);
end
image = double(image);
radius = 4;
[rows, cols] = size(image);

% élément structurant (disque)
[x, y] = meshgrid(-radius:radius, -radius:radius);
se = (x.^2 + y.^2) <= radius^2; 
dilatedI = zeros(rows, cols); % Pour la dilatation (prendre le max)
erodedI = ones(rows, cols) * 255; % Pour l'érosion (prendre le min)

% **DILATATION**
for i = 1:rows
    for j = 1:cols
        max_val = 0;
        for m = -radius:radius
            for n = -radius:radius
                if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols && se(m+radius+1, n+radius+1)
                    max_val = max(max_val, image(i+m, j+n));
                end
            end
        end
        dilatedI(i, j) = max_val;
    end
end

% **ÉROSION **
for i = 1:rows
    for j = 1:cols
        min_val = 255;
        for m = -radius:radius
            for n = -radius:radius
                if i+m > 0 && i+m <= rows && j+n > 0 && j+n <= cols && se(m+radius+1, n+radius+1)
                    min_val = min(min_val, image(i+m, j+n));
                end
            end
        end
        erodedI(i, j) = min_val;
    end
end

% ** gradient morphologique== difference entre dilation and erosion**
gradient_morphologique = dilatedI - erodedI;
nv = uint8(gradient_morphologique);
axes(handles.axes2);
imshow(nv);




% --------------------------------------------------------------------
function gausien3X3_Callback(hObject, eventdata, handles)
% hObject    handle to gausien3X3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/16)*[1 2 1 ;2 4 2 ; 1 2 1];
for x = 2 : n-1
    for y = 2 : m-1
    f=image(x-1:x+1,y-1:y+1);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.axes2);
    subimage(b);

% --------------------------------------------------------------------
function gausien5X5_Callback(hObject, eventdata, handles)
% hObject    handle to gausien5X5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/256)*[1 4 6 4 1 ; 4 16 24 16 4 ; 6 24 36 24 6 ; 4 16 24 16 4 ; 1 4 6 4 1];
for x = 3 : n-2
    for y = 3 : m-2
  f=image(x-2:x+2,y-2:y+2);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.axes2);
subimage(b);

% --------------------------------------------------------------------
function moyenneur3X3_Callback(hObject, eventdata, handles)
% hObject    handle to moyenneur3X3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/9)*[1 1 1 ; 1 1 1 ; 1 1 1 ];
for x = 2 : n-1
    for y = 2 : m-1
     f=image(x-1:x+1,y-1:y+1);
      v=f.*H;
      b(x,y)=sum(v(:));
      
    end 
end
b=uint8(b);
axes(handles.axes2);
subimage(b);


% --------------------------------------------------------------------
function moyenneur5x5_Callback(hObject, eventdata, handles)
% hObject    handle to moyenneur5x5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/25)*[1 1 1 1 1 ; 1 1 1 1 1 ; 1 1 1 1 1 ; 1 1 1 1 1 ; 1 1 1 1 1];
for x = 3 : n-2
    for y = 3 : m-2
     f=image(x-2:x+2,y-2:y+2);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.axes2);
     subimage(b);
