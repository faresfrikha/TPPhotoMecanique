clear; close all
% Domain definition
x = 0:500; %limite du domaine suivant x
y = 0:500; %limite du domaine suivant y
% Loading image de ref
Img = double(imread('Img_Speckle_Def_02.tif')); %Importer l'image de référence

Marquage_x = @(x) interp1(-1:201,Img(100,100 + (-1:201)),x,'spline');
Marquage_y = @(y) interp1(-1:201,Img(100 + (-1:201) ,100 ),y,'spline');


% Signal definition
Signal_Ref_x = Marquage_x(x);
Signal_Ref_y = Marquage_y(y);
Signal_Ref_2D = [Signal_Ref_x;Signal_Ref_y];
figure(1); subplot 211; hold on; box on; plot(x,Signal_Ref_x); xlim([0,200]); ylim([0,255]); xlabel('x [px]'); title('Signal x  [GV]'); legend('Reference signal'); subplot 212; imagesc(Signal_Ref_x); colormap(gray); set(gca,'ytick',[]); xlabel('x [px]');
figure(2); subplot 211; hold on; box on; plot(y,Signal_Ref_y); xlim([0,200]); ylim([0,255]); xlabel('y [px]'); title('Signal y [GV]'); legend('Reference signal'); subplot 212; imagesc(Signal_Ref_y); colormap(gray); set(gca,'ytick',[]); xlabel('y [px]');

% Loading image def
Img = double(imread('Img_Speckle_Def_04.tif')); %Importer l'image déformée

Marquage_x = @(x) interp1(-1:201,Img(100,100 + (-1:201)),x,'spline');
Marquage_y = @(y) interp1(-1:201,Img(100 + (-1:201) ,100 ),y,'spline');

% Signal definition
Signal_Def_x = Marquage_x(x);
Signal_Def_y = Marquage_y(y);
Signal_Def_2D = [Signal_Def_x;Signal_Def_y];
figure(1); subplot 211; plot(x,Signal_Def_x);  legend('Reference signal','Deformed signal')
figure(2); subplot 211; plot(y,Signal_Def_y);  legend('Reference signal','Deformed signal')
ImgSize = numel(Signal_Def_x);
x = 1:ImgSize;
y = 1:ImgSize;
% Kinematic defition 
Subset_Size = 17; %On fixe la taille de 'subset' en 17
Subset_xy    = -(Subset_Size-1)/2:(Subset_Size-1)/2; %Choix de fonction de forme 1 pour un mvmt de translation 2 pour une déformation 
[Subset_x,Subset_y] = meshgrid(Subset_xy);
Subset_Phi1x = ones(1,Subset_Size);
Subset_Phi1y = zeros(1,Subset_Size);
Subset_Phi2x = zeros(1,Subset_Size);
Subset_Phi2y = ones(1,Subset_Size);
Subset_Phix = [Subset_Phi1x(:).';Subset_Phi2x(:).'];
Subset_Phiy = [Subset_Phi1y(:).';Subset_Phi2y(:).'];

% Computing the displacement

DICu  = nan(2,ImgSize);
Gradf = gradient(Signal_Ref_2D);
for iPixel = Subset_Size:ImgSize-Subset_Size
    
       Loci = iPixel + Subset_xy;
       
       
        L    = Gradf(1,Loci).*Subset_Phix + Gradf(2,Loci).*Subset_Phiy ;
        M    = L*L.';
        Locu_x = zeros(1,Subset_Size);
        Locu_y = zeros(1,Subset_Size);
        for iLoop = 1:10
        BackDeformed_x = interp1(x,Signal_Def_x,x(Loci)+Locu_x);
        BackDeformed_y = interp1(y,Signal_Def_y,y(Loci)+Locu_y);
        b_x            = L*(Signal_Ref_x(Loci) - BackDeformed_x).';
        b_y            = L*(Signal_Ref_y(Loci) - BackDeformed_y).';
        dDoF_x         = M\b_x;
        dDoF_y         = M\b_y;
        Locu_x         = Locu_x + dDoF_x.'*Subset_Phix;
        Locu_y         = Locu_y + dDoF_y.'*Subset_Phiy;
        if norm(dDoF_x) < 1e-3; break; 
        elseif norm(dDoF_y) < 1e-3; break;
        end
        end
        DICu(1,iPixel) = Locu_x((Subset_Size-1)/2);
        DICu(2,iPixel) = Locu_y((Subset_Size-1)/2);
   
end
figure(4); plot(x,DICu(1,:)); legend('x Retrieved by DIC')
figure(5); plot(y,DICu(2,:)); legend('y Retrieved by DIC')
