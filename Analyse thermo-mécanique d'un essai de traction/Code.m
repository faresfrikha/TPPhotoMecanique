clear all; close all; clc;

% ---force/allongement-------------
load DONNEES_temps_force_allongement.mat;
figure('name','Force/Allongement')
plot(allongement,force);
xlabel('Allongement en mm');
ylabel('Force en N');
grid;
axis([min(allongement) 1.1*max(allongement) min(force) 1.1*max(force)]) ;
% ------température mesurée---
load DONNEES_temperature.mat;

%---------------------------------------
%------theta----------- 
size=size(temperature);
theta = zeros(size(1),size(2));
for i=1:size(1)
    for j=1:size(2)
        theta(i,j) = (temperature(i,j)) - (temperature(i,1)) ;
    end
end

%-----filtrage de theta------%
taille_filtre_temps=17;
taille_filtre_espace=27;
filtre=ones(taille_filtre_espace,taille_filtre_temps)/(taille_filtre_espace*taille_filtre_temps);
theta_filtre=conv2(theta,filtre,'same');
%---enlever l'effet du bord
theta_filtre(1:round((taille_filtre_espace-1)/2),1:size(2))=NaN*ones(size(2),round((taille_filtre_espace-1)/2))';
theta_filtre(size(1)-round((taille_filtre_espace-1)/2)+1:size(1),1:size(2))=NaN*ones(size(2),round((taille_filtre_espace-1)/2))';
theta_filtre(1:size(1),1:round((taille_filtre_temps-1)/2))=NaN*ones(round((taille_filtre_temps-1)/2),size(1))';
theta_filtre(1:size(1),size(2)-round((taille_filtre_temps-1)/2)+1:size(2))=NaN*ones(round((taille_filtre_temps-1)/2),size(1))';
%---afichage de theta filtrée--------
z_en_mm=linspace(0,40,640);
figure('name','theta filtré')
surface(temps,z_en_mm,theta_filtre(:,1:500),'EdgeColor','none');
axis([min(temps) max(temps) min(z_en_mm) max(z_en_mm) min(min(theta_filtre)) max(max(theta_filtre))]) ;
view(0,90);
xlabel('temps en seconde');
ylabel('z en mm');
colorbar

%---afichage de theta sans filtrage-------
z_en_mm=linspace(0,40,640);
figure('name','theta non filtré')
surface(temps,z_en_mm,theta,'EdgeColor','none');
axis([min(temps) max(temps) min(z_en_mm) max(z_en_mm) min(min(theta)) max(max(theta))]) ;
view(0,90);
xlabel('temps en seconde');
ylabel('z en mm');
colorbar

%theta à prendre en considération
theta_f = theta_filtre;

%dériv à l'ordre 1 en temps de theta
d1theta = zeros(size(1),size(2));
dt = 0.01; %pas temprel en sec
%-----points de bords------
for i=1:size(1)
d1theta(i,1) = (1/(2*dt))*(-theta_f(i,3) + 4*theta_f(i,2) - 3*theta_f(i,1));
end
for i=1:size(1)
d1theta(i,size(2)) = (1/(2*dt))*(3*theta_f(i,size(2)) - 4*theta_f(i,size(2)-1) + theta_f(i,size(2)-2));
end
%-----points intermédiaires------
for i=1:size(1)
    for j=2:(size(2)-1)
d1theta(i,j) = (1/(2*dt))*(theta_f(i,j+1) - theta_f(i,j-1));
    end
end

%----filtrage de d1theta-----%
taille_filtre_t1=41;
taille_filtre_z1=51;
filtre1=ones(taille_filtre_z1,taille_filtre_t1)/(taille_filtre_z1*taille_filtre_t1);
d1theta_filtre=conv2(d1theta,filtre1,'same');
%---enlever l'effet du bord
d1theta_filtre(1:round((taille_filtre_espace1-1)/2),1:size(2))=NaN*ones(size(2),round((taille_filtre_espace1-1)/2))';
d1theta_filtre(size(1)-round((taille_filtre_z1-1)/2)+1:size(1),1:size(2))=NaN*ones(size(2),round((taille_filtre_z1-1)/2))';
d1theta_filtre(1:size(1),1:round((taille_filtre_t1-1)/2))=NaN*ones(round((taille_filtre_t1-1)/2),size(1))';
d1theta_filtre(1:size(1),size(2)-round((taille_filtre_t1-1)/2)+1:size(2))=NaN*ones(round((taille_filtre_t1-1)/2),size(1))';
%d1theta à prendre en considération
d1theta_f = d1theta_filtre;

%---afichage de d1theta sans filtrage
z_en_mm=linspace(0,40,640);
figure('name','d1theta/dt non filtré')
surface(temps,z_en_mm,d1theta,'EdgeColor','none');
axis([min(temps) max(temps) min(z_en_mm) max(z_en_mm) min(min(d1theta)) max(max(d1theta))]) ;
view(0,90);
xlabel('temps en seconde');
ylabel('z en mm');
colorbar

%---afichage de d1theta filtré
z_en_mm=linspace(0,40,640);
figure('name','d1theta/dt filtré')
surface(temps,z_en_mm,d1theta_f,'EdgeColor','none');
axis([min(temps) max(temps) min(z_en_mm) max(z_en_mm) min(min(d1theta_f)) max(max(d1theta_f))]) ;
view(0,90);
xlabel('temps en seconde');
ylabel('z en mm');
colorbar

%---dériv à l'ordre 2 en espace de theta----
L = 40*10^-3; %longueur utile de l'éprouvette en m
d2theta = zeros(size(1),size(2));
dz = L/640; %pas en espace (m)
%---points de bords----------
for j=1:size(2)
   d2theta(1,j) = (theta_f(3,j) - 2*theta_f(2,j) + theta_f(1,j))/(dz^2);
end
for j=1:size(2)
   d2theta(size(1),j) = (theta_f(size(1)-2,j) - 2*theta_f(size(1)-1,j) + theta_f(size(1),j))/(dz^2);
end
%-----points intermédiaires------
    for i=2:(size(1)-1)
        for j=1:size(2)
        d2theta(i,j) = (theta_f(i+1,j) - 2*theta_f(i,j) + theta_f(i-1,j))/(dz^2);
        end
    end
 

%-----filtrage de d2theta------
taille_filtre_t2=51;
taille_filtre_z2=71;
filtre2=ones(taille_filtre_z2,taille_filtre_t2)/(taille_filtre_z2*taille_filtre_t2);
d2theta_filtre=conv2(d2theta,filtre2,'same');
%---enlever l'effet du bord
d2theta_filtre(1:round((taille_filtre_z2-1)/2),1:size(2))=NaN*ones(size(2),round((taille_filtre_z2-1)/2))';
d2theta_filtre(size(1)-round((taille_filtre_z2-1)/2)+1:size(1),1:size(2))=NaN*ones(size(2),round((taille_filtre_z2-1)/2))';
d2theta_filtre(1:size(1),1:round((taille_filtre_t2-1)/2))=NaN*ones(round((taille_filtre_t2-1)/2),size(1))';
d2theta_filtre(1:size(1),size(2)-round((taille_filtre_t2-1)/2)+1:size(2))=NaN*ones(round((taille_filtre_t2-1)/2),size(1))';

%d2theta à prendre en considération
d2theta_f = d2theta_filtre;

%---afichage de d2theta sans filtrage -----
z_en_mm=linspace(0,40,640);
figure('name','d2theta/dz2 non filtré')
surface(temps,z_en_mm,d2theta,'EdgeColor','none');
axis([min(temps) max(temps) min(z_en_mm) max(z_en_mm) min(min(d2theta)) max(max(d2theta))]) ;
view(0,90);
xlabel('temps en seconde');
ylabel('z en mm');
colorbar

%-----afichage de d2theta filtrée------ 
z_en_mm=linspace(0,40,640);
figure('name','d2theta/dz2 filtré')
surface(temps,z_en_mm,d2theta_f,'EdgeColor','none');
axis([min(temps) max(temps) min(z_en_mm) max(z_en_mm) min(min(d2theta_f)) max(max(d2theta_f))]) ;
view(0,90);
xlabel('temps en seconde');
ylabel('z en mm');
colorbar

%------source de chaleur s-----
rho = 7800;
c = 480;
tho = 12.1;
lamda = 40;
s = zeros(size(1),size(2));
for i=1:size(1)
    for j=1:size(2)
s(i,j) = (rho*c)*(d1theta_f(i,j) + (1/tho)*theta_f(i,j) ) - lamda*d2theta_f(i,j) ;
    end
end

%----plot du source de chaleur------
z_en_mm=linspace(0,40,640);
figure('name','source de chaleur avec filtrage')
surface(temps,z_en_mm,s,'EdgeColor','none');
axis([min(temps) max(temps) min(z_en_mm) max(z_en_mm) min(min(s)) max(max(s))]) ;
view(0,90);
xlabel('temps en seconde');
ylabel('z en mm');
colorbar

%---plot de la source de chaleur pr un point p en fct du temps----
p = 400;
figure('name','s en pt')
plot(temps,s(p,:))
xlabel('temps en seconde');
ylabel('s');
grid;
axis([min(temps) max(temps) min(s(p,:)) max(s(p,:))]) ;






