close all
clear all

%% param

nx = 64;
nz = 64;
L = 10;
H = 10;

%% hor grid

x_r = repmat(((1:1:nx)-0.5) * L/nx,[nz 1])';
x_u = repmat((0:1:nx) * L/nx,[nz 1])';
x_w = repmat(((1:1:nx)-0.5) * L/nx,[nz+1 1])';

%% vert grid

% z_r = repmat(((1:1:nz)-0.5) * H/nz,[nx 1]);
% z_u = repmat(((1:1:nz)-0.5) * H/nz,[nx+1 1]);
% z_w = repmat((0:1:nz) * H/nz,[nx 1]);

ncdir = '/local/tmp/1/ducousso/CROCO-NH/Tank_NHMG_1x1_1em3_2period_track_pb_modified_firststep_noUV_ADV/';
ncdir = '/local/tmp/1/ducousso/CROCO-NH/Tank_NHMG_1x1_1.5em1_2period_track_pb_modified_firststep_noUV_ADV/';
sc_r = ncread([ncdir 'tank_his.nc'],'sc_r');
sc_w = ncread([ncdir 'tank_his.nc'],'sc_w');
Cs_r = ncread([ncdir 'tank_his.nc'],'Cs_r');
Cs_w = ncread([ncdir 'tank_his.nc'],'Cs_w');
hc = ncread([ncdir 'tank_his.nc'],'hc');
h = ncread([ncdir 'tank_his.nc'],'h');
h = squeeze(h(2:end-1,32));
%zeta or not?
zeta = h*0;
%or
zeta = ncread([ncdir 'tank_his.nc'],'zeta'); zeta = squeeze(zeta(2:end-1,32,1));

hinv=1./h;
h2=(h+hc);
h2inv=1./h2;
%z_r
cff=hc*sc_r;
for k=1:nz
    z0=hc*sc_r(k)+Cs_r(k)*h;
    z_r(:,k)=z0.*h./(h2) + zeta.*(1.+z0.*h2inv);
end
z_r = z_r+H;
%z_w
cff=hc*sc_w;
for k=1:nz
    z0=hc*sc_w(k)+Cs_w(k)*h;
    z_w(:,k)=z0.*h./(h2) + zeta.*(1.+z0.*h2inv);
end
z_w(:,nz+1)=zeta;
z_w = z_w+H;

z_u = cat(1,z_r,z_r(end,:));

%% analytic solution

%k = 2*pi/(2*L);   % half wavelength
%k = 2*pi/L;       % one wavelength
k = 2*pi/(0.5*L); % two wavelength
k = 2*pi/(0.25*L); % four wavelength

step = 1:900;
t = step * 0.002;

grav = 9.81;
%amp = 1e-3;
amp = 1.5e-1;

omega_nh = sqrt(grav*k*tanh(k*H));

for kt=1:length(step)
   
   zeta_nh(:,kt) = amp * cos(k*x_r(:,1)) * cos(omega_nh*t(kt));
   zeta_nh(:,kt) = zeta_nh(:,kt) ...
       + amp^2 * k * (2+cosh(k*H))/(4*sinh(k*H)^3) * cos(2*k*x_r(:,1)) * cos(2*omega_nh*t(kt));
             
   u_nh(:,:,kt) = amp * (cosh(k*z_u)/sinh(k*H)) .* sin(k*x_u) * sin(omega_nh*t(kt));
   u_nh(:,:,kt) = u_nh(:,:,kt) ...
       + 3/4 * amp^2 * k * omega_nh * (cosh(2*k*z_u))/(sinh(k*H)^4) .* sin(2*k*x_u) * sin(2*omega_nh*t(kt));
            
   w_nh(:,:,kt) =-amp * omega_nh/sinh(k*H) * cos(k*x_w).*sinh(k*z_w) * sin(omega_nh*t(kt));
   
   %U_nh = amp * omega_nh/sinh(k*H) * sin(k*x_u(:,1))*sinh(k*H)/k * sin(omega_nh*t);
   %ubar_nh = amp * omega_nh/sinh(k*H) * sin(k*x_u(:,1))*sinh(k*H)/k/H * sin(omega_nh*t);
   
   dt_zeta_nh(:,kt) =-amp * omega_nh*cos(k*x_r(:,1)) * sin(omega_nh*t(kt));
   dt_u_nh(:,:,kt) = amp * omega_nh^2/sinh(k*H) * sin(k*x_u).*cosh(k*z_u) * cos(omega_nh*t(kt));
   dt_w_nh(:,:,kt) =-amp * omega_nh^2/sinh(k*H) * cos(k*x_w).*sinh(k*z_w) * cos(omega_nh*t(kt));

end
  
% omega_h = sqrt(grav*H*k^2);
% zeta_h = amp * cos(k*x_r(:,1)) * cos(omega_h*t);
% u_h = amp * grav*k/omega_h * sin(k*x_u).*(z_u./z_u) * sin(omega_h*t);
% w_h =-amp * grav*k^2/omega_h * cos(k*x_w).*z_w * sin(omega_h*t);
% U_h = amp * grav*k/omega_h * sin(k*x_u(:,1))*H * sin(omega_h*t);

%% plot vertical grid

figure;
plot(1:nz,z_r(33,:)-H,'o');
hold on
plot((1:nz+1)-0.5,z_w(33,:)-H,'*');

%% plot time serie

figure;
plot(t,squeeze(u_nh(33,64,:)),'o');title('analytic u');

%% plot snapshot

figure;
plot(x_r(:,1),zeta_nh,'o');title('analytic zeta');

figure;
subplot(1,2,1);
plot(x_u(:,1),squeeze(u_nh(:,end)),'o');title('analytic u');
%plot(squeeze(u_nh(33,:)),squeeze(z_u(33,:))-H,'o');title('analytic u');
subplot(1,2,2);
plot(x_w(:,1),squeeze(w_nh(1,:)),'o');title('analytic w');
%plot(squeeze(w_nh(1,:)),squeeze(z_w(1,:))-H,'o');title('analytic w');

figure;
subplot(1,2,1);
x_r_pc = cat(1,cat(1,zeros(1,nz),x_r),L*ones(1,nz));
x_r_pc = cat(2,x_r_pc,x_r_pc(:,end));
z_w_pc = cat(1,cat(1,z_w(1,:),z_w),z_w(end,:)) ;
u_nh_pc = cat(2,cat(1,u_nh,zeros(1,nz)),zeros(nx+2,1));
pcolor(x_r_pc,z_w_pc-H,double(u_nh_pc));shading flat;
axis xy;axis equal;xlim([0 L]);ylim([-H 0]);colorbar;title('analytic u');
clear x_r_pc z_w_pc u_nh_pc
subplot(1,2,2);
x_u_pc = cat(2,cat(2,x_u,x_u(:,1)),x_u(:,end));
z_r_pc = cat(2,cat(2,zeros(nx,1),z_r),H*ones(nx,1));
z_r_pc = cat(1,z_r_pc,z_r_pc(end,:));
w_nh_pc = cat(2,cat(1,w_nh,zeros(1,nz+1)),zeros(nx+1,1));
pcolor(x_u_pc,z_r_pc-H,double(w_nh_pc));shading flat;
clear x_u_pc z_r_pc w_nh_pc
axis xy;axis equal;xlim([0 L]);ylim([-H 0]);colorbar;title('analytic w');

figure;
subplot(1,2,1);
x_r_pc = cat(1,cat(1,zeros(1,nz),x_r),L*ones(1,nz));
x_r_pc = cat(2,x_r_pc,x_r_pc(:,end));
z_w_pc = cat(1,cat(1,z_w(1,:),z_w),z_w(end,:)) ;
dt_u_nh_pc = cat(2,cat(1,dt_u_nh,zeros(1,nz)),zeros(nx+2,1));
pcolor(x_r_pc,z_w_pc-H,double(dt_u_nh_pc));shading flat;
axis xy;axis equal;xlim([0 L]);ylim([-H 0]);colorbar;title('analytic d_t(u)');
clear x_r_pc z_w_pc dt_u_nh_pc
subplot(1,2,2);
x_u_pc = cat(2,cat(2,x_u,x_u(:,1)),x_u(:,end));
z_r_pc = cat(2,cat(2,zeros(nx,1),z_r),H*ones(nx,1));
z_r_pc = cat(1,z_r_pc,z_r_pc(end,:));
dt_w_nh_pc = cat(2,cat(1,dt_w_nh,zeros(1,nz+1)),zeros(nx+1,1));
pcolor(x_u_pc,z_r_pc-H,double(dt_w_nh_pc));shading flat;
clear x_u_pc z_r_pc dt_w_nh_pc
axis xy;axis equal;xlim([0 L]);ylim([-H 0]);colorbar;title('analytic d_t(w)');

figure;
plot(x_u(:,1),U_nh,'*');
figure;
plot(x_u(:,1),ubar_nh,'*');

% figure;
% plot(zeta_h,'*');title('analytic h zeta');
% figure;
% subplot(1,2,1);imagesc(squeeze(u_h(:,:))');axis xy;axis equal;axis tight;colorbar;title('analytic h u');
% subplot(1,2,2);imagesc(squeeze(w_h(:,:))');axis xy;axis equal;axis tight;colorbar;title('analytic h w');
% figure;
% plot(U_h,'*');
% hold on; plot(sum(u_h(:,:).*repmat((z_w(1,2:end)-z_w(1,1:end-1)),[nx+1 1]),2),'o')
