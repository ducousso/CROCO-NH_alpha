close all
clear all

%% directory

ncdir = '/local/tmp/1/ducousso/CROCO-NH/Tank_NHMG_1x1_0.000001_0.5period_track_pb_modified_firststep/';
ncdir = '/local/tmp/1/ducousso/CROCO-NH/Tank_NHMG_1x1_0.000001_0.5period_track_pb_modified_firststep_noUV_ADV_2/';
ncdir = '/local/tmp/1/ducousso/CROCO-NH/Tank_NHMG_1x1_0.000001_0.5period_track_pb_modified_firststep_nofreesurface/';

ncdir = '/local/tmp/1/ducousso/CROCO-NH/Tank_NHMG_1x1_0.001_0.5period_track_pb_modified_firststep/';

step = 1;

%% predictor stage

% rhs and advance

fl_u = ncread([ncdir 'fl_u_000_0000' num2str(step) '.nc'],'u');
fl_w = ncread([ncdir 'fl_w_000_0000' num2str(step) '.nc'],'w');
figure; 
subplot(1,2,1); imagesc(squeeze(fl_u(:,32,1:end))); axis xy; axis equal; axis tight; colorbar;
subplot(1,2,2); imagesc(squeeze(fl_w(:,32,2:end-1))); axis xy; axis equal; axis tight; colorbar;

% nh projection

so_uin = ncread([ncdir 'so_uin_000_0000' num2str(1+(step-1)*2) '.nc'],'uin');
so_win = ncread([ncdir 'so_win_000_0000' num2str(1+(step-1)*2) '.nc'],'win');
figure; 
subplot(1,2,1);imagesc(squeeze(so_uin(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uin');
subplot(1,2,2);imagesc(squeeze(so_win(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so win');
% figure;
% %so_uin_diff = squeeze(so_uin(:,32,1:end)) - flipdim(squeeze(so_uin(:,32,1:end)),2);
% so_uin_diff = squeeze(so_uin(:,32,1:end)) + flipdim(squeeze(so_uin(:,32,1:end)),2);
% subplot(1,2,2); imagesc(so_uin_diff); axis xy; axis equal; axis tight; colorbar;
% so_win_diff = squeeze(so_win(:,32,2:end-1)) + flipdim(squeeze(so_win(:,32,2:end-1)),2);
% subplot(1,2,2); imagesc(so_win_diff); axis xy; axis equal; axis tight; colorbar;
% figure;
% subplot(1,2,1);plot(so_uin_diff(32,:),'o-')
% subplot(1,2,2);plot(so_win_diff(32,:),'o-')

so_cA = ncread([ncdir 'so_cA_000_0000' num2str(1+(step-1)*2) '.nc'],'cA');
so_b = ncread([ncdir 'so_b_000_0000' num2str(1+(step-1)*2) '.nc'],'b');
so_p = ncread([ncdir 'so_p_000_0000' num2str(1+(step-1)*2) '.nc'],'p');
so_r = ncread([ncdir 'so_r_000_0000' num2str(1+(step-1)*2) '.nc'],'r');
figure; 
subplot(3,3,1); imagesc(squeeze(so_cA(6,:,32,2:end  ))); axis xy; axis equal; axis tight; colorbar;
subplot(3,3,4); imagesc(squeeze(so_cA(7,:,32,2:end  ))); axis xy; axis equal; axis tight; colorbar;
subplot(3,3,5); imagesc(squeeze(so_cA(1,:,32,2:end-1))); axis xy; axis equal; axis tight; colorbar;
subplot(3,3,7); imagesc(squeeze(so_cA(8,:,32,2:end  ))); axis xy; axis equal; axis tight; colorbar;
subplot(3,3,8); imagesc(squeeze(so_cA(2,:,32,2:end-1))); axis xy; axis equal; axis tight; colorbar;
figure; 
subplot(1,3,1); imagesc(squeeze(so_b(:,32,2:end-1))); axis xy; axis equal; axis tight; colorbar;
subplot(1,3,2); imagesc(squeeze(so_p(:,32,2:end-1))); axis xy; axis equal; axis tight; colorbar;
subplot(1,3,3); imagesc(squeeze(so_r(:,32,2:end-1))); axis xy; axis equal; axis tight; colorbar;
% figure
% %so_b_diff = squeeze(so_b(:,32,2:end-1)) + flipdim(squeeze(so_b(:,32,2:end-1)),2);
% so_b_diff = squeeze(so_b(:,32,2:end-1)) - flipdim(squeeze(so_b(:,32,2:end-1)),2);
% subplot(1,3,1); imagesc(so_b_diff); axis xy; axis equal; axis tight; colorbar;
% %so_p_diff = squeeze(so_p(:,32,2:end-1)) + flipdim(squeeze(so_p(:,32,2:end-1)),2);
% so_p_diff = squeeze(so_p(:,32,2:end-1)) - flipdim(squeeze(so_p(:,32,2:end-1)),2);
% subplot(1,3,2); imagesc(so_p_diff); axis xy; axis equal; axis tight; colorbar;
% %so_r_diff = squeeze(so_r(:,32,2:end-1)) + flipdim(squeeze(so_r(:,32,2:end-1)),2);
% so_r_diff = squeeze(so_r(:,32,2:end-1)) - flipdim(squeeze(so_r(:,32,2:end-1)),2);
% subplot(1,3,3); imagesc(so_r_diff); axis xy; axis equal; axis tight; colorbar;
% figure;
% subplot(1,3,1);plot(so_b_diff(64,:),'o-')
% subplot(1,3,2);plot(so_p_diff(64,:),'o-')
% subplot(1,3,3);plot(so_r_diff(64,:),'o-')

so_uout = ncread([ncdir 'so_uout_000_0000' num2str(1+(step-1)*2) '.nc'],'uout');
so_wout = ncread([ncdir 'so_wout_000_0000' num2str(1+(step-1)*2) '.nc'],'wout');
so_bout = ncread([ncdir 'so_bout_000_0000' num2str(1+(step-1)*2) '.nc'],'bout');
figure;
subplot(1,3,1);imagesc(squeeze(so_uout(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uout');
subplot(1,3,2);imagesc(squeeze(so_wout(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so wout');
subplot(1,3,3);imagesc(squeeze(so_bout(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so bout');
% figure;
% %so_uout_diff = squeeze(so_uout(:,32,1:end)) - flipdim(squeeze(so_uout(:,32,1:end)),2);
% so_uout_diff = squeeze(so_uout(:,32,1:end)) + flipdim(squeeze(so_uout(:,32,1:end)),2);
% subplot(1,3,1); imagesc(so_uout_diff); axis xy; axis equal; axis tight; colorbar;
% %so_wout_diff = squeeze(so_wout(:,32,2:end-1)) + flipdim(squeeze(so_wout(:,32,2:end-1)),2);
% so_wout_diff = squeeze(so_wout(:,32,2:end-1)) - flipdim(squeeze(so_wout(:,32,2:end-1)),2);
% subplot(1,3,2); imagesc(so_wout_diff); axis xy; axis equal; axis tight; colorbar;
% %so_bout_diff = squeeze(so_bout(:,32,2:end-1)) + flipdim(squeeze(so_bout(:,32,2:end-1)),2);
% so_bout_diff = squeeze(so_bout(:,32,2:end-1)) - flipdim(squeeze(so_bout(:,32,2:end-1)),2);
% subplot(1,3,3); imagesc(so_bout_diff); axis xy; axis equal; axis tight; colorbar;
% figure;
% subplot(1,3,1);plot(so_uout_diff(64,:),'o-')
% subplot(1,3,2);plot(so_wout_diff(64,:),'o-')
% subplot(1,3,3);plot(so_bout_diff(64,:),'o-')

% bt2bc coupling

co_uin = ncread([ncdir 'co_uin_000_0000' num2str(1+(step-1)*2) '.nc'],'uin');
co_win = ncread([ncdir 'co_win_000_0000' num2str(1+(step-1)*2) '.nc'],'win');
figure;
subplot(1,2,1);imagesc(squeeze(co_uin(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred co uin');
subplot(1,2,2);imagesc(squeeze(co_win(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred co win');

co_uf_bar = ncread([ncdir 'co_uf_bar_000_0000' num2str(1+(step-1)*2) '.nc'],'uf_bar');
coin_uf1in = ncread([ncdir 'coin_uf1in_000_0000' num2str(1+(step-1)*2) '.nc'],'uf1in');
coin_uf2in = ncread([ncdir 'coin_uf2in_000_0000' num2str(1+(step-1)*2) '.nc'],'uf2in');
coin_ufin = ncread([ncdir 'coin_ufin_000_0000' num2str(1+(step-1)*2) '.nc'],'ufin');
coin_ufout = ncread([ncdir 'coin_ufout_000_0000' num2str(1+(step-1)*2) '.nc'],'ufout');
%coin_diff_uf = ncread([ncdir 'coin_diff_uf_000_0000' num2str(1+(step-1)*2) '.nc'],'diff_uf');
figure; 
subplot(5,1,1); plot(squeeze(co_uf_bar(32,1:end)),'*');title('pred co uf bar');
subplot(5,1,2); plot(squeeze(coin_uf1in(32,2:end)),'*');title('pred co uf1 in');
subplot(5,1,3); plot(squeeze(coin_uf2in(32,2:end)),'*');title('pred co uf2 in');
subplot(5,1,4); plot(squeeze(coin_ufin(32,2:end)),'*');title('pred co uf in');
subplot(5,1,5); plot(squeeze(coin_ufout(32,1:end)),'*');title('pred co uf out');
%subplot(6,1,4); plot(squeeze(coin_diff_uf(32,1:end)),'*');title('pred co diff uf');

co_uout = ncread([ncdir 'co_uout_000_0000' num2str(1+(step-1)*2) '.nc'],'uout');
co_wout = ncread([ncdir 'co_wout_000_0000' num2str(1+(step-1)*2) '.nc'],'wout');
figure;
subplot(1,2,1);imagesc(squeeze(co_uout(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred co uout');
subplot(1,2,2);imagesc(squeeze(co_wout(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred co wout');

%% corrector stage

% rhs and advance

% nh projection

so_uin = ncread([ncdir 'so_uin_000_0000' num2str(1+2+(step-1)*2) '.nc'],'uin'); %care
so_win = ncread([ncdir 'so_win_000_0000' num2str(1+2+(step-1)*2) '.nc'],'win'); %care
figure; 
subplot(1,2,1);imagesc(squeeze(so_uin(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('corr so uin');
subplot(1,2,2);imagesc(squeeze(so_win(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('corr so win');

so_cA = ncread([ncdir 'so_cA_000_0000' num2str(1+2+(step-1)*2) '.nc'],'cA'); %care
so_b = ncread([ncdir 'so_b_000_0000' num2str(1+2+(step-1)*2) '.nc'],'b'); %care
so_p = ncread([ncdir 'so_p_000_0000' num2str(1+2+(step-1)*2) '.nc'],'p'); %care
so_r = ncread([ncdir 'so_r_000_0000' num2str(1+2+(step-1)*2) '.nc'],'r'); %care
figure; 
subplot(3,3,1); imagesc(squeeze(so_cA(6,:,32,2:end  ))); axis xy; axis equal; axis tight; colorbar;
subplot(3,3,4); imagesc(squeeze(so_cA(7,:,32,2:end  ))); axis xy; axis equal; axis tight; colorbar;
subplot(3,3,5); imagesc(squeeze(so_cA(1,:,32,2:end-1))); axis xy; axis equal; axis tight; colorbar;
subplot(3,3,7); imagesc(squeeze(so_cA(8,:,32,2:end  ))); axis xy; axis equal; axis tight; colorbar;
subplot(3,3,8); imagesc(squeeze(so_cA(2,:,32,2:end-1))); axis xy; axis equal; axis tight; colorbar;
figure; 
subplot(1,3,1); imagesc(squeeze(so_b(:,32,2:end-1))); axis xy; axis equal; axis tight; colorbar;
subplot(1,3,2); imagesc(squeeze(so_p(:,32,2:end-1))); axis xy; axis equal; axis tight; colorbar;
subplot(1,3,3); imagesc(squeeze(so_r(:,32,2:end-1))); axis xy; axis equal; axis tight; colorbar;

so_uout = ncread([ncdir 'so_uout_000_0000' num2str(1+2+(step-1)*2) '.nc'],'uout');%care
so_wout = ncread([ncdir 'so_wout_000_0000' num2str(1+2+(step-1)*2) '.nc'],'wout');%care
%so_bout = ncread([ncdir 'so_bout_000_0000' num2str(2+(step-1)*2) '.nc'],'bout');
figure;
subplot(1,2,1);imagesc(squeeze(so_uout(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('corr so uout');
subplot(1,2,2);imagesc(squeeze(so_wout(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('corr so uout');
%subplot(1,3,3); imagesc(squeeze(so_bout(:,32,2:end-1))); axis xy; axis equal; axis tight; colorbar;

% bt2bc coupling

co_uin = ncread([ncdir 'co_uin_000_0000' num2str(2+(step-1)*2) '.nc'],'uin');
co_win = ncread([ncdir 'co_win_000_0000' num2str(2+(step-1)*2) '.nc'],'win');
figure;
subplot(1,2,1);imagesc(squeeze(co_uin(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('corr co uin');
subplot(1,2,2);imagesc(squeeze(co_win(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('corr co win');

co_uf_bar = ncread([ncdir 'co_uf_bar_000_0000' num2str(2+(step-1)*2) '.nc'],'uf_bar');
coin_uf1in = ncread([ncdir 'coin_uf1in_000_0000' num2str(2+(step-1)*2) '.nc'],'uf1in');
coin_uf2in = ncread([ncdir 'coin_uf2in_000_0000' num2str(2+(step-1)*2) '.nc'],'uf2in');
coin_ufin = ncread([ncdir 'coin_ufin_000_0000' num2str(2+(step-1)*2) '.nc'],'ufin');
figure; 
subplot(4,1,1); plot(squeeze(co_uf_bar(32,1:end)),'*');title('pred co uf bar');
subplot(4,1,2); plot(squeeze(coin_uf1in(32,2:end)),'*');title('pred co uf1 in');
subplot(4,1,3); plot(squeeze(coin_uf2in(32,2:end)),'*');title('pred co uf2 in');
subplot(4,1,4); plot(squeeze(coin_ufin(32,2:end)),'*');title('pred co uf in');

co_uout = ncread([ncdir 'co_uout_000_0000' num2str(2+(step-1)*2) '.nc'],'uout');
co_wout = ncread([ncdir 'co_wout_000_0000' num2str(2+(step-1)*2) '.nc'],'wout');
figure;
subplot(1,2,1);imagesc(squeeze(co_uout(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('corr co uout');
subplot(1,2,2);imagesc(squeeze(co_wout(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('corr co wout');

%% finalized step

figure;
subplot(1,2,1);imagesc(squeeze(u(1:end,32,:,step+1))');axis xy;axis equal;axis tight;colorbar;title('step u');
subplot(1,2,2);imagesc(squeeze(w(2:end-1,32,:,step+1))');axis xy;axis equal;axis tight;colorbar;title('step w');

%%
