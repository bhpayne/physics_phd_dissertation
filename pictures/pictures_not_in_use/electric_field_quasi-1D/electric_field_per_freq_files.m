% 20091207
% Ben Payne, Alexey Yamilov
% input-channel resolved analysis for a single frequency, single gain, single rlz

% NOTE: f90 compiled using
% gfortran -L/usr/lib -llapack quasi1d_vary_active_jwpw_in_body.f90
% 1000 freq took 1.5 hours on /media/primary

clc
clear
format long;
tic;

set(0,'DefaultAxesPosition',[0.05 0.05 0.9 0.95])
set(0,'defaultaxesfontsize',20);
set(0,'defaultaxesfontname','times');
set(0,'defaultaxeslinewidth',1);
set(0,'defaultaxesbox','on');
set(0,'defaulttextfontname','times');
set(0,'defaulttextfontsize',20);
set(0,'defaulttextlinewidth',1);
set(0,'defaultlinelinewidth',2);
set(0,'defaultlinemarkersize',10);

% system parameters
L=200;
W=10.25;
abstep=10;

% variables that are functions of the system parameters
n_open=floor(2*W);
n_closed=0;
n_interior = n_open+n_closed;
nz = L*abstep+1;
y = 0:W/floor(abstep*W):W;

freq = load('out_frequencies.dat');
ge = load('out_uf_gE.dat');

AB_comp = zeros(nz,n_interior,n_open);
tic;
for fli=1:size(freq,1)
    AB_raw = zeros(nz,4*n_interior);
    AB_raw_comp = zeros(nz*n_interior,n_interior);
    if (fli<10)
    AB_raw(1:nz*n_interior,1:(n_interior*4)) = load(['out_AB_freq_0000000',num2str(fli),'.dat']);
    elseif (fli<100)
    AB_raw(1:nz*n_interior,1:(n_interior*4)) = load(['out_AB_freq_000000',num2str(fli),'.dat']);
    elseif (fli<1000)
    AB_raw(1:nz*n_interior,1:(n_interior*4)) = load(['out_AB_freq_00000',num2str(fli),'.dat']);
    elseif (fli<10000)
    AB_raw(1:nz*n_interior,1:(n_interior*4)) = load(['out_AB_freq_0000',num2str(fli),'.dat']);
    else
       disp('loading index problem. pausing')
       pause
    end
    disp(fli)
    toc
    AB_raw_comp(1:nz*n_interior,1:n_interior) = ...
                AB_raw(1:nz*n_interior,             1:n_interior)  +1i*AB_raw(1:nz*n_interior,  n_interior+1:2*n_interior)+ ...
                AB_raw(1:nz*n_interior,2*n_interior+1:3*n_interior)+1i*AB_raw(1:nz*n_interior,3*n_interior+1:4*n_interior);
    for inputchan=1:n_open
        %[1+nz*(inputchan-1) nz+nz*(inputchan-1)]
        AB_comp(1:nz,1:n_interior,inputchan) = AB_raw_comp(1+nz*(inputchan-1):nz+nz*(inputchan-1),1:n_interior);
    end
    AB_raw = zeros(nz,4*n_interior);
    AB_raw_comp = zeros(nz*n_interior,n_interior);
    %clear AB_raw AB_raw_comp;

    k_para(1,1:n_interior) = sqrt( (2*pi*freq(fli))^2 - ((1:n_interior)*pi/W).^2 );
    chi = zeros(n_interior,size(y,2));
    chi(1:n_interior,1:size(y,2)) = sqrt(2/W)*sin([1:n_interior]'*y*pi/W);
    E = zeros((L*10+1),size(y,2),n_open);
    [xm,ym] = meshgrid(0:.1:L,y);

    for inputchan=1:n_open,
      for z=1:(L*10+1)
          E(z,:,inputchan) = squeeze(AB_comp(z,1:n_interior,inputchan))*chi(1:n_interior,:);
      end
      %disp(inputchan)
    end

    E_pw = zeros((L*10+1),size(y,2));
    for inputchan=1:n_open,
      E_pw=E_pw+E(:,:,inputchan)*(mod(inputchan,2)/inputchan);
    end

    figure('Position',[0 -50 1200 400]); hold on;
    subplot(2,1,1);pcolor(xm,ym,abs(E_pw(:,:)').^2); axis equal tight;shading flat;drawnow; 
    title(['field intensity, freq = ',num2str(freq(fli))]); xlabel('L/\lambda'); ylabel('W/\lambda');
    subplot(2,1,2);plot(freq,ge(:,1)); hold on; axis([.99 1.01 1 3]);
    plot(freq(fli),ge(fli,1),'MarkerSize',20,'Marker','.','Color',[1 0 0]);
    xlabel('frequency'); ylabel('g')

    if (fli<10)
      saveas(gcf, ['E_field_000',num2str(fli)], 'png')
    elseif (fli<100)
      saveas(gcf, ['E_field_00',num2str(fli)], 'png')
    elseif (fli<1000)
      saveas(gcf, ['E_field_0',num2str(fli)], 'png')
    elseif (fli<10000)
      saveas(gcf, ['E_field_',num2str(fli)], 'png')
    else
      disp('saving figure problem. pausing')
      pause
    end
    close all 
    
end % fli