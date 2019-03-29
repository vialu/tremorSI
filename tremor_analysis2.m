%% Data extraction tremor

%%%%%% this script is for Extracting tremor data.%%%%%%%

%% Extracting and preprocesing
%
addpath ('/Users/vialundurragaf2/Documents/MATLAB/neurospec20'); % set path to load neurospec toolbox
addpath ('/Users/vialundurragaf2/Documents/MATLAB/fieldtrip-20181119');
addpath ('/Users/vialundurragaf2/Documents/MATLAB/stefslon-exportToPPTX-b2a106d 2');

direc = '//Volumes/shares/DIRFS1/Protocol 93-N-0202/Tremor Studies/';
folder = 'NHP-36/myoclonusStudy_11272018/tremor/';

study = {'rest.eeg','posture.eeg','posture1.eeg','posture15.eeg','posture2.eeg','armsUpStraightOut.eeg','restRotatedAcc.eeg','standing.eeg'};

for sti = 1:length(study);

data_name=strcat(direc,folder,study{sti});


cfg                      = [];
% cfg.channel              = {'FC3'; 'FC1' ;'FCz' ;'FC2' ;'FC4';
% 'C3'; 'C1'; 'Cz'; 'C2'; 'C4';
% 'CP3'; 'CP1'; 'CPz'; 'CP2'; 'CP4'};
cfg.channel              ='all';
cfg.dataset              = data_name;
cfg.bpfilter             ='yes';
cfg.bpfreq               = [2 150];
% cfg.polyremoval          = 'no';
% cfg.dftfilter            = 'yes';
% cfg.dftfreq              = [60 120 180];
data_all                 = ft_preprocessing(cfg);

cfg = [];
cfg.resamplefs = 1000;
data_all = ft_resampledata(cfg,data_all);

% after runing this part, be sure to select the correct channels for EMG
% and ACC in the following section

%%

cfg                      = [];
cfg.channel              = [1:17];
cfg.dataset              = data_name;
cfg.bpfilter             = 'yes';
cfg.bpfreq               = [2 200] ;
% cfg.bsfilter             ='yes';
% cfg.bsfreq               = [58 62];
% cfg.dftfilter            = 'yes';
% cfg.dftfreq              = [60 120 180];
data_eeg                 = ft_preprocessing(cfg);
cfg = [];
cfg.resamplefs = 1000;
data_eeg = ft_resampledata(cfg,data_eeg);






cfg                      = [];
cfg.channel              = [20:23,26:33];
cfg.dataset              = data_name;
cfg.bpfilter             = 'yes';
cfg.bpfreq               = [2 200] ;
% cfg.bsfilter             ='yes';
% cfg.bsfreq               = [58 62];
% cfg.dftfilter            = 'yes';
% cfg.dftfreq              = [60 120 180];
data_emg                 = ft_preprocessing(cfg);
cfg = [];
cfg.resamplefs = 1000;
data_emg = ft_resampledata(cfg,data_emg);


cfg                      = [];
cfg.channel              = [18:19,24:25];
cfg.dataset              = data_name;
cfg.lpfilter             = 'yes';
cfg.lpfreq               = 30 ;
% cfg.bsfilter             ='yes';
% cfg.bsfreq               = [58 62];
% cfg.dftfilter            = 'yes';
% cfg.dftfreq              = [60 120 180];
data_acc                 = ft_preprocessing(cfg);
cfg = [];
cfg.resamplefs = 1000;
data_acc = ft_resampledata(cfg,data_acc);



 cfg                 = [];
 cfg.method          = 'finite';
 cfg.elecfile        = 'electrodes.set';
 cfg.elec            = ft_read_sens(cfg.elecfile);
 
 data_eeg.elec = cfg.elec;
 
 % 
% cfg=[];
% cfg.demean = 'yes';
% cfg.viewmode = 'vertical';
% cfg.layout = 'acticap-64ch-standard2.mat';
% art = ft_databrowser(cfg, data_acc);


%% EMG reconstruction ( +/- channel substraction)

em1 = reshape(data_emg.trial{1,1},2, []);
em2 = em1(1,:)-em1(2,:);
em3 = reshape(em2,6,[]);

data_emg.trial{1,1} = em3;
data_emg.label =   {'Flexor_r';'extensor_r'; 'flexor_l';'extensor_l'; 'bicep'; 'triceps' };

cfg                      = [];
cfg.demean               = ' yes';
cfg.rectify              = 'yes';
data_emg                 = ft_preprocessing(cfg, data_emg);

em1 = reshape(data_acc.trial{1,1},2, []);
em2 = em1(1,:)-em1(2,:);
em3 = reshape(em2,2,[]);

data_acc.trial{1,1} = em3;
data_acc.label =   {'ACC_r','ACC_l'};


em4 = [data_acc.trial{1,1}; data_emg.trial{1,1}];


data_emg.trial{1,1} = em4;
data_emg.label =   {'ACC_r';'ACC_l';'flexor_r';'extensor_r'; 'flexor_l';'extensor_l';'bicep'; 'triceps' };

acc_ch = {'ACC_r';'ACC_l'};

% cfg=[];
% % cfg.mychanscale         = .1;
% % cfg.mychan              = acc_ch;
% cfg.viewmode = 'vertical';
% ft_databrowser(cfg, data_emg);
% set(gcf,'color','w');

data = em4(:,1:60000);
%% Set path to export

cd(strcat(direc,folder));




%% With fft
% 
% Fs = 1000;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L = 60000;
% t=(0:L-1)*T; 
% 
% 
% for ti = 1:6;
%     dataF(ti,:) = fft(data(ti,:));
% end
% 
% frex = Fs*(0:(L/2))/L;
% 
% 
% P2 = abs(dataF/L);
% P1 = P2(:,1:L/2+1);
% P1(:,2:end-1) = 2*P1(:,2:end-1);
% 
% seg =[1 25] ;
% 
% 
% figure;
% subplot(3,2,1)
% plot(frex,P1(1,:))
% title(data_emg.label(1))
% set(gca,'xtick',1:2:25);
% xlim(seg)
% subplot(3,2,2)
% plot(frex,P1(2,:))
% title(data_emg.label(2))
% set(gca,'xtick',1:2:25);
% xlim(seg)
% subplot(3,2,3)
% plot(frex,P1(3,:))
% title(data_emg.label(3))
% set(gca,'xtick',1:2:25);
% xlim(seg)
% subplot(3,2,4)
% plot(frex,P1(5,:))
% title(data_emg.label(5))
% set(gca,'xtick',1:2:25);
% xlim(seg)
% subplot(3,2,5)
% plot(frex,P1(4,:))
% title(data_emg.label(4))
% set(gca,'xtick',1:2:25);
% xlim(seg)
% subplot(3,2,6)
% plot(frex,P1(6,:))
% title(data_emg.label(6))
% xlim(seg)
% set(gca,'xtick',1:2:25);
% set(gcf,'color','w');
% suptitle(strcat(study(sti),'_fft'));
% 
% saveas(gcf,char(strcat(study(sti),'_','fft.png')));

%% Wavelet analysis

npnts = 30000;

% specify frequency range
min_freq = 1;
max_freq = 30;
num_frex = 60;

srate = 1000;
% define frequency and time ranges
frex = linspace(min_freq,max_freq,num_frex);


% parameters for complex Morlet wavelets
wavtime  = -4:1/srate:4-1/srate;
half_wav = (length(wavtime)-1)/2;
cycRange = [15 20];
nCycles  = logspace(log10(cycRange(1)),log10(cycRange(end)),num_frex);


% FFT parameters
nWave = length(wavtime);
nData = npnts * 1;
nConv = nWave + nData -1;

% and create wavelets
cmwX = zeros(length(frex),nConv);
for fi=1:length(frex)
    s          = nCycles(fi)/(2*pi*frex(fi)); % FWHM
    cmw        = exp(2*pi*1i*frex(fi).*wavtime).* exp(-wavtime.^2./(2*s^2));% create wavelet
    tempX      = fft(cmw,nConv);% frequency spectrum
    cmwX(fi,:) = tempX/max(tempX);% amplitude-normalized
end

%

nchans = 6;

fCof = zeros(nchans,length(frex),npnts);



for chani =1:nchans
    
    dataX = fft(data(chani,:),nConv);
    
    for fi=1:length(frex)
        as = ifft(dataX.* cmwX(fi,:));
        as = as (half_wav+1 : end - half_wav);
      
    fCof(chani,fi,:) = as;
    
    end
end



figure;
subplot(3,2,1)
contourf(1:30000,frex,abs(squeeze(fCof(1,:,:))),40,'linecolor','none');
title(data_emg.label(1))
subplot(3,2,2)
contourf(1:30000,frex,abs(squeeze(fCof(2,:,:))),40,'linecolor','none');
title(data_emg.label(2))
subplot(3,2,3)
contourf(1:30000,frex,abs(squeeze(fCof(3,:,:))),40,'linecolor','none');
title(data_emg.label(3))
subplot(3,2,4)
contourf(1:30000,frex,abs(squeeze(fCof(5,:,:))),40,'linecolor','none');
title(data_emg.label(5))
subplot(3,2,5)
contourf(1:30000,frex,abs(squeeze(fCof(4,:,:))),40,'linecolor','none');
title(data_emg.label(4))
subplot(3,2,6)
contourf(1:30000,frex,abs(squeeze(fCof(6,:,:))),40,'linecolor','none');
title(data_emg.label(6))
set(gcf,'color','w');
suptitle(strcat(study(sti),'_','wavelet spectrogarm'));


saveas(gcf,char(strcat(study(sti),'_','wavelet spectrogarm.png')));


spec = squeeze(mean((abs(fCof).^2),3));
seg =[1 25] ;

figure;
subplot(3,2,1)
plot(frex,spec(1,:));
set(gca,'xtick',1:2:25);
xlim(seg)
title(data_emg.label(1))
subplot(3,2,2)
plot(frex,spec(2,:));
set(gca,'xtick',1:2:25);
xlim(seg)
title(data_emg.label(2))
subplot(3,2,3)
plot(frex,spec(3,:));
set(gca,'xtick',1:2:25);
xlim(seg)
title(data_emg.label(3))
subplot(3,2,4)
plot(frex,spec(5,:));
set(gca,'xtick',1:2:25);
xlim(seg)
title(data_emg.label(5))
subplot(3,2,5)
plot(frex,spec(4,:));
set(gca,'xtick',1:2:25);
xlim(seg)
title(data_emg.label(4))
subplot(3,2,6)
plot(frex,spec(6,:));
set(gca,'xtick',1:2:25);
xlim(seg)
title(data_emg.label(6))
set(gcf,'color','w');
suptitle(strcat(study(sti),'_','wavelet frequencies'));

saveas(gcf,char(strcat(study(sti),'_','wavelet frequencies.png')));




% Sampling rate
    rate = 1000;

% Define the segment length for the FFT as power of 2. Ej 2^10 (1024 points)
    seg_pwr =11; 

nchans = 6;

neu_coh= zeros(length(nchans),length(nchans),1024,5);

for chani=1:nchans
    for chanj=1:nchans
     chan1 = em4(chani,:);
     chan2 = em4(chanj,:);  

     [f,t,cl,sc]=sp2a2_m1(0,chan1',chan2',rate,seg_pwr);
     
     neu_coh(chani,chanj,:,:) = f;
        
    end
end

ic = repmat(cl.ch_c95, 1024,1); 


col = 1:6;
squ = reshape((1:36),6,6);
row = squ(1,:);
lb  = {'ACC_r';'ACC_l';'flexor_r';'extensor_r'; 'flexor_l';'extensor_l'};

% all by all

neu_plot = reshape(neu_coh,36,1024,5);

limity = [0 1];


figure;
for ii = 1:36;
    subplot(6,6,ii)
    plot(f(:,1),neu_plot(ii,:,4))
    hold on
    plot(f(:,1),ic)
    ylim(limity)
    xlim ([0 30])
    if ismember(ii,col)
        for ti = 1:6
            title(lb(ii));
        end
    elseif ismember(ii,row)
        
        ylabel(lb(find(row ==ii)));
        
    end
   
end
suptitle(strcat(study(sti),'_','coherence'));

saveas(gcf,char(strcat(study(sti),'_','coherence.png')));

end
%% pictures to powerPoint
% add exportToPPTX toolbox


exportToPPTX('new');
files=dir('*.png');

	
    for zz=1:length(files)
        slideNum = exportToPPTX('addslide');
        fprintf('Added slide %d\n',slideNum);
        io_output=imread(files(zz).name);
        exportToPPTX('addpicture',io_output)
    end

newFile = exportToPPTX('saveandclose','results_output');

