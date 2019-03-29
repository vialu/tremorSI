%mobile version of Viking Data Analysis (version 2)
%made by Patrick McGurrin, HMCS -- 4/17/2018
clear all; clc; %close all

subjNum = 33;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%everything below is automated once the file name is entered

fprintf(1,'\n-----------------------------------------------------------------\n')
fprintf(1,'\n                           Mobile Viking')
% fprintf(1,'\n                     coded by Patrick McGurrin')
% fprintf(1,'\n                           HMCS, NINDS')
fprintf(1,'\n------------------------------------------------------------------\n')
fprintf(1,'\n                                                                  \n')

%load in the txt file into matlab
if subjNum > 10
    subjCount = num2str(subjNum);
else
    subjCount = strcat('0',num2str(subjNum));
end

if ismac == 1
    data_direc = strcat('/Volumes/shares/DIRFS1/Protocol 17-N-0035/00_SCR_PHE data/TS000',subjCount,'/EMG_Acc_data/txtFiles/');
else
    data_direc = strcat('\\nindsdirfs\Shares\DIRFS1\Protocol 17-N-0035\00_SCR_PHE data\TS000',subjCount,'\EMG_Acc_data\txtFiles\');
end
cd(data_direc);

fileNames = {'Posture'   ;...  %1
            'Posture1'   ;... %2
            'Posture15'   ;... %3
            'Posture2'};     %4 
        
for condi = 1:4
    %create text file name
    DatName = strcat('TS000',subjCount,'_',fileNames{condi,1});

    VikingDat = load(strcat(DatName,'.txt')); clear DatName

    %channels 1-3: right accelerometer, right FCR, right ECR
    dataRight{condi,1} = VikingDat(:,2:4);

    %channels 1-3: left accelerometer, left FCR, left ECR
    dataLeft{condi,1} = VikingDat(:,5:7);
end; clear condi

%how many seconds is your trial?
dataTime = 30; 

%channel "zero" - time stamp inserted by Viking
time = VikingDat(:,1);clear VikingDat
sampleRate = length(time)/dataTime;
%%

data = dataLeft{1,1}(:,1);
figure; plot(data)
srate = 1000;
npnts = 30000;
time = 1:npnts;

frex = linspace(0,srate/2, floor(npnts/2) +1);
spec = abs(fft(data));
spec = spec(1:npnts/2+1);

figure; plot(frex,smooth(spec,20));
xlim([1 15])
%%
lpFilt = designfilt('lowpassiir','FilterOrder',2, ...
         'PassbandFrequency',10,'PassbandRipple',0.2, ...
         'SampleRate',1000);


d = designfilt('bandpassiir','FilterOrder',2, ...
    'HalfPowerFrequency1',3,'HalfPowerFrequency2',7, ...
    'SampleRate',1000);

y = filtfilt(lpFilt,data);
y1 = filtfilt(d,y);

figure;
subplot(311);
plot(data);
title('data')
subplot(312);
plot(y);
title('lPass')
subplot(313);
plot(y1);
title('bPass')




zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);  

zCross = zci(y1);
interval = diff(zCross);
i_freq = 1./interval;

deltaF = diff(i_freq);

baZ  = repmat(0,length(deltaF),1);
dfIQ = iqr(deltaF);
df   = dfIQ/2;
up   = repmat(df,length(deltaF),1);
lo   = repmat(df*-1,length(deltaF),1);

c = polyfit(i_freq(1:end-1),deltaF,1);
y_est = polyval(c,i_freq(1:end-1));

figure;
subplot(2,3,[2,3,5,6])
plot(i_freq(1:end-1),deltaF,'.')
xlabel('instantaneus frequency')
ylabel('delta frequency')
hold on
plot(i_freq(1:end-1),baZ, 'b')
hold on
plot(i_freq(1:end-1),up,'g')
hold on
plot(i_freq(1:end-1),lo,'g')
hold on
plot(i_freq(1:end-1),y_est,'r')
subplot(2,3,[1,4])
histfit(deltaF,60)
view([270 90])
