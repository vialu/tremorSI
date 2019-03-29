clear all; clc; %close all

if ismac == 1
    data_direc = strcat('/Volumes/shares/DIRFS1/Protocol 17-N-0035/PD_ET_algorithm_data/');
else
    data_direc = strcat('\\nindsdirfs\Shares\HMCS/DIRFS1\Protocol 17-N-0035\PD_ET_algorithm_data\');
end
cd(data_direc);

load('stabilityIndexData_allSubjects');

%% plot ALl subjects
side = 1;
    %1 left
    %2 right
    
condi = 2;    
        %{'rest'     ;... %1
        %'posture'   ;... %2
        %'posture1'   ;... %3
        %'posture15'   ;... %4
        %'posture2'};       %5

for ploti = 1:2
    
    figure();

    for subji = 1:30

        i_freq = i_frex{side,1}{subji,condi};
        deltaF = DeltaF{side,1}{subji,condi};
        
        lo = LO{side,1}{subji,condi};
        up = UP{side,1}{subji,condi};
        baZ = BAz{side,1}{subji,condi};
        y_est = Y_est{side,1}{subji,condi};


        if ploti == 1
            subplot(6,5,subji)
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
                %ylim([-.01 .01]);

        else
            subplot(6,5,subji)
                histfit(deltaF,60)
                view([270 90])
        end



    end; clear subji
end; clearvars -except filtDat fileNames MissingData i_frex DeltaF BAz UP LO Y_est

%%  Extract IQR data


side = 1;
    %1 left
    %2 right
    
condi = 2;    
        %{'rest'     ;... %1
        %'posture'   ;... %2
        %'posture1'   ;... %3
        %'posture15'   ;... %4
        %'posture2'};       %5
             
for li = 1:2        
for subji = 1:30        
        
d1 = DeltaF{li,1}{subji,condi}; 
d2 = iqr(d1);

if li==1
iqr_left(subji) = d2;
else
iqr_right(subji) = d2;
end

end
end


figure
subplot(211)
plot(iqr_right,'.')
subplot(212)
plot(iqr_left,'.')


%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot single subject
subji = 1;

side = 1;
    %1 left
    %2 right
    
condi = 2;    
        %{'rest'     ;... %1
        %'posture'   ;... %2
        %'posture1'   ;... %3
        %'posture15'   ;... %4
        %'posture2'};       %5

i_freq = i_frex{side,1}{subji,condi};
deltaF = DeltaF{side,1}{subji,condi};

lo = LO{side,1}{subji,condi};
up = UP{side,1}{subji,condi};
baZ = BAz{side,1}{subji,condi};
y_est = Y_est{side,1}{subji,condi};


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
        ylim([-.01 .01]);

    subplot(2,3,[1,4])
        histfit(deltaF,60)
        view([270 90])
        
        
clearvars -except filtDat fileNames MissingData i_frex DeltaF BAz UP LO Y_est
