%%
% <latex> 
% \title{BE 521: Homework 6 \\{\normalsize Spike sorting}\\{\normalsize Spring 2019}} 
% \author{59 points} 
% \date{Due: Thursday, 03/5/2020 11:59pm} 
% \maketitle \textbf{Objective:} Detect and cluster spikes
% </latex>

%% 
% <latex> 
% \begin{center} \author{Andrew Clark \\
%   \normalsize Collaborators: Alex Silva, Vishal Tien \\}
% \end{center} 
% </latex>

%% 
% <latex>
% \subsection*{Overview}
% In this homework, you will do some basic spike sorting using two different datasets. The first (\verb|I521_A0006_D001|) is from a crayfish neuromuscular junction, a good model for human central nervous system synapses\footnote{The sampling rate of this data is 2000 Hz, which is adequate for this homework's instructional purposes but usually inadequate for real spike sorting, which often uses sampling frequencies on the order of 20 kHz.}. Specifically, the data contains two simultaneous recordings: an extracellular recording from the third nerve (channel \verb|nerve|) of a crayfish abdominal ganglion, which contains six spontaneously active motor neurons, and an intracellular recording from the superficial flexor muscle (channel \verb|muscle|) innervated by this nerve. You will attempt to discern relationships between the classes of spike waveforms you extract from the motor nerve trace and elicited potentials seen in the muscle fiber recording. 
% Then, you will revisit a human intracranial EEG recording (\verb|I521_A0006_D002|) and use some of the techniques you've learned in class to build a more automated spike sorter. 
% Note: While spikes may have positive and negative deflections, we will only focus on positive spikes on this homework for simplicity. 
% \section{Spike Detection and Clustering (38 pts)}
% In this section, you will explore some basic filtering and spike thresholding to ultimately compare spike clusters you pick out by eye to those selected by an automated algorithm.
% \begin{enumerate}
%     \item You can assume that the nerve samples have already been low-pass filtered. Here you will high-pass filter in order to remove signals like slow local field potentials and 60 Hz power line noise. Create a 4th order \textit{elliptic filter} with 0.1 dB of ripple in the passband, a stopband 40 dB lower than the peak value in the passband, and a passband edge frequency of 300 Hz (see Matlab's \verb|ellip| function and make sure your give the edge frequency in the correct normalized form). The statement to create this filter (defined by the filter coefficients \verb|b| and \verb|a|) should look something like
%   \begin{lstlisting}
% 	[b,a]=ellip(n,Rp,Rs,Wp,'high')
%   \end{lstlisting}
%   Clearly specify the denominator and numerator coefficients obtained for your filter function. (2pts)
% </latex>

%%


session = IEEGSession('I521_A0006_D001', 'andrewc', '/Users/andrewclark/Downloads/ieeg_password.bin' );
%%
session.data


nr = ceil((session.data.rawChannels(1).get_tsdetails.getEndTime)/1e6*session.data.sampleRate);
allData_muscle = session.data.getvalues(1:nr,1);
allData_nerve=session.data.getvalues(1:nr,2);
%%
durationInUSec = session.data(1).rawChannels(1).get_tsdetails.getDuration;
durationInSec = durationInUSec / 1e6;

%%

%set sample rate
sample_rate=2000;

%create filter
[b, a] = ellip(4,0.1,40, 300/(sample_rate/2),'high');
%%

%As computed by the code above, the denominator coefficients for my filter
%function are 1.000, -1.7432, 1.6167, -0.6559, 0.1430 the numerator coefficients for my filter function are
% 0.3420, -1.2740, 1.8676, -1.2740, 0.3420. 

%% 
% <latex>
%   \item Using the \verb|filter| function and \verb|filtfilt| function, obtain two different filtered outputs of the nerve signal.
%       \begin{enumerate}
%         \item In a 2x1 subplot, plot the first 50 ms of the unfiltered nerve signal in the top subplot; in the bottom subplot, plot the \verb|filter| output in blue and the \verb|filtfilt| output in red. Use a potential range (y-axis) of -20 to 50 millivolts. (4 pts)
% </latex>

%%

%create 50 ms time vector 
time_ms=1/sample_rate:1/sample_rate:0.05;

time_ms=time_ms*1000;

%select 50 ms of nerve data
fifty_ms_nerve=(allData_nerve(1:100))/1000;
%%

%do the filtering using the ellip function output

%filter function
filter_out=filter(b,a, fifty_ms_nerve);

%filtfilt function
filt_filt_out=filtfilt(b,a,fifty_ms_nerve);

%%

%create subplot

figure
subplot(2,1,1)
plot(time_ms, fifty_ms_nerve)
ylim([-20, 50])
ylabel('Voltage (uV)')
xlabel('Time (milliseconds)')
title('Unfiltered Nerve Signal')
subplot(2,1,2)
plot(time_ms, filter_out,'b')
hold on
plot(time_ms, filt_filt_out,'r')
ylabel('Voltage (uV)')
xlabel('Time (milliseconds)')
title('Filtered Nerve Signal')
legend('Filter Output','FiltFilt Output')
ylim([-20, 50])

%% 
% <latex>
%         \item How is the unfiltered signal different from the filtered signal? What is different about the two filtered (red and blue) signals? (2 pts)
% </latex>

%%

% The unfiltered signal is different from the filtered signal in that the points in the unfiltered signal's waveform have   
% a higher amplitude (including much higher spikes) compared to the filtered signal. Additionally, the
% variation in the amplitudes of the points in the unfiltered wavefrom is
% higher than in the filtered waveform (in the filtered waveform the
% points are closer together to each other on the y-axis compared to the unfiltered waveform). The filtered
% (red and blue) signals are different in that the blue signal, generated
% using the MATLAB filter function, is phase shifted compared to both the
% red signal and the unfiltered data. However the red signal, generated
% with the MATLAB filtfilt function, is not phase shifted and the features
% of the signal in the red signal and the unfiltered data are at the same
% x-locations. Additionally, the red signal preserves the large positive voltage spike seen in the unfiltered
% data, while the blue signal cuts it down significantly. The amplitudes of the points in the red and blue signal are
% different throughout both waveforms (sometimes red is higher sometimes blue is higher). The blue waveform presents
% a significantly lower max voltage spike than the red waveform.  






%generally have a higher amplitude compared to the point in the red signal's waveform, 
% and the spikes in the red signal's waveform are usually higher. 
% 

%lower amplitude, shifted left to right (?)

%filt_filt output better 

%% 
% <latex>
%         \item Briefly explain the mathematical difference between the two filtering methods, and why one method might be more advantageous than the other in the context of spike detection? (5 pts)
% </latex>

%%

% The MATLAB filter function and filtfilt function differ in the way they
% filter input data. Specifically the MATLAB filter function represents a 1-D
% digital filter that uses a rational transfer function to filter the input
% data ("filter 1-D digital filter," 2020).
% However, the MATLAB filtfilt function does zero-phase digital filtering
% of the input data, and it filters the input data in the forward
% direction, reverses this now filtered data and puts it back through the
% filter ("filtfilt Zero phase digital filtering," 2020). Thus the filter created by the 
% filtfilt function also has zero phase distortion, a transfer function equal to the squared magnitude of the original transfer function, and an order that is double 
% the order specified by the transfer function numerator and demoniator
% coefficients inputted ("filtfilt Zero phase digital filtering,"2020). The MATLAB
% filtfilt method might
% be more advatangeous compared to the MATLAB filter method in the context of
% spike detection because the filtfilt method does have zero phase
% distortion, which means the filtfilt method will preserve the features of
% the unfiltered data exactly where they occur in the filtered waveform
% ("filtfilt Zero phase digital filtering,"2020). However, the filter function phase distorts the input data in the filtered output. This is
% important for spike detection because the time at which a spike occurs
% may yield insight on the behavior it is decribing, the cause of the spike, and/or affect its
% detection. As such, it would be best to use the filtfilt method so the input data
% is not phase-distorted when it is filtered in the context of spike detection tasks. 

%References:
%filtfilt Zero phase digital filtering. (2020). Retrieved March 10, 2020, from https://www.mathworks.com/help/signal/ref/filtfilt.html

%filter 1-D digital filter. (2020). Retrieved March 10, 2020, from https://www.mathworks.com/help/matlab/ref/filter.html

%%
% <latex>
%       \end{enumerate}
%         \item Using a spike threshold of +30 mV, calculate the index and value of the peak voltage for each spike in the \textbf{filtered} nerve signal (select the best one). Use these values to plot the first 2.5 seconds of the nerve signal with a red dot above (e.g. 10 mV above) each spike. (Hint: Plot the entire length of the nerve signal with all the spikes marked but then restrict the x-axis using \verb|xlim| to [0, 2.5] seconds) (4 pts)
% </latex>

%%

%run filt_filt on the whole dataset. filt_filt function is better than filter
%function

filt_filt_tot=filtfilt(b,a,(allData_nerve/1000));

%plotting 2.5 seconds of the filtered data 

%create 2.5 second time vector 
time_2=1/sample_rate:1/sample_rate:durationInSec+1/sample_rate;

%time_2=time_2*1000;

threshold = 30;

[pks, locs]=findpeaks(filt_filt_tot, time_2, 'MinPeakHeight', threshold)

%convert locatins to times in milliseconds
%loc_times=time_2(locs);
%%
% figure 
% plot((allData_nerve)/1000)
% figure 
% plot(filt_filt_tot)
%%

%create plot
figure
plot(time_2,filt_filt_tot)
hold on
for i=1:length(locs)
    plot(locs(i),pks(i)+10, 'r.','MarkerSize',10)
    hold on
end
xlim([0 2.5])
ylabel('Voltage (mV)')
xlabel('Time (Seconds)')
title('Filtered Nerve Signal With Peaks Plotted')
legend('Filtered Nerve Signal', 'Peak', 'Location', 'northwest')

%%

% The filtfilt method was used to generated the plots above (the filtfilt method is the best one).  

%% 
% <latex>
%  \item Under the assumption that different cells produce different action potentials with distinct peak amplitudes, decide how many cells you think were recorded (some number between 1 and 6). You may find it helpful to zoom in and pan on the plot you made in question 1.3. You may also find it useful to plot the sorted peak values to gain insight into where ``plateaus'' might be. (No need to include these preliminary plots in the report, though.) Use thresholds (which you well set manually/by eye) to separate the different spikes. Make a plot of the first 2.5 seconds similar to that in 1.3 except now color the spike dots of each group a different color (e.g., \verb|'r.'|,\verb|'g.'|,\verb|'k.'|,\verb|'m.'|).(6 pts)
% </latex>

%%

%plot to visualize thresholds. no need to include it. 
% figure
% yline(30)
% hold on
% yline(55)
% yline(100)
% yline(170)
% plot(time_2,filt_filt_tot)



threshold_1=30;
threshold_2=55;
threshold_3=100;
threshold_4=170;

[pks_1, locs_1]=findpeaks(filt_filt_tot, time_2, 'MinPeakHeight', threshold_1);
[pks_2, locs_2]=findpeaks(filt_filt_tot, time_2, 'MinPeakHeight', threshold_2);
[pks_3, locs_3]=findpeaks(filt_filt_tot, time_2, 'MinPeakHeight', threshold_3);
[pks_4, locs_4]=findpeaks(filt_filt_tot, time_2, 'MinPeakHeight', threshold_4);

figure
%Legend=cell(2,1);
for i=1:length(locs_1)
    h=plot(locs_1(i),pks_1(i)+10,'r.','MarkerSize',10);
    %Legend{1}=strcat('Cell 1 Peaks',num2str(1));
    hf(1)=h(1);
    hold on
 
end
for i=1:length(locs_2)
    k=plot(locs_2(i),pks_2(i)+10, 'g.','MarkerSize',10);
    %Legend{2}=strcat('Cell 2 Peaks',num2str(2));
    hf(2)=k(1);
    hold on
end
for i=1:length(locs_3)
    t=plot(locs_3(i),pks_3(i)+10, 'k.','MarkerSize',10);
    hf(3)=t(1);
    hold on
end
for i=1:length(locs_4)
    r=plot(locs_4(i),pks_4(i)+10, 'm.','MarkerSize',10);
    hf(4)=r(1);
    hold on
end

e=plot(time_2,filt_filt_tot);
hf(5)=e(1);
% %legend('Filtered Nerve Signal')
xlim([0 2.5])
ylabel('Voltage (mV)')
xlabel('Time (Seconds)')
title('Filtered Nerve Signal With Peaks Corresponding to the 4 Cells Plotted')
%legend(Legend)
%legend(p(1),'Cell 1 Peaks', h(1), 'Cell 2 Peaks')
legend(hf,{'Cell 1 Peaks','Cell 2 Peaks','Cell 3 Peaks','Cell 4 Peaks','Filtered Nerve Signal'}, 'Location', 'southoutside')

%% 
% <latex>
%  \item Use Matlab's $k$-means\footnote{Clustering, like $k$-means you are using here, is a form of unsupervised learning.} function (\verb|kmeans|) to fit $k$ clusters (where $k$ is the number of cells you think the recording is picking up) to the 1D data for each spike. 
%   \begin{enumerate}
% 	\item Using the same color order (for increasing spike amplitude) as you did for the thresholds in question 1.4, plot the spike cluster colors as little dots slightly above those you made for question 1.4. The final figure should be a new plot of the nerve voltage and two dots above each spike, the first being your manual label and the second your clustered label, which (hopefully/usually) should be the same color. (4 pts)
% </latex>

%%
%concate all peaks into one vector
peaks_vector=[];

peaks_vector=vertcat(pks_1,pks_2, pks_3,pks_4);

%%

%run matlab's k-means clustering
rng(1)% set the seed
idx=kmeans(peaks_vector,4);% this will group all the data into four groups.

all_locs=vertcat(locs_1',locs_2',locs_3',locs_4');

%create a mtrix with the flag, location, and peak for plotting
k_means_mat=[idx all_locs peaks_vector ];


%%

%create plot
figure
%plot from 1.4 
for i=1:length(locs_1)
    h=plot(locs_1(i),pks_1(i)+10,'r.','MarkerSize',10);
    %Legend{1}=strcat('Cell 1 Peaks',num2str(1));
    
    hf(1)=h(1);
    hold on
 
end
for i=1:length(locs_2)
    k=plot(locs_2(i),pks_2(i)+10, 'g.','MarkerSize',10);
    %Legend{2}=strcat('Cell 2 Peaks',num2str(2));
    hf(2)=k(1);
    hold on
end
for i=1:length(locs_3)
    t=plot(locs_3(i),pks_3(i)+10, 'k.','MarkerSize',10);
    hf(3)=t(1);
    hold on
end
for i=1:length(locs_4)
    r=plot(locs_4(i),pks_4(i)+10, 'm.','MarkerSize',10);
    hf(4)=r(1);
    hold on
end

%plot peaks from k means
for i=1:length(k_means_mat)
    if k_means_mat(i,1)==3
        plot(k_means_mat(i,2),k_means_mat(i,3)+30,'r.','MarkerSize',10);
        hold on
    elseif k_means_mat(i,1) ==4
        plot(k_means_mat(i,2),k_means_mat(i,3)+30, 'g.','MarkerSize',10);
        hold on
    elseif k_means_mat(i,1) == 1
        plot(k_means_mat(i,2),k_means_mat(i,3)+30, 'k.','MarkerSize',10);
        hold on
    elseif k_means_mat(i,1) ==2
        plot(k_means_mat(i,2),k_means_mat(i,3)+30, 'm.','MarkerSize',10);
        hold on
    end
end

e=plot(time_2,filt_filt_tot);
hf(5)=e(1);

xlim([0 2.5])
ylabel('Voltage (mV)')
xlabel('Time (Seconds)')
title({'Filtered Nerve Signal With Peaks Corresponding to the 4 Cells Plotted as Dots.', 'Kmeans Peaks Plotted Above Manual Threshold'})
legend(hf,{'Cell 1','Cell 2','Cell 3 Peaks','Cell 4','Filtered Nerve Signal'}, 'Location', 'southoutside')

%% 
% <latex>
% 	\item Which labeling, your manual ones or the ones learned by clustering) seem best, or do they both seem just as good? (Again, panning over the entire plot may be helpful.) (2 pts)
% </latex>


% As can be seen from the dots in the plot generated in part 1.5a, the clustering
% labels
% are almost the same as the dots created by manual plotting most all of the points. We can
% see however there are some points where the two labels are different (for example at x = 2.05 seconds). 
% While the clustering and manual labels are almost the same (and thus the performance of the two methods is almost the same), it appears that the manual labeling is better at doing the
% spike sorting (assigning the spikes to the correct cell) because the
% kmeans clustering unsupervised learning algorithim can at best provide an
% approximation of the sorting that the human eye can do. As such kmeans will
% not be as good as a human manually sorting the different spikes with
% thresholding. Additionally, manual sorting will not utilize any of the
% assumptions about the input data that are inherent in the k-means
% clustering algorithm, making it a more rigorous way to do labeling. 


%% 
% <latex>
%   \end{enumerate}
%  \item In this question,  you will test the hypothesis that the muscle potential responses are really only due to spikes from a subset of the cells you have identified in the previous two questions. First, plot the first 2.5 seconds of the muscle fiber potential and compare it with that of the nerve. Observe the relationship between spikes and the muscle fiber response. (No need to include this plot and observation in your report.)
%      Now, calculate the maximum muscle fiber potential change\footnote{max voltage - min voltage} in the 25 ms\footnote{Note that this 25 ms window is somewhat ad hoc and is just what seems reasonable by eye for this data. It implies no underlying physiological time scale or standard.} window after each spike (with the assumption that spikes without any/much effect on the muscle fiber potential do not directly innervate it). 
%   \begin{enumerate}
%    \item Using the cell groups you either manually defined or found via $k$-means clustering (just specify which you're using) again with different colors, plot a colored point for each spike where the x-value is the spike amplitude and the y-value is the muscle potential change. (6 pts)
% </latex>

%%

%use the manual groups

%first plot 2.5 seconds of the muscle response and compare it with the
%nerve
% figure
% plot(time_2, allData_nerve/1000)
% xlim([0 2.5])
% figure
% plot(time_2, allData_muscle/1000)
% xlim([0 2.5])

mv_muscle_data=allData_muscle/1000;
%define all the peak times
all_peak_times=k_means_mat(:,2);

max_potential_change=[];%list of potential changes

%now calculate the maximum muscle fiber potential change (max voltage - min voltage) in the 25 ms
%window after each spike 
for j=1:length(all_peak_times)
    %pull out 25ms = 0.025 seconds of muscle data. sample rate is 2000 Hz. so 0.025 sec of data
    %is 50 data points. 
    index_1=find(time_2==all_peak_times(j));%define the first index
    voltage_window=mv_muscle_data(index_1:index_1+49);%define the voltage window
    max_potential_change(j)=max(voltage_window)-min(voltage_window);%calculate max potential change
end
%     
% %create plot with the colors using the k means method 
% figure
% for k=1:length(max_potential_change)
%     if k_means_mat(k,1)==3
%         plot(k_means_mat(k,3),max_potential_change(k),'r.','MarkerSize',10);
%         hold on
%     elseif k_means_mat(k,1) ==4
%         plot(k_means_mat(k,3),max_potential_change(k), 'g.','MarkerSize',10);
%         hold on
%     elseif k_means_mat(k,1) == 1
%         plot(k_means_mat(k,3),max_potential_change(k), 'k.','MarkerSize',10);
%         hold on
%     elseif k_means_mat(k,1) ==2
%         plot(k_means_mat(k,3),max_potential_change(k), 'm.','MarkerSize',10);
%         hold on
%     end
% end
% %plot(k_means_mat(:,3),max_potential_change,'.','MarkerSize',10)
% xlabel('Peak Amplitude (mV)')
% ylabel('Max Muscle Potential Change (mV)')
% title('Max Muscle Potential Change v. Peak Amplitude Sorted Using K-Means Algorithm')

%create plot with the colors using the manual sorting thresholding method 
figure
for i=1:length(locs_1)
    h=plot(pks_1(i),max_potential_change(i),'r.','MarkerSize',10);
    %Legend{1}=strcat('Cell 1 Peaks',num2str(1));
    gf(1)=h(1);
    hold on
 
end
for i=1:length(locs_2)
    k=plot(pks_2(i),max_potential_change(i+155), 'g.','MarkerSize',10);
    %Legend{2}=strcat('Cell 2 Peaks',num2str(2));
    gf(2)=k(1);
    hold on
end
for i=1:length(locs_3)
    t=plot(pks_3(i),max_potential_change(i+268), 'k.','MarkerSize',10);
    gf(3)=t(1);
    hold on
end
for i=1:length(locs_4)
    r=plot(pks_4(i),max_potential_change(i+298), 'm.','MarkerSize',10);
    gf(4)=r(1);
    hold on
end
%plot(k_means_mat(:,3),max_potential_change,'.','MarkerSize',10)
xlabel('Peak Amplitude (mV)')
ylabel('Max Muscle Potential Change (mV)')
title('Max Muscle Potential Change v. Peak Amplitude Sorted Using Manually Defined Thresholds')
legend(gf,{'Cell 1','Cell 2','Cell 3','Cell 4'})



%% 
% <latex>
%    \item Does this plot support the hypothesis that the muscle fiber responses are only due to a subset of the cells. Explain why or why not. (3 pts)
% </latex>
%%

% No, this plot does not support the hypothesis that the muscle fiber responses
% are only due to a subset of cells. This is because we can clearly see
% that each cell (as grouped by manual sorting) creates a
% meaningful change in muscle potential compared to each of the other respective cells.
% Specifically cell 1 creates a change in muscle potential up to 9.6 mV,
% cell 2 creates a change in muscle potential up to 9.5 mV, cell 3 creates a
% change in muscle potential up to 9.4 mV, and cell 4 creates a change in
% muscle potential up to 6.8 mV. Thus, there were 4 cells recorded in this
% experiment and all of them contributed to a relatively large change in
% muscle potential and as such this plot does not support the hypothesis
% that the muscle fiber repsonses are only due to a subset of cells. The plot supports the 
% hypothesis that the muscle fiber responses are due to all cells. 

%% 
% <latex>
%   \end{enumerate}
% \end{enumerate}
% \section{Multivariate Clustering (21 pts)}
% In this section, you will explore similar methods for spikes sorting and clustering but with a different dataset, the human intracranial data in \verb|I521_A0006_D002|, 
% which is a larger dataset of the same recording you saw in \verb|I521_A0001_D001| of Homework 1. 
%   \begin{enumerate}
%    \item Using a threshold six standard deviations above the mean of the signal, detect the spikes in the signal. In addition, extract the waveform from 1 ms before the peak to 1 ms after it with peak value in the middle. (You will end up with a matrix where each row corresponds to the number of data points in 2 ms of signal minus 1 data point. Use the closest integer number of data points for the $\pm$ 1 ms window.) 
% </latex>

%%

%create new session


session2 = IEEGSession('I521_A0006_D002', 'andrewc', '/Users/andrewclark/Downloads/ieeg_password.bin' );
%%
%session2.data;

sample_rate2=32258;
%%
nr2 = ceil((session2.data.rawChannels(1).get_tsdetails.getEndTime)/1e6*session2.data.sampleRate);
%%
allData= session2.data.getvalues(1:nr2,1);

%%
durationInUSec_2 = session2.data(1).rawChannels(1).get_tsdetails.getDuration;
durationInSec_2 = durationInUSec_2 / 1e6;%178.8449 seconds

%%

%figure 
%plot(allData)
%allData_mv=allData/1000;



mean_data=mean(allData);
std_data=std(allData);
threshold_p2 = mean_data + 6*std_data;

%run find peaks function to detect spikes in the signal
[pks_2, locs_2]=findpeaks(allData);
locs_2=locs_2(pks_2>threshold_p2);
pks_2=pks_2(pks_2>threshold_p2);


%[pks_p2, locs_p2]=findpeaks(allData, 'MinPeakHeight', threshold_p2);

%sampling rate is 32258 Hz, so 1 ms of data is 1e-3* 32258= 32.2580, rounds to 32
%datapoints of data

%extract the waveform from 1ms before the peak to 1ms after peak
wave_form_mat=zeros(length(locs_2),65);

for u=1:length(locs_2)
    index_left=locs_2(u)-32;
    index_right=locs_2(u)+32;
    wave_form_mat(u,:)=allData(index_left:index_right);
    %index_1=find(time_2==all_peak_times(j));%define the first index
    %voltage_window=mv_muscle_data(index_1:index_1+49);%define the voltage
    %window
end

%% 
% <latex>
% 	\begin{enumerate}
% 	  \item Plot the waveforms of all the spikes overlaid on each other in the same color. (4 pts)
% </latex>

%%

%create 2.5 second time vector 
%time_2=1/sample_rate:1/sample_rate:durationInSec+1/sample_rate;

time_plot=(1/65:(2/65):2);

%create plot
figure
for s=1:308
    plot(time_plot,wave_form_mat(s,:),'-b')
    hold on
end
xlabel('Time (ms)')
ylabel('Voltage (uV)')
title('Waveforms of All Spikes Overlaid On Each Other')


%% 
% <latex>
% 	  \item Does it looks like there is more than one type of spike? (1 pt)
% </latex>
%%

% Yes it looks like there are two different types of spikes, one type has a higher amplitude than the other. 

%% 
% <latex>
% 	\end{enumerate} 
%    \item For each spike, represent the waveform by its  principal components. Use the \verb|pca| command in Matlab. Intuitively, principal component analysis finds the coordinate system that most reduces the variability in your data. 
% 	\begin{enumerate}
% 	  \item Run principal component analysis on all the spike waveforms and represent your data with the top two principal components. Make a scatterplot of your data in this principal component (PC) space. (3 pts)
% </latex>

%%

% Run PCA on all spike waveforms and represent data with top two principal
% components. 

transposed_mat=wave_form_mat;

[coeff,score,latent,~,explained]=pca(transposed_mat);

%Represent data in the PC space
%PC_space=transposed_mat*coeff(:,1:2);

%%

%create plot
figure
plot(score(:,1),score(:,2),'.','MarkerSize',14)
xlabel('Principal Component 1 (uV)')
ylabel('Principal Component 2 (uV)')
title('Top Two Principal Components of the Data')
% figure 
% bar(latent)

%% 
% <latex>
% 	  \item Each PC also has an associated eigenvalue, representing the amount of variance explained by that PC. This an output of the \verb|PCA| command. Plot the  principal component vs the total (cumulative) percent variance explained. What is the percent variance explained if you include the top two principal components? (3 pts)
% </latex>
%%

explained_plotted=[];
for i=1:length(explained)
    value=explained(1:i);
    explained_plotted(i)=sum(value);
end

figure
%plot(explained,score,'b.','MarkerSize',10)
plot(explained_plotted)
xlabel('Number of Principal Components')
ylabel('Percent Variance Explained (%)')
top_two_explained=explained(1)+explained(2);
title('Principal Component v. Total Percent Variance Explained')
ylim([60 105])

%%

% The top two principal components explain 73.7349% of the total variance.

%% 
% <latex>
% 	  \item Does it look like there is more than one cluster of spikes? (1 pt)
% 	\end{enumerate} 
% </latex>
%%

%Yes, it looks like there are more than one cluster of spikes. It looks like there are two clusters of spikes. 

%% 
% <latex>
%    \item Use the same \verb|kmeans| function as you used before to cluster the spikes based on these two (normalized) features (the waveforms represented by the top two PCs). You will use a slight twist, though, in that you will perform $k$-medians (which uses the medians instead of the mean for the cluster centers) by using the \verb|'cityblock'| distance metric (instead of the default \verb|'sqEuclidean'| distance). Make a plot similar to that in 2.2.a but now coloring the two clusters red and green. (3 pts)
% </latex>
%%

%run k-medians

idx_2=kmeans(score(:,1:2),2,'Distance','cityblock');

%create plot
figure
for i=1:length(score)
    if idx_2(i)==1
        q=plot(score(i,1),score(i,2),'r.','MarkerSize',14);
        qf(1)=q;
        hold on
    elseif idx_2(i)==2
        o=plot(score(i,1),score(i,2),'g.','MarkerSize',14);
        qf(2)=o;
        hold on
    end
end
xlabel('Principal Component 1 (uV)')
ylabel('Principal Component 2 (uV)')
title('Top Two Principal Components of the Spike Data Sorted by K-Medians')
legend(qf,{'Spike Cluster 1','Spike Cluster 2'})


%% 
% <latex>
%   \item Make a plot similar to 2.1 but now coloring the traces red and green according to which cluster they are in. Overlay the mean of the waveforms in each cluster with a thick black line (use the parameter \verb|'LineWidth'| and value \verb|'4'|). (3 pts)
% </latex>

%%

%create plot with colored traces. 

%figure out how many 1s and 2s there are in the idx_2
count_1s=0;
count_2s=0;
for i=1:length(idx_2)
    if idx_2(i)==1
        count_1s=count_1s+1;
    elseif idx_2(i)==2
        count_2s=count_2s+1;
    end
    
end
%%
figure
cluster_1=zeros(count_1s,65);
cluster_2=zeros(count_2s,65);
for s=1:308
    if idx_2(s)==1
        e=plot(time_plot,wave_form_mat(s,:),'-r');
        cluster_1(s,:)=wave_form_mat(s,:);
        kl(1)=e;
        hold on
    elseif idx_2(s)==2
        z=plot(time_plot,wave_form_mat(s,:),'-g');
        %cluster_2=[cluster_2 wave_form_mat(s,:)];
        cluster_2(s,:)=wave_form_mat(s,:);
        kl(2)=z;
        hold on
    end
end

cluster_1_filtered=cluster_1(any(cluster_1,2),:);
cluster_2_filtered=cluster_2(any(cluster_2,2),:);
mean_cluster_1=mean(cluster_1_filtered,1);
hold on
mean_cluster_2=mean(cluster_2_filtered,1);
hold on
q=plot(time_plot,mean_cluster_1,'k-','LineWidth',4);
kl(3)=q;
hold on
plot(time_plot,mean_cluster_2,'k-','LineWidth',4)
xlabel('Time (ms)')
ylabel('Voltage (uV)')
title('Waveforms of All Spikes Overlaid On Each Other With Means Of Each Shown')
legend(kl,{'Spike Cluster 1','Spike Cluster 2','Means Of Each Cluster'})

%% 
% <latex>
%   \item What are some dangers of using the clustering techniques in this homework? (List 3) (3 pts)
% </latex>
%%

% Some of the dangers of using the clustering techniques in this homework
% include: First, k-means clustering assumes that each cluster present in
% the input data should have roughly the same about of datapoints in it
% (Robinson, 2015). This is a danger because if the clusters actually have different numbers of 
% datapoints the k-means clustering unsupervised learning algorithm may
% group the data non-optimally, leading to inaccurate
% conclusions. Second, the behavior of the k-means clustering algorithm is also
% greatly effected by the presence of outliers in the input data, when
% outliers are present the output clusters may be created non-optimally
% leading to a poor representation of the data (Franklin, 2019). 
% A third danger that the input data to the k-means clustering technique used in this
% homework should actually have some inherehent clusters in it or else
% k-means clustering will not applicable; that is if k-means clustering is
% run on a uniform data (non-clustered) dataset the k-means clustering algorithm will partition
% the data into clusters, even though the clusters would be meaningless in
% this case (Kim, 2015).



%References:
%Franklin, J. (2019, November 16). Effect of outliers on K-Means algorithm using Python. Retrieved March 10, 2020, from https://medium.com/analytics-vidhya/effect-of-outliers-on-k-means-algorithm-using-python-7ba85821ea23

%Kim, K. (2015, February). How to understand the drawbacks of K-means. Retrieved March 10, 2020, from https://stats.stackexchange.com/questions/133656/how-to-understand-the-drawbacks-of-k-means

%Robinson, D. (2015, January 16). K-means clustering is not a free lunch . Retrieved March 10, 2020, from http://varianceexplained.org/r/kmeans-free-lunch/

%% 
% <latex> 
% \end{enumerate}
% \end{document}
% </latex>
