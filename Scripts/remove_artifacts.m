function [index_artifacts, new_signal] = remove_artifacts(time,signal, window_time, method, plot_i)

% Inputs: 
%   * signal: samples x channels.
%   * window_time: positive value.
%   * method: 0: zero padding, 1: linear extension , 2:quadratic extension
%   * plot: 0: don't show graphical results, 1: show

new_signal=signal;
average_abs = mean(abs(signal),2);
mov_av = movmean(average_abs,window_time);
threshold = median(mov_av)/0.62; 
index_artifacts = mov_av > threshold;

switch method 
    case 0 %zero padding
        new_signal(index_artifacts,:)=0;
    case 1 %modify linear approx
        new_signal(index_artifacts,:)=0;
    case 2 %modify quadratic approx.
        new_signal(index_artifacts,:)=0;
end

if plot_i==1
    figure()
    
    subplot(4,1,1)
    plot(time,signal)
    
    subplot(4,1,2)
    plot(time,mov_av)
    
    subplot(4,1,3)
    plot(time,class)
    
    subplot(4,1,4)
    plot(time,new_signal)
end

end