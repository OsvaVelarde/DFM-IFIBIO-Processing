function refsignal = function_referencing(signal,REFcfg)

    % Inputs:
    % signal
    % method
    % structElectrode

    % Outputs:
    % refsignal

    %Argument completion ------------------------------------------------------
    if ~isfield(REFcfg, 'electrode'),
        error('MATLAB:function_referencing','Error in the ELECTRODE configuration');
    end

    if ~isfield(REFcfg, 'method'),
        method = 'hardware';
    else
        method = REFcfg.method;
    end

    %----------------------------------------------------------------
    Electrode = REFcfg.electrode;
    NumberChannels = size(Electrode,1)*size(Electrode,2);
    
    switch lower(method), %Switch for Referencing selection.
        case 'car', %Implement the Common Average Referencing (CAR).

            car = mean(signal,2);                                               %Compute the CAR for the pre-configured electrodes.
            refsignal = signal - kron(ones(1,NumberChannels) ,car);             %Referencing the signals to the CAR.

            disp('Common Average Referencing');

        case 'br', %Implement the Bipolar Referencing (BR). 
            
            refsignal = signal;

            for ii = 1:size(Electrode,1)-1
                for jj = 1:size(Electrode,2)
                    refsignal(:,Electrode(ii,jj)) = signal(:,Electrode(ii,jj)) - signal(:,Electrode(ii+1,jj));
                end
            end
            
            disp(...
                'Bipolar Referencing: La señal se resta hacia abajo. El ultimo contacto está a Hardware.');
            
        otherwise,
            refsignal = signal;
            disp('Referencing is HARDWARE by default');

    end %Switch for Referencing selection.
                
end


