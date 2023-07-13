function [cell_type] = func_get_cell_type_WhisperSwitch(waveform, whisper_flag)


if whisper_flag == 0
    % old recording system (Guo & Li 2014 criteria), sampling rate=19531
    waveform_tmp = mean(waveform);
    waveform_tmp = waveform_tmp/norm(waveform_tmp);
    [wave_min i_peak_min] = min(waveform_tmp);
    [wave_max i_peak_max] = max(waveform_tmp(i_peak_min:end));
    spk_width_tmp = i_peak_max;
    
    if spk_width_tmp > 9
        cell_type = 1;
    elseif spk_width_tmp < 7
        cell_type = 2;
    else
        cell_type = 0;
    end


elseif whisper_flag==1

    % whisper system, switched on 3/1/15 (Guo & Li 2014 criteria), sampling rate=25000
    waveform_tmp = mean(waveform);
    waveform_tmp = waveform_tmp/norm(waveform_tmp);
    [wave_min i_peak_min] = min(waveform_tmp);
    [wave_max i_peak_max] = max(waveform_tmp(i_peak_min:end));
    spk_width_tmp = i_peak_max;
    
    if spk_width_tmp > round(9/19531*25000)
        cell_type = 1;
    elseif spk_width_tmp < round(7/19531*25000)
        cell_type = 2;
    else
        cell_type = 0;
    end
    
else
    error('')
end