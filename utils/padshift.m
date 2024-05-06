function data_shifted = padshift(data, shift)
    % data (1,length)
    if shift >= 0
        data_shifted = double(ones(size(data))) * data(1);
        data_shifted(shift+1:end) = data(1:end-shift);
    else
        data_shifted = double(ones(size(data))) * data(end);
        data_shifted(1:end+shift) = data(-shift+1:end);
    end
end