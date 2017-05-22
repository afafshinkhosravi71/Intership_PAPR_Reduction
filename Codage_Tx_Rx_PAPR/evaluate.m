function [mse, mae, SNR, PSNR] = evaluate (img1, img2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes metrics 
% Metrics mean squared error MSE, mean absolute error MAE, 
% signal to noise ratio SNR, peak signal to noise ratio PSNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


temp = double(img1);
y = double(img2);
[l, c, ch] = size(temp);

% MSE : Mean squared error
mse = 0;
for k = 1:ch
    for i = 1:l
        for j = 1:c
            mse = mse + (y(i,j,k)-temp(i,j,k))^2;
        end
    end
end
mse = mse/(l*c*ch);
fprintf('MSE: Mean Squared Error %f\n',mse);

% MAE: Mean absolute error
mae = 0;
for k = 1:ch
    for i = 1:l
        for j = 1:c
            mae = mae + abs(y(i,j,k)-temp(i,j,k));
        end
    end
end
mae = mae/(l*c*ch);
fprintf('MAE: Mean Absolute Error %f\n',mae);


%SNR and PSNR %signal to noise ratio %peak signal to noise ratio
num = 0;
den = 0;
for k = 1:ch
    for i = 1:l
        for j = 1:c
            den = den + (y(i,j,k)-temp(i,j,k))^2;
        end
    end
end
for k = 1:ch
    for i = 1:l
        for j = 1:c
            num = num + temp(i,j,k)^2;
        end
    end
end
SNR = 20*log10(sqrt(num)/sqrt(den));
PSNR = 20*log10(255/sqrt(mse));
fprintf('SNR: Signal to Noise Ratio %f dB \n',SNR);
fprintf('PSNR: Peak Signal to Noise Ratio %f dB \n',PSNR);