function Xl_m = nnMpEvenRotationMatrix(euler_m, l)
% Usage: Xl_m = nnMpEvenRotationMatrix(euler_m, l)
%
% nnMpEvenRotationMatrix calculates the rotation matrix that transforms
% signal coefficients (real valued and in nnMpEven format) to the frames 
% described by the euler angles in each row of euler_m. If the euler angles
% give the orientation of the lab frame relative to each FF, then Xl_m gives
% the transformation from each fiber-frame to the lab-frame [See 
% Eqs. (B.20) and (B.30) in the 'Response Function' Notability document 
% for details].
%
% The input parameters are:
%
%   euler_m: Nf by 3 matrix of Euler angles ([alpha, beta, gamma] in 
%       one row per fiber-frame) giving the orientation of the lab frame,
%       relative to the fiber frame, or more generally from the initial
%       frame to the frame described by the Euler angles.
%   l: order of signal coefficients.
%
% Output variables are:
%
%   Xl_m: (2l+1) by (2l+1)*Nf matrix that transforms fiber-frame signal 
%       coefficients (in real-valued, nnMpEven format) to the lab-frame
%       (same format), or more generally from the original frame to the
%       frames specified by the Euler angles in each row of euler_m.
%
% nnMpEvenRotationMatrix calls the functions:
%   D_v = wignerD2(l, mPrime, m, alpha_v, beta_v, gamma_v)

% Derived parameters:
Nf = size(euler_m, 1);

idx = 1:(2*l+1);
isEven = mod(idx,2)==0;
isOdd = ~isEven;

% Loop over fiber-frames:
Xl_m = [];
for fIndex = 1:Nf
    alpha = euler_m(fIndex, 1);
    beta = euler_m(fIndex, 2);
    gamma = euler_m(fIndex, 3);
    % Calculate Rl_m for current fiber-frame:
    Rl_m = zeros(2*l+1, 2*l+1);
    % Loop over matrix rows:
    for row = 1:(2*l+1)
        % Find current m value:
        m = ceil((row-1)/2);
        % Loop over matrix columns:
        for col = 1:(2*l+1)
            % Find current mPrime value:
            mPrime = ceil((col-1)/2);
            % Find current Wigner rotation matrix elements:
            DlmPm = wignerD2(l, mPrime, m, alpha, beta, gamma);
            DlnmPm = wignerD2(l, -mPrime, m, alpha, beta, gamma); % -mPrime.
            % Calculate Rl_m:
            if (row == 1 && col == 1)
                % Eq. (B.21):
                Rl_m(1,1) = DlmPm;
            elseif (row == 1 && col > 1)    % For clarity.
                % Eq. (B.22):
                Rl_m(1, col) = 2*real(DlmPm)*isEven(col) ...
                    + 2*imag(DlmPm)*isOdd(col);
            elseif (row > 1 && isEven(row) && col == 1) % For clarity.
                % Eq. (B.25),  simplified:
                Rl_m(row, 1) = real(DlmPm);
            elseif (row > 1 && isEven(row) && col > 1)  % For clarity.
                % Eq. (B.26), simplified:
                Rl_m(row, col) = real(DlmPm + ((-1)^mPrime)*DlnmPm)*isEven(col) ...
                    + imag(DlmPm - ((-1)^mPrime)*DlnmPm)*isOdd(col);
            elseif (row > 1 && isOdd(row) && col == 1)  % For clarity.
                % Eq. (B.27), simplified:
                Rl_m(row, 1) = -imag(DlmPm);
            elseif (row > 1 && isOdd(row) && col > 1)   % For clarity.
                % Eq. (B.28), simplified:
                Rl_m(row, col) = -imag(DlmPm + ((-1)^mPrime)*DlnmPm)*isEven(col) ...
                    + real(DlmPm - ((-1)^mPrime)*DlnmPm)*isOdd(col); 
            end % End of row,col IF/ELSEIF.                
        end     % End of col FOR.
    end         % End of row FOR.
    % Append Rl_m for current fiber-frame to Xl_m:
    Xl_m = [Xl_m, Rl_m];
end     % End of fIndex FOR.
