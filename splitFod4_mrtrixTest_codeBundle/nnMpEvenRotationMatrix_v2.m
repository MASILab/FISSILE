function Xl_m = nnMpEvenRotationMatrix_v2(euler_m, l)

persistent fac
if isempty(fac)
    fac = factorial(0:20);
end

idx = 1:(2*l+1);
isEven = mod(idx,2)==0;
isOdd = ~isEven;

alpha = euler_m(1, 1);
beta = euler_m(1, 2);
gamma = euler_m(1, 3);
Xl_m = zeros(2*l+1, 2*l+1);

% Loop over matrix rows:
for row = 1:(2*l+1)
    % Find current m value:
    m = ceil((row-1)/2);
    % Loop over matrix columns:
    for col = 1:(2*l+1)
        % Find current mPrime value:
        mPrime = ceil((col-1)/2);
        % Find current Wigner rotation matrix elements:

        kappa_v = max(0, m-mPrime):min(l+m, l-mPrime);  % Must have at least one element (0).
        dSum_v = zeros(size(beta));
        term_v = zeros(size(beta));
        for kappaIndex = 1:length(kappa_v)
            kappa = kappa_v(kappaIndex);
            if kappaIndex == 1
                term_v = (-1)^kappa * cos(beta/2).^(2*l+m-mPrime-2*kappa) .* ...
                (-sin(beta/2)).^(mPrime-m+2*kappa) / ...
                (fac(l-mPrime-kappa+1)*fac(l+m-kappa+1)*...
                fac(kappa+mPrime-m+1)*fac(kappa+1));
            else
                term_v = term_v .* (-tan(beta/2).^2) .* ((l-mPrime-kappa+1) * ...
                (l+m-kappa+1) / (kappa+mPrime-m) / kappa);
            end
            dSum_v = dSum_v + term_v;
        end
        d_v = sqrt(fac(l+m+1)*fac(l-m+1)*fac(l+mPrime+1) * ...
            fac(l-mPrime+1)) * dSum_v;
        DlmPm = exp(-1i*mPrime*alpha) .* d_v .* exp(-1i*m*gamma);
        
        mPrime = -mPrime;
        kappa_v = max(0, m-mPrime):min(l+m, l-mPrime);  % Must have at least one element (0).
        dSum_v = zeros(size(beta));
        term_v = zeros(size(beta));
        for kappaIndex = 1:length(kappa_v)
            kappa = kappa_v(kappaIndex);
            if kappaIndex == 1
                term_v = (-1)^kappa * cos(beta/2).^(2*l+m-mPrime-2*kappa) .* ...
                (-sin(beta/2)).^(mPrime-m+2*kappa) / ...
                (fac(l-mPrime-kappa+1)*fac(l+m-kappa+1)*...
                fac(kappa+mPrime-m+1)*fac(kappa+1));
            else
                term_v = term_v .* (-tan(beta/2).^2) .* ((l-mPrime-kappa+1) * ...
                (l+m-kappa+1) / (kappa+mPrime-m) / kappa);
            end
            dSum_v = dSum_v + term_v;
        end
        d_v = sqrt(fac(l+m+1)*fac(l-m+1)*fac(l+mPrime+1) * ...
            fac(l-mPrime+1)) * dSum_v;
        DlnmPm = exp(-1i*mPrime*alpha) .* d_v .* exp(-1i*m*gamma);

        mPrime = -mPrime;
        if (row == 1 && col == 1)
            % Eq. (B.21):
            Xl_m(1, col) = real(DlmPm);
        elseif (row == 1 && col > 1)    % For clarity.
            % Eq. (B.22):
            Xl_m(1, col) = 2*real(DlmPm)*isEven(col) ...
                    + 2*imag(DlmPm)*isOdd(col);
        elseif (row > 1 && isEven(row) && col == 1) % For clarity.
            % Eq. (B.25),  simplified:
            Xl_m(row, col) = real(DlmPm);
        elseif (row > 1 && isEven(row) && col > 1)  % For clarity.
            % Eq. (B.26), simplified:
            Xl_m(row, col) = real(DlmPm + ((-1)^mPrime)*DlnmPm)*isEven(col) ...
                    + imag(DlmPm - ((-1)^mPrime)*DlnmPm)*isOdd(col);
        elseif (row > 1 && isOdd(row) && col == 1)  % For clarity.
            % Eq. (B.27), simplified:
            Xl_m(row, col) = -imag(DlmPm);
        elseif (row > 1 && isOdd(row) && col > 1)   % For clarity.
            % Eq. (B.28), simplified:
            Xl_m(row, col) = -imag(DlmPm + ((-1)^mPrime)*DlnmPm)*isEven(col) ...
                    + real(DlmPm - ((-1)^mPrime)*DlnmPm)*isOdd(col); 
        end % End of row,col IF/ELSEIF.           
    end
end