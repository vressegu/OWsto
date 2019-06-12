function delta_t = fct_choice_time_step(model,w)
% Compute the spectrum of a function
%

dx_mesure = 1e4;
fw=fft2(w);
U(:,:,1)= fct_U_k(fw(:,:,1),model,dx_mesure);
U(:,:,2)= fct_U_k(fw(:,:,2),model,dx_mesure);
U=sqrt(sum(U(:)^2)/prod(model.grid.MX))
delta_t=dx_mesure/U

%%
    function U = fct_U_k(ft,model,dx_m)
        
        MX=model.grid.MX;
        dX=model.grid.dX;
        
        if any(size(ft)~=MX)
            error('wrong size');
        end
        % [Mx,My]=size(ft);
        % MX=[Mx My];
        % M=floor(sqrt(Mx*My));
        % M=min(MX/8);
        M=min(MX/4);
        
        persistent idxref MXref dXref
        
        if ~exist('MXref','var') ||  isempty (MXref) || any(MXref ~= MX) ...
                || any(dXref ~= dX)
            MXref=MX;
            dXref=dX;
            %% Define wave number
            if any( mod(MX,2)~=0)
                error('the number of grid points by axis need to be even');
            end
            PX=MX/2;
            ft(PX(1),:)=0;
            ft(:,PX(2))=0;
            kx=1/(MX(1)*dX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
            ky=1/(MX(2)*dX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
            %     kx=1/Mx*[ 0:(PX(1)-1) 0 (1-PX(1)):-1];
            %     ky=1/My*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
            [kx,ky]=ndgrid(kx,ky);
            k=log(kx.^2+ky.^2)/log(10)/2; clear kx ky
            % k=sqrt(kx.^2+ky.^2); clear kx ky
            k=k(:);
            
            %% Time step choice
            % dx_m = 1e4;
            k_mesure = -log(2*dx_m)/log(10);
            % idx = sparse( bsxfun(@le,kidx(:,1:end-1), k_mesure ) );
            % idx = idx & sparse( bsxfun(@lt,k_mesure, kidx(:,2:end)  ) );
            idx = sparse( bsxfun(@le,k_mesure, k ) );
            keyboard;
            idx = idx & sparse( bsxfun(@lt,k, k_mesure) );
            idxref=reshape(idx,MX);
        end
        spectrum2 = idxref .* ft;
        
        contour(abs(spectrum2));
        keyboard;
        
        U=ifft2(spectrum2);
        
        norm_table(imag(U))/norm_table(real(U))
        keyboard;
        
        U=real(U);
    end

end