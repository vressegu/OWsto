function w = fct_mirror_on_border(w,nbp)
% Recplicate value ouside the domain assuming double periodicity
%

w=w(:,:,1)+1i*w(:,:,2);
w=[w(end-nbp+1:end,:);w;w(1:nbp,:)];
w=[repmat(w(:,1),[1 nbp]) w repmat(w(:,end),[1 nbp])];
% w=reshape(w,[(MX(2)+2*nbp)*(MX(1)+2*nbp)]);
% w=[w;w;w];
% w=[w(:,end-Py+1:end) w w(:,1:Py)];
% % w=reshape(w,[(MX(2)+2*Py)*3*MX(1) ]);
% w=w(:);
w(:,:,2)=imag(w);
w(:,:,1)=real(w(:,:,1));
