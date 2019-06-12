function temp=fct_real_sst_erwan(model)
% Filter the SQG temperature field
% 

load(['/Users/vressegu/Documents/manu/matlab_files/erwan1sqg' ...
    num2str(model.grid.MX(1)) 'all.mat'],'temp');
temp=temp(:,:,model.id_t);
temp=temp-mean(mean(temp));
