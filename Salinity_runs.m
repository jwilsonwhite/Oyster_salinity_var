

name = 'Salinity_stdev_multiplied';

Y = [0 0.05 0.1];
for i= 1:length(Y)
    for run = 1:2000
        TS_sal = Create_Mock_Salinity(Y(i));   
    
    save_tmp = strcat('Salinity_ts/',name,num2str(Y(i)),'_run_num_',num2str(run),'.mat');
 %   save(save_tmp,"TS_sal")

    end

end