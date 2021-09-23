    %fid = fopen('ttime_1919_0325_1e24.dat','r');
    fid = fopen('ttime_1e4_1e24_0417.dat','r');
    line = fgetl(fid);
    str = sscanf(line, '%s');
    ttime2 = zeros(1e4*8, 1);
    ii = 0;

    for num_ii = 1 : 8
        for num_jj = 1 : totalnum
            line = fgetl(fid);
            str = sscanf(line, '%f');
            ttime(num_ii, num_jj) = str;
        end
    end
    disp('!');
    line = fgetl(fid);
    fclose(fid)
    %% check
    
    %totaln = 1919;
    %totaln = 5179;
    ttime(1, totaln)
    %ttime(1, totaln+1:totaln+10)
    
    %% remove the tailer
    
    %ttime(:,totaln+1:1e4)=[];
    
    %%
    % the form ttime
    % //  
    %ttime(1,n) = start_thi;
    %ttime(2,n) = final_thi;
    %ttime(3,n) = time;   
    %ttime(4,n) = W;
    %ttime(5,n) = flight_time;
    %ttime(6,n) = final_fi;
    %ttime(7,n) = start_plotT(td,fd);
    %ttime(8,n) = rest_time;
    
    % ouou part
    table_state = zeros(4, 1);
    % [ouou photoionized 
    ouou = 0;
    for num = 1 : totaln
        if(ttime(2, num)==0),
            if(ttime(4, num)==0)
                ouou = ouou + 1;
            else
                ouou = ouou + ttime(4, num);
            end
        end
    end
    table_state(1) = ouou./totaln;
    W_mat = zeros(totaln, 1);
    W_mat(:) = ttime(4,:);
    %hist(W_mat)
    
    % photoionized part
    W_loss = 0;
    for num = 1 : totaln
        W_loss = W_loss + (1 - ttime(4, num));
    end
    disp(W_loss./totaln)
    table_state(2) = W_loss./totaln;
    
    % southern polar caps
    S_trap = 0;
    N_trap = 0;
    for num = 1 : totaln
        if(ttime(2, num) < 0)
            %southern polar caps
            S_trap = S_trap + ttime(4, num);
        else
            N_trap = N_trap + ttime(4, num);
        end
    end
    disp(N_trap./totaln)
    disp(S_trap./totaln)
    table_state(4) = S_trap./totaln;
    table_state(3) = N_trap./totaln;
    
    disp((table_state(:))./sum(table_state))
    
    
