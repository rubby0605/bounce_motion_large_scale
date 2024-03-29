        
        td = 1; 
        fd = 1;
        ht = 1;
        nn = count_num(td,fd,ht);
        if(nn < 8),disp('!'),end
        
        max_v = max_v_mat(td, fd, ht);
        min_v = min_v_mat(td, fd, ht);
        p2(1:num_poly+1) = poly_mat(td,fd,ht,1:num_poly+1);
        xx1 = linspace(min_v, max_v, ceil(nn/2));
        vv2 = polyval(p2, xx1);
        plot(xx1, vv2, 'kp');hold on;
        vv3=zeros(1e4, 1);
  %%      用二分法  應該改成牛頓法!!
        for nn2 = 1 : 1e3
            df = polyval(p2, max_v);
            ranking = rand * (nn-df) + df;
            hold on;
           % plot(xx1, ones(size(xx1)).*ranking,'r-');
            hold on;
            
            f1 =  polyval(p2, max_v) - ranking;
            v1 = max_v;
            f0 = polyval(p2, min_v) - ranking;
            v0 = min_v;
                for num_dev = 1 : 1e2
                    v2 = (v1 + v0)/2;
                    f2 = polyval(p2, v2) - ranking;
                    if(f2 * f1 < 0),f0 = f2; v0 = v2;
                    else f1 = f2; v1 = v2;
                    end
                    if(abs(v1-v0) < 0.01)
                        break;
                    end
                end
              %  v2
            vv3(nn2)=v2;
            plot(v2, ranking,'r');hold on;
        end
        %
        h2 = hist(vv3(1:1000),ceil(nn/2));
        h3 = zeros(ceil(nn/2),1);
        for nn3 = 1 : ceil(nn/2)
            h3(nn3) = sum(h2(nn3:ceil(nn/2)));
        end
        %plot(xx1, h3./max(h3).*50,'g-');