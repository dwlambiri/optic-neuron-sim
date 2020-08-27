function lines = test(m, main_r)
    [vv, cc] = voronoin([2*main_r*(rand(1,m)-0.5)' 2*main_r*(rand(1,m)-0.5)']);
    
    total = 0;
    for q=1:length(cc)
        total = total + length(cc{q})-1;
    end
    
    lines = zeros(4,total);
    good = vv < main_r & vv > -main_r;
    good = good(:,1) & good(:,2);
    iter = 1;
    for q=1:length(cc)
        face = cc{q};
        for w=1:length(face)-1
           if good(face(w)) == 0 || good(face(w+1))==0
               continue;
           end 
           lines(1,iter) = vv(face(w),1);
           lines(2,iter) = vv(face(w),2);
           lines(3,iter) = vv(face(w+1),1);
           lines(4,iter) = vv(face(w+1),2);
           iter = iter+1;
        end
    end

    rem = 2;
    
    num = iter;
    obstacle_r = 20*ones(1,num);
    
    fprintf("Number of obstacles = %d\n", num);
end