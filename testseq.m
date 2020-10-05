function testseq(noPlanes, simIterations)

totalPlanes = noPlanes + 2;
head = 2;
for i=0:simIterations
    %tic
    fprintf("LEVEL %d [%d, %d, %d, %d]\n", 1, 1+mod(head-2, totalPlanes), 1+mod(head, totalPlanes), 1+mod(head, totalPlanes), 1+mod(head+1, totalPlanes));    
    for j = 1: noPlanes-2
        fprintf("LEVEL %d [%d, %d, %d, %d %d]\n", j, 1+mod(head-2+j, totalPlanes), 1+mod(head-2+j, totalPlanes), 1+mod(head-1+j, totalPlanes), 1+mod(head+j, totalPlanes), 1+mod(head+j+1, totalPlanes));
    end
    fprintf("LEVEL %d [%d, %d, %d, %d %d]\n\n\n", noPlanes, 1+mod(head-2+noPlanes-1, totalPlanes), 1+mod(head-2+noPlanes-1, totalPlanes), 1+mod(head-1+noPlanes-1, totalPlanes), 1+mod(head+noPlanes-1, totalPlanes), 1+mod(head+noPlanes-1, totalPlanes));

    head = mod(head-2, totalPlanes);
end

end