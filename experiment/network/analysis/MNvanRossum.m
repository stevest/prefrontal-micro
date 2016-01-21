% MultiNeuron van Rossum:

function output = MNvanRossum(trains_a,trains_b,tau,c)
% input spike trains are cells
    output = [];
    output = d(trains_a,trains_b,tau,c);
    return;
end

function out = d(trains_a,trains_b,tau,c)
out = [];

big_p = size(trains_a,1);

d=0;

%R_p
for p=1:big_p
    d = d + big_r(trains_a(p),tau);
    d = d + big_r(trains_b(p),tau);
    d = d - (2*big_r2(trains_a(p),trains_b(p),tau));
end

r_pq=0;

for p=1:big_p-1
    for q=p+1:big_p
        r_pq = r_pq + big_r2(trains_a(p),trains_a(q),tau) ;
        r_pq = r_pq + big_r2(trains_b(p),trains_b(q),tau);
        r_pq = r_pq - big_r2(trains_a(p),trains_b(q),tau);
        r_pq = r_pq - big_r2(trains_b(p),trains_a(q),tau);
    end
end

d = d + (2*c*r_pq);

out = sqrt(d);
return ;

end





function total = big_r(u,tau)

u_size = length(u{:});

total = [];
if(u_size==0)
    return;
end

total = u_size;

for i=1:u_size
    for j=i+1:u_size
        total = total + (2*exp(-(u{1}(j)-u{1}(i))/tau) ) ;
    end
end

return;

end




function total = big_r2(u_a,u_b,tau)

u_a_size = length(u_a{:});
u_b_size = length(u_b{:});

total = [];
if(u_a_size==0 || u_b_size==0)
    return;
end

total=0;

for i=1:u_a_size
    for j=1:u_b_size
        total = total + exp(-abs(u_a{1}(i)-u_b{1}(j))/tau);
    end
end

return;

end