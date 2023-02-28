% The code converts the Jsave and A to synaptic weight and connectivity matrices.

A=zeros(N);
J=zeros(N);

for m1 = 1:size(AB,1) % J Pre , X row
    for m2 = 1:N % J Post, X column
        if AB(m1,m2)
            A(m2,AB(m1,m2))=1;
            J(m2,AB(m1,m2)) = Jsave(m1,m2);
        end
    end
end

%%

