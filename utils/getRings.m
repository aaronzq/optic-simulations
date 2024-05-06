function rings=getRings(X1,Y1,radius,nRing)

rings = zeros(size(X1,1),size(X1,2),nRing);

for i=1:nRing
    if i==1
        rings(:,:,i)=(X1.^2+Y1.^2)<(i/nRing)*(radius)^2;
    else
        temp = sum(rings,3);
        rings(:,:,i)=(X1.^2+Y1.^2)<(i/nRing)*(radius)^2;      
        rings(:,:,i) = rings(:,:,i) & ~temp;
    end
end

end