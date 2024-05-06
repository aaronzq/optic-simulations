function fout = rotation_rep(fin)
n = length(fin);
[rx,ry]=meshgrid( -(n-1)/2:(n-1)/2, -(n-1)/2:(n-1)/2);
r = hypot( rx , ry );
halffin = fin(round((n-1)/2+1):n);
vq = interp1( mod(n-1,2)/2:mod(n-1,2)/2+length(halffin)-1, halffin , r(:) , 'linear' , 0 );
fout=r; fout(:)=vq;
end