function distances=distance_between_points(x,y,z)
%this function calculates the distance between every possible comination of 
%points described each by (x y z) in 3-D space
    n=numel(x);
    X1=repmat(x,1,n);
    X2=repmat(x',n,1);
    Y1=repmat(y,1,n);
    Y2=repmat(y',n,1);
    Z1=repmat(z,1,n);
    Z2=repmat(z',n,1);

    distances=sqrt((X1-X2).^2+(Y1-Y2).^2+(Z1-Z2).^2);
end