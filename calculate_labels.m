function F=calculate_labels(W,Y,alpha)
    F=(1-alpha)*pinv(eye(size(W,1))-alpha*W)*Y;
end