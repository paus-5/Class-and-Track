function v = feed_forward(lambda,B,P,X,n)
B1_vec = diag(B(1:n,:));
x_trunc = X(1:n);
v = max(-1,min(-1./lambda.*B1_vec.*(P.*x_trunc ),0.2));
end