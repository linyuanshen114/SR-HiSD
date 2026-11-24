function R = onsager_rotation_matrix(a,t)
a = a(:)/norm(a);
S = [0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
R = eye(3)+sin(t)*S+(1-cos(t))*S*S;
end

