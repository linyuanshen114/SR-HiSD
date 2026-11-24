function Y = onsager_project(Y,U,w)
Y = Y-U*(U'*Y*w);
end

