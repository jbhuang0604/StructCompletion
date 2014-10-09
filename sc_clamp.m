function v = sc_clamp(v, vLowerB, vUpperB)

% SC_CLAMP: clamp value v

v(v<vLowerB) = vLowerB;
v(v>vUpperB) = vUpperB;

end