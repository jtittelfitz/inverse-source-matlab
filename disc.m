function d = disc (c,r,X)
    x = c(1); y = c(2);
    d = (X.x1 - x).^2 + (X.x2 - y).^2 <= r^2;
end