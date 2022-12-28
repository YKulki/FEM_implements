s = tf('s')
G1 = 100*(1+5*s) / (s^4)*(1+s);
nyquist(G1, 'red')