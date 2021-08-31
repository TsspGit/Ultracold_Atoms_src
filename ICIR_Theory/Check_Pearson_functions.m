disp('taylora2f1')
taylora2f1(0.1,0,0.2,0,0.3,0,0.5,0,0.0001)
disp('taylorb2f1')
taylorb2f1(0.1,0,0.2,0,0.3,0,0.5,0,0.0001)
disp('singlefraction2f1')
singlefraction2f1(-1000, 0, -2000, 0, -4000.1, 0, -0.5, 0, 0.00001) 
% not a good method for large a, b
disp('buhringa2f1')
buhringa2f1(5, 0, 7.5, 0, 2.5, 0, 5, 0, 1, 0.00001)
