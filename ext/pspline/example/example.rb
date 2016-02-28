require 'pspline'
include PSPLINE
nn = [4, 6]
jj = [3, 2]
ss = [3, 2]
puts "# 2 variable parametric Reisenfeld interpolation"
S = [0, 1, 2, 3, 4]
T = [0, 1, 2, 3, 4, 5]
XYZ = [ [[ 0.0, 5.0, 9.0],[ 0.0, 5.0, 6.0],[ 0.0, 0.5, 6.0],[ 0.0, 0.5, 2.0],[ 0.0, 3.0, 1.0],[ 0.0, 3.0, 0.0]],
        [[ 5.0,-2.5, 6.5],[ 5.0,-2.5, 3.5],[ 0.5,-0.25,5.75],[0.5,-0.25,1.75],[3.0,-1.5,-0.5],[ 3.0,-1.5,-1.5]],
        [[ 0.0,-5.0, 9.0],[ 0.0,-5.0, 6.0],[ 0.0,-0.5, 6.0],[ 0.0,-0.5, 2.0],[ 0.0,-3.0, 1.0],[ 0.0,-3.0, 0.0]],
        [[-5.0, 2.5,11.5],[-5.0, 2.5, 8.5],[-0.5,0.25,6.25],[-0.5,0.25,2.25],[-3.0, 1.5, 2.5],[-3.0, 1.5, 1.5]] ]
Jbn = ARGV[0].to_i
Dp = 8
Rs = Pspline.new(XYZ, [S, T], nn, jj, ss)
for i in 0..3
    for j in 0..5
        printf("<%d,%d>", S[i], T[j])
        u = XYZ[i][j]
        printf(" % .7f % .7f %10.7f\n", u[0], u[1], u[2])
    end
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
    print "\n"
else
    printf ", Jbn = %d\n", Jbn
end
vv = [S, T]
N = Rs.plot(vv, Dp, [Jbn,0]) do |a, b|
    if (a[1] == 0.0)
        printf("%.3f %.3f % .7f % .7f %10.7f\n", a[0], a[1], b[0], b[1], b[2])
    end
end
