%(1)Generate the secret key
x0_l = 0.437;
u0_l = 3.785;
x0_s = [0.364,0.785,0.293];
epsilon = 0.2582;
Iinput=imread('lena512.bmp');
[m,n] = size(Iinput);
% figure;imhist(Iinput);

%(2)Image DNA encoding and substitution 
I = ones(m,4*n);
for i = 1:m
    for j = 1:n
        num2decomposed = Iinput(i,j);
        for z = 1:4
            rem = mod(num2decomposed,4);
            I(i,4*(j-1)+(5-z)) = rem;
            num2decomposed = floor(num2decomposed/4);
        end
    end
end

L = ones(1,4*m*n);
L(1) = x0_l;
for i=1:4*m*n-1
    L(i+1)=logisticChaoticMap(L(i),u0_l);
end
L = reshape(L,m,4*n);

%依据DNA加减法规则，拟定：C=0，A=1，T=2，G=3
L1 = ones(1,4*m*n);
for i=1:4*m*n
    l = mod(floor(L(i)*256),4);
    switch l
        case 0
            L1(i) = 1;
        case 1
            L1(i) = 0;
        case 2
            L1(i) = 3;
        case 3
            L1(i) = 2;
    end
end
L1 = reshape(L1,m,4*n);

DNArules = [1,1,0,0,3,3,2,2;
            0,3,1,2,1,2,0,3;
            3,0,2,1,2,1,3,0;
            2,2,3,3,0,0,1,1];
        
I1 = ones(m,4*n);
for i = 1:m
    e = mode(floor(I(m,4*n)*256),8);
    switch e
        case 0
            for j = 1:4*n
                switch I(i,j)
                    case 0
                        I1(i,j) = 1;
                    case 1
                        I1(i,j) = 0;
                    case 2
                        I1(i,j) = 3;
                    case 3
                        I1(i,j) = 2;
                end
            end
        case 1
            for j = 1:4*n
                switch I(i,j)
                    case 0
                        I1(i,j) = 1;
                    case 1
                        I1(i,j) = 3;
                    case 2
                        I1(i,j) = 0;
                    case 3
                        I1(i,j) = 2;
                end
            end
        case 2
            for j = 1:4*n
                switch I(i,j)
                    case 0
                        I1(i,j) = 0;
                    case 1
                        I1(i,j) = 1;
                    case 2
                        I1(i,j) = 2;
                    case 3
                        I1(i,j) = 3;
                end
            end
        case 3
            for j = 1:4*n
                switch I(i,j)
                    case 0
                        I1(i,j) = 0;
                    case 1
                        I1(i,j) = 2;
                    case 2
                        I1(i,j) = 1;
                    case 3
                        I1(i,j) = 3;
                end
            end
        case 4
            for j = 1:4*n
                switch I(i,j)
                    case 0
                        I1(i,j) = 3;
                    case 1
                        I1(i,j) = 1;
                    case 2
                        I1(i,j) = 2;
                    case 3
                        I1(i,j) = 0;
                end
            end
        case 5
            for j = 1:4*n
                switch I(i,j)
                    case 0
                        I1(i,j) = 3;
                    case 1
                        I1(i,j) = 2;
                    case 2
                        I1(i,j) = 1;
                    case 3
                        I1(i,j) = 0;
                end
            end
        case 6
            for j = 1:4*n
                switch I(i,j)
                    case 0
                        I1(i,j) = 2;
                    case 1
                        I1(i,j) = 0;
                    case 2
                        I1(i,j) = 3;
                    case 3
                        I1(i,j) = 1;
                end
            end
        case 7
            for j = 1:4*n
                switch I(i,j)
                    case 0
                        I1(i,j) = 2;
                    case 1
                        I1(i,j) = 3;
                    case 2
                        I1(i,j) = 0;
                    case 3
                        I1(i,j) = 1;
                end
            end
    end
end

IDNA = ones(m,4*n);
IDNA = mod(I1+L1,4);

%(3)Image entropy computation
A = tabulate(reshape(IDNA,1,4*m*n));
p = reshape(A(:,2),1,4);
H = 0;
for i = 1:4
    H = H + p(i)*log2(1/p(i));
end

if H-floor(H) ~= 0
    x0_H = H-floor(H);
else
    x0_H = x0_l;
end

x_H = ones(1,20);
x_H(1) = x0_H;
for i=1:20
    x_H(i+1)=u0_l*x_H(i)*(1-x_H(i));
end
x20_H = x_H(20);
Hbinary = my_dec2bin(x20_H,64);
Hdecimal = my_bin2dec(Hbinary);

x_l = ones(1,10);
x_l(1) = x0_l;
for i=1:10
    x_l(i+1)=u0_l*x_l(i)*(1-x_l(i));
end
x10_l = x_l(10);

Hcipher = my_bitxor( my_dec2bin(x10_l,64),Hbinary);

%(4)Image permutation
u_s = 3.75 + 0.25*Hdecimal;
x1 = ones(1,4*n);
x2 = ones(1,4*n);
x3 = ones(1,4*n);
x1(1) = (1-epsilon)*logisticChaoticMap(x0_s(1),u0_l)+(epsilon/2)*(logisticChaoticMap(x0_s(2),u0_l)+logisticChaoticMap(x0_s(3),u0_l));
x2(1) = (1-epsilon)*logisticChaoticMap(x0_s(2),u0_l)+(epsilon/2)*(logisticChaoticMap(x0_s(3),u0_l)+logisticChaoticMap(x0_s(1),u0_l));
x3(1) = (1-epsilon)*logisticChaoticMap(x0_s(3),u0_l)+(epsilon/2)*(logisticChaoticMap(x0_s(1),u0_l)+logisticChaoticMap(x0_s(2),u0_l));
for i = 2:4*n
    x1(i) = (1-epsilon)*logisticChaoticMap(x1(i-1),u0_l)+(epsilon/2)*(logisticChaoticMap(x2(i-1),u0_l)+logisticChaoticMap(x3(i-1),u0_l));
    x2(i) = (1-epsilon)*logisticChaoticMap(x2(i-1),u0_l)+(epsilon/2)*(logisticChaoticMap(x3(i-1),u0_l)+logisticChaoticMap(x1(i-1),u0_l));
    x3(i) = (1-epsilon)*logisticChaoticMap(x3(i-1),u0_l)+(epsilon/2)*(logisticChaoticMap(x1(i-1),u0_l)+logisticChaoticMap(x2(i-1),u0_l));
end
R = x1(:,1:m);
C = x3;
[R1,IndR] = sort(R);
[C1,IndC] = sort(C);

IDNA1 = ones(m,4*n);
for i = 1:m
    for j = 1:4*n
        IDNA1(i,j) = IDNA(IndR(i),IndC(j));
    end
end

%(5)Image decoding
IDNA1dec = ones(m,4*n);
for i = 1:m
    for j = 1:4*n
        switch IDNA1(i,j)
            case 0
                IDNA1dec(i,j) = 1;
            case 1
                IDNA1dec(i,j) = 0;
            case 2
                IDNA1dec(i,j) = 3;
            case 3
                IDNA1dec(i,j) = 2;
        end    
    end
end

Icipher = ones(m,n);
sign = 1;
num = 0;
for i = 1:m
    for j = 1:4*n
        switch mod(j,4)
            case 1
                num = num + IDNA1(i,j)*64;
            case 2
                num = num + IDNA1(i,j)*16;
            case 3
                num = num + IDNA1(i,j)*4;
            case 0
                num = num + IDNA1(i,j)*1;
                Icipher(i,sign) = num;
                num = 0;
                sign = sign + 1;
        end
    end
    sign = 1;
end

Icipher = mat2gray(Icipher);
imshow(Icipher);
figure;imhist(Iinput);
hold on
figure;imhist(Icipher);