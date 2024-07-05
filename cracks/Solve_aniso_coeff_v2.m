clc
close all

%===============================================================================
% Read elastic constants
datafile = "./bench/data/results.txt";
data = readlines(datafile) ;  
c11 = split(data(18));
c12 = split(data(19));
c44 = split(data(20));


c_11 = double(c11(5));
c_12 = double(c12(5));
c_44 = double(c44(5));
surf100 = double(data(28));
surf110 = double(data(31)); 

folder_name = "lefm_coeffs";
disp("Crack sysem includes 1: (100)[010]; 2: (100)[011]; 3: (110)[001]; 4: (110)[1-10]")
% Define elastic constants
C=[c_11,c_12,c_12,0,0,0;
   c_12,c_11,c_12,0,0,0;
   c_12,c_12,c_11,0,0,0;
   0,0,0,c_44,0,0;
   0,0,0,0,c_44,0;
   0,0,0,0,0,c_44;];

status = mkdir(folder_name);
% Write the parameter solution 
% Loop for four crack systems and output the results.
for cs=1:4
%------------------------------------
% Choose crystal orientation of the crack system 
% system 1: (100)[010]
% system 2: (100)[011]
% system 3: (110)[001]
% system 4: (110)[1-10]
%------------------------------------
    if cs == 1
        a1=[0 0 1];
        a2=[1 0 0];
        a3=[0 1 0];
        surfE = surf100;
    elseif cs == 2
        a1=[0 -1 1];
        a2=[1 0 0];
        a3=[0 1 1];
        surfE = surf100;
    elseif cs == 3
        a1=[1 -1 0];
        a2=[1 1 0];
        a3=[0 0 1];
        surfE = surf110;
    elseif cs == 4
        a1=[0 0 -1];
        a2=[1 1 0];
        a3=[1 -1 0];
        surfE = surf110;
    end
    % calculation
    [ s,p,q,K_I,G_I ] = aniso_disp_solution( C,a1,a2,a3,surfE );

    fileID = fopen(folder_name+"/Aniso_paras."+int2str(cs),'w');
    fprintf(fileID,'%30s \r\n','Paramters list real and imag for s p q');
    fprintf(fileID,'%5s %12.8f %12.8f\r\n','s1=',real(s(1)),imag(s(1)));
    fprintf(fileID,'%5s %12.8f %12.8f\r\n','s2=',real(s(2)),imag(s(2)));
    fprintf(fileID,'%5s %12.8f %12.8f\r\n','p1=',real(p(1)),imag(p(1)));
    fprintf(fileID,'%5s %12.8f %12.8f\r\n','p2=',real(p(2)),imag(p(2)));
    fprintf(fileID,'%5s %12.8f %12.8f\r\n','q1=',real(q(1)),imag(q(1)));
    fprintf(fileID,'%5s %12.8f %12.8f\r\n','q2=',real(q(2)),imag(q(2)));
    fprintf(fileID,'%5s %12.8f %10s\r\n','G_I=',G_I,'J*m^2');
    fprintf(fileID,'%5s %12.8f %10s\r\n','K_I=',K_I,'MPa*m^1/2');
    fclose(fileID);
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
    % Write the parameter for lammps implementation
    fileID = fopen(folder_name+"/lefm_paras.CrackSystem_"+int2str(cs),'w');
    fprintf(fileID,'%30s \r\n','#Paramters list real and imag for s p q');
    fprintf(fileID,'%30s \r\n','#--------------------------------------');
    fprintf(fileID,'%30s %12.8f \r\n','variable        s1_real equal',real(s(1)));
    fprintf(fileID,'%30s %12.8f \r\n','variable        s1_imag equal',imag(s(1)));
    fprintf(fileID,'%30s %12.8f \r\n','variable        s2_real equal',real(s(2)));
    fprintf(fileID,'%30s %12.8f \r\n','variable        s2_imag equal',imag(s(2)));
    fprintf(fileID,'%30s %12.8f \r\n','variable        p1_real equal',real(p(1)));
    fprintf(fileID,'%30s %12.8f \r\n','variable        p1_imag equal',imag(p(1)));
    fprintf(fileID,'%30s %12.8f \r\n','variable        p2_real equal',real(p(2)));
    fprintf(fileID,'%30s %12.8f \r\n','variable        p2_imag equal',imag(p(2)));
    fprintf(fileID,'%30s %12.8f \r\n','variable        q1_real equal',real(q(1)));
    fprintf(fileID,'%30s %12.8f \r\n','variable        q1_imag equal',imag(q(1)));
    fprintf(fileID,'%30s %12.8f \r\n','variable        q2_real equal',real(q(2)));
    fprintf(fileID,'%30s %12.8f \r\n','variable        q2_imag equal',imag(q(2)));
    fprintf(fileID,'%30s \r\n','#--------------------------------------');
    fprintf(fileID,'%5s %12.8f %10s\r\n','#G_I=',G_I,'J*m^2');
    fprintf(fileID,'%5s %12.8f %10s\r\n','#K_I=',K_I,'MPa*m^1/2');
    fprintf(fileID,'%30s \r\n','#--------------------------------------');
    fclose(fileID);
end

exit;
%-----------------------------------------------------------------------
%===============================================================================
% Define the function used to solve the LEFM crack tip coefficients.
function [ s,p,q,K_I,G_I ] = aniso_disp_solution( C,a1,a2,a3,surfE )
    % Define the compliance matrix
    S=inv(C);
    % transform the compliance matrix from default corrdinates to current ones by rotation.
    % default corrdinates alinged with crystal orientation: [100][010][001]
    % current corrdinates alinged with crystal orientation: [abc][def][ghi]
    % define rotation matrix Q
    Q=[a1/sqrt(dot(a1,a1));a2/sqrt(dot(a2,a2));a3/sqrt(dot(a3,a3))];
    % define the transformation matrix in Voigt notation
    K1=[Q(1,1).^2,Q(1,2).^2,Q(1,3).^2;
        Q(2,1).^2,Q(2,2).^2,Q(2,3).^2;
        Q(3,1).^2,Q(3,2).^2,Q(3,3).^2];
    K2=[Q(1,2)*Q(1,3),Q(1,3)*Q(1,1),Q(1,1)*Q(1,2);
        Q(2,2)*Q(2,3),Q(2,3)*Q(2,1),Q(2,1)*Q(2,2);
        Q(3,2)*Q(3,3),Q(3,3)*Q(3,1),Q(3,1)*Q(3,2)];
    K3=[Q(2,1)*Q(3,1),Q(2,2)*Q(3,2),Q(2,3)*Q(3,3);
        Q(3,1)*Q(1,1),Q(3,2)*Q(1,2),Q(3,3)*Q(1,3);
        Q(1,1)*Q(2,1),Q(1,2)*Q(2,2),Q(1,3)*Q(2,3)];
    K4=[Q(2,2)*Q(3,3)+Q(2,3)*Q(3,2),Q(2,3)*Q(3,1)+Q(2,1)*Q(3,3),Q(2,1)*Q(3,2)+Q(2,2)*Q(3,1);
        Q(3,2)*Q(1,3)+Q(3,3)*Q(1,2),Q(3,3)*Q(1,1)+Q(3,1)*Q(1,3),Q(3,1)*Q(1,2)+Q(3,2)*Q(1,1);
        Q(1,2)*Q(2,3)+Q(1,3)*Q(2,2),Q(1,3)*Q(2,1)+Q(1,1)*Q(2,3),Q(1,1)*Q(2,2)+Q(1,2)*Q(2,1)];
    K=[K1,2*K2;
        K3,K4];
    % rotate the compliance matrix
    S_star=inv(K).'*S*inv(K);
    %--------------------------------------------------------------------
    % General elastic solution
    % Define the coeffience of fourth order PDE for plane strain problem
    b_11=(S_star(1,1)*S_star(3,3)-S_star(1,3).^2)/S_star(3,3);
    b_22=(S_star(2,2)*S_star(3,3)-S_star(2,3).^2)/S_star(3,3);
    b_66=(S_star(6,6)*S_star(3,3)-S_star(3,6).^2)/S_star(3,3);
    b_12=(S_star(1,2)*S_star(3,3)-S_star(1,3)*S_star(2,3))/S_star(3,3);
    b_21=b_12;
    b_16=(S_star(1,6)*S_star(3,3)-S_star(1,3)*S_star(3,6))/S_star(3,3);
    b_61=b_16;
    b_26=(S_star(2,6).*S_star(3,3)-S_star(2,3).*S_star(3,6))/S_star(3,3);
    b_62=b_26;
    % define matrix b for computing G and K.
    b=zeros(6,6);
    b(1,1)=b_11;
    b(2,2)=b_22;
    b(6,6)=b_66;
    b(1,2)=b_12;
    b(2,1)=b_12;
    b(1,6)=b_16;
    b(6,1)=b_16;
    b(2,6)=b_26;
    b(6,2)=b_26;
    B=sqrt((b_11*b_22/2)*(sqrt(b_22/b_11)+((2*b_12+b_66)/(2*b_11))));
    disp(B);
    K_I=sqrt(2*surfE*(1/(B*1000))); % 1000 is due to the unit Gpa(C)->Mpa(K).
    G_I=2*surfE;
    %--------------------------------------------------------------------
    % Solve the characteristic equation ('b_11*x^4-2*b_16*x^3+(2*b_12+b_66)*x^2-2*b_26*x+b_22=0')
    coefvct = [b_11 -2.*b_16 2.*b_12+b_66 -2.*b_26 b_22];     % Coefficient Vector
    rt = roots(coefvct);
    s = rt;
    s(imag(s)<0) = [];
    if real(s(1))<real(s(2))
        exchange_s=s(2);
        s(2)=s(1);
        s(1)=exchange_s;
    end
    % Define some constants that will be used for the stress functions
    p(1)=b_11.*s(1).^2+b_12-b_16.*s(1);
    p(2)=b_11.*s(2).^2+b_12-b_16.*s(2);
    q(1)=b_12.*s(1)+b_22/s(1)-b_26;
    q(2)=b_12.*s(2)+b_22/s(2)-b_26;
end
%===============================================================================
