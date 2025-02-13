clc
clear


% number of users and targets
K=4;
M=2;


%transmit power and noise
Pt=db2pow(10-30); %dBm
noise_c=db2pow(0-30); %dBm
noise_s=db2pow(0-30); %dBm
L=30;
kappa=2*L/noise_s;

tolerance=1e-5;

I_in=5;
I_out=100;

SR_all=zeros(I_out,I_in);
CRB_all=zeros(I_out,I_in);
Time_all=zeros(I_out,I_in);
Obj_all=zeros(I_out,I_in);


for k_par=150
    %:I_out*I_in

    channel=mod((k_par-1),I_out)+1;
    weight=floor((k_par-1)/I_out)+1;

    disp([channel,weight])
    delta_s=1;
    delta_c=0.25;

    %antenna for transmit
    Nth=4;
    Ntv=2^(weight+1);
    Nt=Nth*Ntv;

    %antenna for receiver
    Nrh=5;
    Nrv=2^(weight+1);
    Nr=Nrh*Nrv;


    num_sensing_streams=3*M;
    rng(1,'twister')
    alpha=0.1*(1+0.2*randn(M,1)).*exp(1j*2*pi*rand(M,1));
    theta=-pi/3+2*pi/3*rand(M,1);
    phi=-pi/3+2*pi/3*rand(M,1);

    H_all=1/sqrt(2)*(randn(I_out,Nt,K)+1j*randn(I_out,Nt,K));

    H=squeeze(H_all(channel,:,:));


    tic
    [A,dAtheta,dAphi]=construct_steer_matrix_and_derivative_steer_matrix(theta,phi,Nth,Ntv);
    [B,dBtheta,dBphi]=construct_steer_matrix_and_derivative_steer_matrix(theta,phi,Nrh,Nrv);

    U=diag(alpha);
    RS=[H,A,dAtheta,dAphi];
    RS_p=pinv(RS);
    RS_co=RS'*RS;
    H_e=RS'*H;

    A_e=RS'*A;
    dAtheta_e=RS'*dAtheta;
    dAphi_e=RS'*dAphi;


    Wc=delta_c*H./vecnorm(H);

    %Wc=randn(Nt,K)+1j*randn(Nt,K);
    % Ws=delta_s*repmat(eigvector,1,num_sensing_streams);
    Ws=randn(Nt,num_sensing_streams)+1j*randn(Nt,num_sensing_streams);

    % Ws=initial_Ws(L,noise_s,Nt,num_sensing_streams,A,dAtheta,dAphi,B,dBtheta,dBphi,U);
    W=[Wc,delta_s*Ws];


    W=W*sqrt(Pt/trace(W*W'));


    P=RS_p*W;

    FIM=calculateFIM(L,noise_s,P*P',A_e,dAtheta_e,dAphi_e,B,dBtheta,dBphi,U);
    P_last=P;
    Con=[];

    for count=1:4000

        T_k=sum(square_abs(H_e'*P(:,1:K)),2)+noise_c*ones(K,1);
        alpha_k=T_k./(T_k-square_abs(diag(H_e'*P(:,1:K))))-1;
        beta_k=sqrt(1+alpha_k).*diag(H_e'*P(:,1:K))./T_k;
        Sigma1=diag(sqrt(1+alpha_k).*beta_k);
        Sigma2=diag(square_abs(beta_k));

        CRBM=inv(FIM);

        % PHI=CRBM*CRBM;
        Q=construct_matrixQ(L,noise_s,CRBM*CRBM,A_e,dAtheta_e,dAphi_e,B,dBtheta,dBphi,U);
        % rank(Q+Q')
        C1=[delta_c*H_e*Sigma1,zeros(size(RS,2),num_sensing_streams)];
        C2=0.5*delta_s*(Q+Q')-delta_c*H_e*Sigma2*H_e';

        %dom eigvalue

        mu=abs(eigs(H_e*Sigma2*H_e',1,'LM'));
        C2=delta_c*mu*RS_co+C2;

       for iter=1:20
        Linear=C1+C2*P;
        P=RS_co\Linear;
        P=projection_ellipsoid(RS, RS_p, P, Pt);
       end
       

        FIM=calculateFIM(L,noise_s,P*P',A_e,dAtheta_e,dAphi_e,B,dBtheta,dBphi,U);

        T_k=sum(square_abs(H_e'*P),2)+noise_c*ones(K,1);

        rate_p=log(T_k./(T_k-square_abs(diag(H_e'*P))));
        % CRBM=delta_s*trace(inv(FIM));
        obj=delta_c*sum(log(T_k./(T_k-square_abs(diag(H_e'*P)))))-delta_s*trace(inv(FIM));
        Con=[Con;[sum(log(T_k./(T_k-square_abs(diag(H_e'*P))))),-trace(inv(FIM)),obj]];

        if norm(RS*P-RS*P_last)<tolerance
            break
        else
            P_last=P;
        end

    end

    W=RS*P;
    FIM=calculateFIM(L,noise_s,W*W',A,dAtheta,dAphi,B,dBtheta,dBphi,U);
    CRB=inv(FIM);
    SR=sum(log(T_k./(T_k-square_abs(diag(H'*W(:,1:K))))));

    CRB_all(k_par)=trace(inv(FIM));
    SR_all(k_par)=SR;
    Time_all(k_par)=toc;
    Obj_all(k_par)=delta_c*SR-delta_s*trace(inv(FIM));

end

% save data_SCA_Ns=Nt

toc
figure(1)
plot(Con(:,1))
xlabel('iteration')
ylabel('Sum Rate (nat/s/Hz)')
figure(2)
plot(Con(:,2))
xlabel('iteration')
ylabel('trace of inverse of FIM')
figure(3)
plot(Con(:,3))
xlabel('iteration')
ylabel('total objective value')


function Ws=initial_Ws(L,noise_s,Nt,num_sensing_streams,A,dAtheta,dAphi,B,dBtheta,dBphi,U)

Ws=randn(Nt,num_sensing_streams)+1j*randn(Nt,num_sensing_streams);


for iter=1:5
    FIM=calculateFIM(L,noise_s,Ws,A,dAtheta,dAphi,B,dBtheta,dBphi,U);
    CRBM=inv(FIM);
    Q=construct_matrixQ(L,noise_s,CRBM*CRBM,A,dAtheta,dAphi,B,dBtheta,dBphi,U);
    Ws=(Q+Q')*Ws;
    Ws=Ws*sqrt(1/trace(Ws*Ws'));
end


end

function X = projection_ellipsoid(A, A_p, X0, Pt)
    
W=A*X0;
W=W*sqrt(Pt/trace(W*W'));
X=A_p*W;
end

function Q=construct_matrixQ(L,noise_s,Phi,A,dAtheta,dAphi,B,dBtheta,dBphi,U)
M=size(U,1);

phi11=Phi(1:M,1:M);
phi12=Phi(1:M,M+1:2*M);
phi13=Phi(1:M,2*M+1:3*M);
phi14=Phi(1:M,3*M+1:4*M);
phi22=Phi(M+1:2*M,M+1:2*M);
phi23=Phi(M+1:2*M,2*M+1:3*M);
phi24=Phi(M+1:2*M,3*M+1:4*M);
phi33=Phi(2*M+1:3*M,2*M+1:3*M);
phi34=Phi(2*M+1:3*M,3*M+1:4*M);
phi44=Phi(3*M+1:4*M,3*M+1:4*M);
% Compute Q11

Q11 = (A * U' * (phi11 .* (dBtheta' * dBtheta)) * U * A')  + ...
    (dAtheta * U' * (phi11.* (B' * dBtheta)) * U * A')  + ...
    (A * U' * (phi11.* (dBtheta' * B)) * U * dAtheta' )  + ...
    (dAtheta * U' * (phi11.* (B' * B)) * U * dAtheta') ;

% Compute F12
Q12 = (A * U' * (2*phi12.* (dBtheta' * dBphi)) * U * A' )  + ...
    ( dAtheta * U' * (2*phi12.* (B' * dBphi)) * U * A')  + ...
    (A * U' * (2*phi12.* (dBtheta' * B)) * U * dAphi') + ...
    (dAtheta * U' * (2*phi12.* (B' * B)) * U * dAphi') ;

% Compute F13
Q13 = (A * U' * ((2*phi13+2j*phi14).* (dBtheta' * B)) * A' )  + ...
    (dAtheta * U' * ((2*phi13+2j*phi14).* (B' * B)) * A') ;

% Compute F22
Q22 = (A * U' * (phi22.* (dBphi' * dBphi)) * U * A' )  + ...
    (dAphi * U' * (phi22.* (B' * dBphi)) * U * A')  + ...
    (A * U' * (phi22.* (dBphi' * B)) * U * dAphi' )  + ...
    (dAphi * U' * (phi22.* (B' * B)) * U * dAphi' ) ;

% Compute F23
Q23 = (A * U' * ((2*phi23+2j*phi24).* (dBphi' * B)) *A')  + ...
    (dAphi * U' * ((2*phi23+2j*phi24).* (B' * B)) *A' ) ;

% Compute F33
Q33 = ( A * ((phi33+phi44+2j*phi34).* (B' * B))*A') ;

Q=2*L/noise_s*(Q11+Q12+Q13+Q22+Q23+Q33);
end


function FIM=calculateFIM(L,noise_s,Rx,A,dAtheta,dAphi,B,dBtheta,dBphi,U)

% Inputs
% A: Matrix A
% U: Matrix U
% dAtheta: Partial derivative of A with respect to theta
% dAphi: Partial derivative of A with respect to phi
% dBtheta: Partial derivative of B with respect to theta
% dBphi: Partial derivative of B with respect to phi
% B: Matrix B
% Rx: Matrix Rx


% Compute F11
F11 = (U * A' * Rx * A * U').' .* (dBtheta' * dBtheta) + ...
    (U * A' * Rx * dAtheta * U').' .* (B' * dBtheta) + ...
    (U * dAtheta' * Rx * A * U').' .* (dBtheta' * B) + ...
    (U * dAtheta' * Rx * dAtheta * U').' .* (B' * B);

% Compute F12
F12 = (U * A' * Rx * A * U').' .* (dBtheta' * dBphi) + ...
    (U * A' * Rx * dAtheta * U').' .* (B' * dBphi) + ...
    (U * dAphi' * Rx * A * U').' .* (dBtheta' * B) + ...
    (U * dAphi' * Rx * dAtheta * U').' .* (B' * B);

% Compute F13
F13 = (A' * Rx * A * U').' .* (dBtheta' * B) + ...
    (A' * Rx * dAtheta * U').' .* (B' * B);

% Compute F22
F22 = (U * A' * Rx * A * U').' .* (dBphi' * dBphi) + ...
    (U * A' * Rx * dAphi * U').' .* (B' * dBphi) + ...
    (U * dAphi' * Rx * A * U').' .* (dBphi' * B) + ...
    (U * dAphi' * Rx * dAphi * U').' .* (B' * B);

% Compute F23
F23 = (A' * Rx * A * U').' .* (dBphi' * B) + ...
    (A' * Rx * dAphi * U').' .* (B' * B);

% Compute F33
F33 = (A' * Rx * A).' .* (B' * B);


M=size(U,1);
FIM=zeros(4*M);


FIM(1:M,1:M)=real(F11);
FIM(1:M,M+1:2*M)=real(F12);
FIM(1:M,2*M+1:3*M)=real(F13);
FIM(1:M,3*M+1:4*M)=-imag(F13);
FIM(M+1:2*M,M+1:2*M)=real(F22);
FIM(M+1:2*M,2*M+1:3*M)=real(F23);
FIM(M+1:2*M,3*M+1:4*M)=-imag(F23);
FIM(2*M+1:3*M,2*M+1:3*M)=real(F33);
FIM(2*M+1:3*M,3*M+1:4*M)=-imag(F33);
FIM(3*M+1:4*M,3*M+1:4*M)=real(F33);

FIM=triu(FIM)+triu(FIM)'-FIM.*eye(4*M);

FIM=2*L/noise_s*FIM;
end




function [A,dAtheta,dAphi]=construct_steer_matrix_and_derivative_steer_matrix(theta,phi,Mx,My)
% input:
% theta - Elevation angle (in radians)
% phi   - Azimuth angle (in radians)
% Mx    - Number of elements along x-axis
% My    - Number of elements along y-axis

% Generate indices for x-axis and y-axis
ix = (0:Mx-1)'; % Column vector for x-axis indices
iy = (0:My-1)'; % Column vector for y-axis indices

% Number of steering vectors
M=length(theta);
% Preallocate outputs
A=zeros(Mx*My,M);
dAtheta=zeros(Mx*My,M);
dAphi=zeros(Mx*My,M);
% Compute steering vectors along x and y axes
for m=1:M
    % Steering vectors for x and y dimensions
    ax =1/sqrt(Mx)* exp(1j * pi * ix * sin(theta(m)) * sin(phi(m)));
    ay =1/sqrt(My)* exp(1j * pi * iy * cos(phi(m)));
    % Derivatives of ax and ay
    daxtheta=1j*pi*ix*cos(theta(m))*sin(phi(m)).*ax;

    daxphi=1j*pi*ix*sin(theta(m))*cos(phi(m)).*ax;
    dayphi=-1j*pi*iy*sin(phi(m)).*ay;
    % Compute the full UPA steering vector using Kronecker product
    A(:,m) = kron(ay, ax); % Kronecker product combines ax and ay
    dAtheta(:,m)=kron(ay,daxtheta);
    dAphi(:,m)=kron(ay,daxphi)+kron(dayphi,ax);
end

end
