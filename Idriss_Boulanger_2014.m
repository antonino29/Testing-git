function [Ic,Fs,PL,epsV,cedV,gamma_max,Delta_LDI,Delta_LD,spessore,Delta_LPI,Delta_LSI,Delta_LSN,wfspessore,CSR75,CRR75,zass]=Idriss_Boulanger_2014(amax,M,lat,long,z0,a,gamma,gamma_sat,water,conf_geom,configurazione,p,z,qc,f_s,u,n,truncName,folderCPTout,folderCPTout1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo del Fattore di Sicurezza (Fs), dell'Indice del Potenziale di
% Liquefazione LPI, dell'Indice di Severità di Liquefazione LSI, dei
% Cedimenti Verticali del terreno e delle Espansioni Laterali del terreno
% a partire da dati di prove CPT secondo le indicazioni suggerite in Seed(2010).
% Nello specifico:
%   - Fs secondo Idriss & Boulanger (2008, 2014);
%   - PL secondo Boulanger & Idriss (2014, 2015);
%   - CedV secondo Yoshimine (2006);
%   - LD secondo Idriss & Boulanger (2008).
%
%N.B.: metodo specifico per prova penetrometrica statica con punta elettrica.
%
% (Bozzoni Francesca, 19/11/2013)
% (Famà Antonino, 15/3/2017) controllo e aggiunta accorgimenti sulla base del
% confronto con CLiq
% (Famà Antonino 08/06/2018 aggiunta correzione Ic per prove CPTu messa a punto da Meisina et al 2018
% e modifiche per figura Fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Inizializzazione delle variabili
    sigV_eff = zeros(n,1);
    sigV_t = zeros(n,1);
    u0 = zeros(n,1);
    
    qt = zeros(n,1);
    Rf = zeros(n,1);
    
    errore_n_stress=zeros(n,1);
    n_stress_calcolo=zeros(n,1);
    Cn_calcolo= zeros(n,1);
    F_r=zeros(n,1);
    Q=zeros(n,1);
    Ic_calcolo=zeros(n,1);
    Qtn = zeros(n,1);
    FC = zeros(n,1);
    Ic = zeros(n,1);
    n_stress = zeros(n,1);

    qc1n=zeros(n,1);
    qc1ncs=zeros(n,1);
    delta_qc1n=zeros(n,1);
    alfa_rd = zeros(n,1);
    beta_rd=zeros(n,1);
    rd = zeros(n,1);
    CSR = zeros(n,1);
    MSF = zeros (n,1);
    CSR_eq = zeros(n,1);
    CSR75 = zeros(n,1);
    CRR75 = zeros(n,1);
    Csig = zeros(n,1);
    Ksig = zeros(n,1);
    CRR = zeros(n,1);
    Fs = zeros(n,1);
    deltaU = zeros(n,1);
    PL = zeros(n,1);
    F_alfa = zeros(n,1);
    gamma_lim = zeros(n,1);
    epsV = zeros(n,1);
    cedV = zeros(n,1);
    gamma_max = zeros(n,1);
    Delta_LDI = zeros(n,1);
    Delta_LD = zeros(n,1);
    Delta_LPI = zeros(n,1);
    Delta_LSI = zeros(n,1);
    Delta_LSN = zeros(n,1);
    spessore = zeros(n,1);
        
    zass=z0-z;% AF 08/06/2018 aggiunta quota assoluta
    spessore(1) = z(2)-z(1);
    
    
    for ii=2:n
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %qt resistenza di punta corretta con rapporto aree a
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        qt(ii)=qc(ii)+u(ii)*(1-a);
        Rf(ii)=(f_s(ii)./qt(ii))*100;
        z_calcolo(ii)=z(ii);
        
        %Calcolo del peso per unità di volume da prove CPT da Robertson e
        %Cabal 2010 %Antonino Famà 27/02/2018
        if gamma(ii)==0
            gamma(ii)=9.807*(0.27*log10(Rf(ii))+0.36*(log10(qt(ii)/0.101325))+1.236);
        end
        if gamma_sat(ii)==0
            gamma_sat(ii)=gamma(ii);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calcolo Sforzo Verticale Efficace sig1v e Sforzo Verticale Totale sigv
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        spessore(ii) = z_calcolo(ii)-z_calcolo(ii-1);
        if z(ii)<water
            sigV_t(ii) = z(ii)*gamma(ii);
            u0(ii)= 0;
        else
            sigV_t(ii) = water*gamma(ii)+(z(ii)-water)*gamma_sat(ii);
            u0(ii)= (z(ii)-water)*9.806;
        end
        sigV_eff(ii) = sigV_t(ii)-max((z(ii)-water),0)*9.81;
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Valutazione tipo di terreno - Soil Behaviour Type norm. (Robertson, 2009)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        n_stress_calcolo(ii)=0.5;
        errore_n_stress(ii)=0.05;
        while abs(errore_n_stress(ii))>0.00001
            n_stress_0(ii)=n_stress_calcolo(ii);
            Cn_calcolo(ii)=min(1.7,(100/sigV_eff(ii))^n_stress_calcolo(ii));
%             Cn_calcolo(ii)=(100/sigV_eff(ii))^n_stress_calcolo(ii);
            Q(ii)=((qt(ii)*1000-sigV_t(ii))/100)*Cn_calcolo(ii);
            F_r(ii)=(1000*f_s(ii)/(1000*qt(ii)-sigV_t(ii)))*100;
            Ic_calcolo(ii)=((3.47-log10(Q(ii)))^2+(1.22+log10(F_r(ii)))^2)^0.5; %Ic indice di comportamento del terreno
            n_stress_calcolo(ii)=min(1,0.381*Ic_calcolo(ii)+0.05*(sigV_eff(ii)/100)-0.15);
%             n_stress_calcolo(ii)=0.5;
            errore_n_stress(ii)=n_stress_0(ii)-n_stress_calcolo(ii);
        end
        n_stress(ii)=n_stress_calcolo(ii);
        Ic(ii)=Ic_calcolo(ii);
        Qtn(ii)=Q(ii);
%AF 16/6/2017 Aggiunta limitazione ai valori di Ic, per valori di Rf e Q
%negativi erano stati riscontrati valori complessi di Ic
        if Q(ii)<0 || F_r(ii)<0
                Ic(ii)=5;
        elseif f_s(ii)==0  %aggiunta limitazione sul valore di Ic nel caso di fs=0, questo generava valori NaN in Ic (65535 nei files excel)
           Ic(ii)=4;
        end

        %stima FC con Boulanger & Idriss (2014)
        %verifiche con Cfc variabile, in assenza di prove di lab. o campioni di terreno, per
        %verificare l'effettiva sensibilità di FC alla variazione
        Cfc=0;
        FC(ii) = (80*(Ic(ii)+Cfc))-137;  
        if FC(ii)>=100    %FC deve essere compreso tra lo 0% e il 100%
           FC(ii)=100;
        end
        if FC(ii)<=0
           FC(ii)=0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calcolo fattore di sicurezza =CRR/CSR
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
            
            %Calcolo CSR(Cycling Stress Ratio)
            
            %rd con Idriss (1999) come suggerito in Idriss & Boulanger (2008, 2014)
            alfa_rd(ii) = -1.012-1.126*sin((z(ii)/11.73)+5.133);
            beta_rd(ii) = 0.106+0.118*sin((z(ii)/11.28)+5.142);
            rd(ii) = exp(alfa_rd(ii)+(beta_rd(ii)*M));
         
            CSR(ii) = 0.65*amax*rd(ii)*(sigV_t(ii)/sigV_eff(ii));
            
            %Calcolo CRR (Cycling Resistance Ratio) da risultati prove CPT 
            %procedura di normalizzazione iterativa come in Idriss & Boulanger (2008, 2014)
            errore(ii)=0.05;
            qc1n_primo(ii)=(qt(ii)*1000/101.325);
            m_primo(ii)=1.338-0.249*max(21,min(254,(qc1n_primo(ii))))^0.264;
            m_calcolo(ii)= m_primo(ii);
            while abs(errore(ii))>0.00001
                m_0(ii)=m_calcolo(ii);
                Cn_calcolo(ii)=min(1.7,((101.325/sigV_eff(ii))^m_calcolo(ii)));
                qc1n_calcolo(ii)=(qt(ii)*1000/101.325*Cn_calcolo(ii));
                
                %correzione in base al contenuto di fine FC; Idriss & Boulanger (2014):
                delta_qc1n_calcolo(ii) = (11.9+qc1n_calcolo(ii)/14.6)*exp(1.63-(9.7/(FC(ii)+2))-(((15.7/(FC(ii)+2))^2)));
                qc1ncs_calcolo(ii)=qc1n_calcolo(ii)+delta_qc1n_calcolo(ii);
                
                m_calcolo(ii)=1.338-0.249*(max(21,(min(254,(qc1ncs_calcolo(ii))))))^0.264;
                errore(ii)=m_0(ii)-m_calcolo(ii);
            end
            m(ii)=m_calcolo(ii);
            Cn(ii)=min(1.7,Cn_calcolo(ii));
            qc1n(ii)=qt(ii)*1000/101.325*Cn(ii);
            delta_qc1n(ii)=delta_qc1n_calcolo(ii);
            if Ic(ii) > 2.6
            delta_qc1n(ii)=0;
            else
            delta_qc1n(ii)=delta_qc1n_calcolo(ii);
            end
            qc1ncs(ii)=qc1n(ii)+delta_qc1n(ii);
      
            %MSF con Boulanger & Idriss (2014)
            MSFmax(ii) = min(2.2,(1.09+((qc1ncs(ii)/180)^3)));
            MSF(ii) = 1+((MSFmax(ii)-1)*(8.64*exp((-M)/4)-1.325));
            
            
            %Idriss & Boulanger (2008) consigliano Idriss & Boulanger
            %(2004); CRR75 aggiornato a Boulanger & Idriss (2014):
            
            %AF 24/7/2018 ATTENZIONE IL METODO NON LIMITA il valore di
            %CRR75 in funzione di qc1ncs (211)
            if qc1ncs(ii)>=211
                CRR75(ii)= 2;
            elseif qc1ncs(ii)<211
                CRR75(ii) = exp((qc1ncs(ii)/113)+((qc1ncs(ii)/1000)^2)-((qc1ncs(ii)/140)^3)+((qc1ncs(ii)/137)^4)-2.80);
            end
            
            %Calcolo Fs
            %Idriss & Boulanger (2008) restrizione di Ksig a 1.1 come in Boulanger & Idriss (2014)
            Csig(ii)=min(0.3,1/(37.3-8.27*((min(qc1n(ii),211))^0.264)));
            Ksig(ii)=min(1.1,(1-Csig(ii)*log(sigV_eff(ii)/101.325)));
            CSR_eq(ii)=CSR(ii)/MSF(ii);
            CSR75(ii) = CSR(ii)/(MSF(ii)*Ksig(ii)); 
            CRR(ii) = min(2,(CRR75(ii)*MSF(ii)*Ksig(ii)));
%             if CRR(ii)/CSR(ii)>2
%                 Fs(ii) = 2;  %AF 07/06/18 modifica da Fs = 2 Fs = 0 solamente per un miglioramento nella figura 
%             else
                Fs(ii) = CRR75(ii)/CSR75(ii);
%             end
          
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Limite casi non liquefacibili
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if  Ic(ii) > 2.6 || z(ii)<= water
%                 CSR75(ii)=2;
%                 CRR75(ii)=4;
                Fs(ii)=100; %AF 07/06/18 modifica da Fs = 2 Fs = 0 solamente per un miglioramento nella figura
                PL(ii)=0;
            end
            
            %Calcolo DeltaU Ishihara 1985
            if Fs(ii)>=1 
            deltaU(ii)=0.8*Fs(ii)^-5.574;
            else
            deltaU(ii)=1;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calcolo PL  come in Boulanger & Idriss (2014, 2015)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if Fs(ii)==0 || Fs(ii)==2
                PL(ii)=0;
            else
            X(ii)=-(((qc1ncs(ii)/113)+((qc1ncs(ii)/1000)^2)-((qc1ncs(ii)/140)^3)+((qc1ncs(ii)/137)^4)-2.60-log(CSR75(ii)))/0.20);
            PL(ii) = normcdf(X(ii),0,1); %probabilità di liquefazione
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calcolo dei Cedimenti Verticali del terreno - Settlements (Zhang et al.,2002, basata su Ishihara&Yoshimine, 1992)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if qc1ncs(ii)==0
                gamma_lim(ii)=0;
            else
                gamma_lim(ii)=max(0,(1.859*(2.163-0.478*(qc1ncs(ii)^0.264))^3));
            end
            
            if qc1ncs(ii)==0
                F_alfa(ii)=0;
            else
                F_alfa(ii)=-11.74+8.34*((max(69,qc1ncs(ii)))^0.264)-1.371*((max(69,qc1ncs(ii)))^0.528);
            end
            
            if Fs(ii)==0
                gamma_max(ii)=0;
            elseif Fs(ii)>2
                gamma_max(ii)=0;
            elseif Fs(ii)<=F_alfa(ii)
                gamma_max(ii)=gamma_lim(ii);
            else
                gamma_max(ii)=min(gamma_lim(ii),(0.035*(1-F_alfa(ii))*((2-Fs(ii))/(Fs(ii)-F_alfa(ii)))));
            end
            
            if qc1ncs(ii)==0
                epsV(ii)=0;
            else
                epsV(ii)=1.5*(exp(2.551-1.147*(max(21,qc1ncs(ii)))^0.264))*(min(0.08,gamma_max(ii)));
            end
            epsV(ii)=epsV(ii)*100;
            gamma_lim(ii)=gamma_lim(ii)*100;
            gamma_max(ii)=gamma_max(ii)*100;
            cedV(ii)=epsV(ii)*spessore(ii);
            
                      
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calcolo delle Espansioni Laterali con Zhang et al. (2004), come implementato come in Idriss & Boulanger (2008);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Delta_LDI(ii) = gamma_max(ii)*spessore(ii);
            
            % NOTA: Idriss e Boulanger (2008) suggeriscono estrema cautela
            % nell'adozione della procedura semplificata per passare da LDI a LD;
            % il software Cliq (http://www.geologismiki.gr/Products/CLiq.html) assume che LDI=LD.
            % Si sceglie di lasciare la stima di LD (come in Zhang, 2004) per
            % rendere confrontabili i risultati ottenuti con le altre procedure.
            
            if configurazione == 1
                Delta_LD(ii)=Delta_LDI(ii)*(conf_geom+0.2);
            else
                Delta_LD(ii)=6*(Delta_LDI(ii))*((conf_geom)^(-0.8));
            end
    end
 
    
    a=find(z(:)<=20);
    a_max = max(a);

    for jj=1:a_max
        w(jj) = 10-0.5*z(jj);
        wf(jj) = 1-0.05*z(jj);
        w_M(jj) = 25.56/z(jj);
        if w_M(jj)>25.56
            w_M(jj)=25.56;
        else
            w_M(jj)=w_M(jj);
        end
        m_M(jj)=exp(5/(25.56*(1-Fs(jj))))-1;
        
%         if Fs(jj)>1
%             F(jj) = 0;
%         elseif Fs(jj)==0
%             F(jj) = 0;
%         else
            F(jj) = 1/Fs(jj)*(0.3/Ic(jj));
%         end
        Delta_LPI(jj) = F(jj)*w_M(jj)*spessore(jj);
        Delta_LSI(jj) = PL(jj)*wf(jj)*spessore(jj);
        wfspessore(jj) = wf(jj)*spessore(jj);
    end
% Calcolo LSN non limitato ai primi 20m di profondità
for zz=1:length(z)
Delta_LSN(zz)=(epsV(zz)/(z(zz)*100)*spessore(zz));
end   
      
    %salvataggio file relativo all'amax(jj) analizzata
%         truncName = listing_CPT(kk).name(1:strfind(listing_CPT(kk).name,'.')-1);
        nome_out = [truncName,'_amax=',num2str(amax),'_M=',num2str(M),'_Idriss_Boulanger2014','.xlsx'];
        fileOut = [folderCPTout '/' folderCPTout1 '/' nome_out];
        Testata = {'z (mpc)','z (mslm)','qc (MPa)','fs (MPa)','u (MPa)','qt (MPa)','Rf(%)','gamma (kN/m3)','FC(%)','sigV_tot (kPa)','u0 (kpa)','sigV_eff (kPa)',...
            'Ic (-)','Cn','n','Qtn (MPa)','Fr(%)','rd','CSR','MSF','CSR,eq','Ksig','CSR*','qc1n','delta_qc1n','qc1ncs','CRR75','Fs(-)','ur/sigmaeff','PL(%)','epsv(%)','gamma_lim(%)','gamma_max(%)','cedv(cm)','spessore','Delta_LPI(-)','Delta_LSI(-)','Delta_LSN(-)','delta_LDI(cm)','delta_LD(cm)'};
        out=[z,zass,qc,f_s,u,qt,Rf,gamma,FC,sigV_t,u0,sigV_eff,Ic,Cn_calcolo,n_stress,Qtn,F_r,rd,CSR,MSF,CSR_eq,Ksig,CSR75,qc1n,delta_qc1n,qc1ncs,CRR75,Fs,deltaU,PL,epsV,gamma_lim,gamma_max,cedV,spessore,Delta_LPI,Delta_LSI,Delta_LSN,Delta_LDI,Delta_LD];
        xlswrite(fileOut,out,'Foglio1','A2');
        xlswrite(fileOut,Testata,'Foglio1','A1');

% close(bar)
% 
