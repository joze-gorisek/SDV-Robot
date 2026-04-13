clc; clear; close all;                              % počisti Command Window, izbriše spremenljivke, zapre okna 

%{
ZICNA KOMUNIKACIJA

UPORABLJAMO:
Nivojska predstavitev: NRZ
Zaščitno kodiranje: CRC-16
Kanalsko kodiranje: Hamming
%}
 
%https://drmoazzam.com/matlab-code-bpsk-modulation-and-demodulation-with-explanation
%https://dsp.stackexchange.com/questions/22316/simulating-bpsk-costas-loop-in-matlab

%NASTAVITVE
global T R U_max frekvenca vzorci

T = 0.1;       %Cas posiljanja enega bita
R = 1/T;        % [bit/s] bitna hitrost
U_max = 1;      % [V] maksimalna generirana napetost
frekvenca = 20; % [Hz] frekvenca nosilnega signala
vzorci = 20;   % Stevilo vzorcev za vsak poslan bit. To potrebujemo za simuliranje realne zvezne napetosti




fprintf('\n==========================================\n');
fprintf('   NASTAVITVE KOMUNIKACIJSKEGA SISTEMA   \n');
fprintf('==========================================\n');
fprintf('Trajanje enega bita (T):      %.3f s\n', T);
fprintf('Bitna hitrost (R):            %.0f bit/s\n', R);
fprintf('Maksimalna napetost (U_max):  %.1f V\n', U_max);
fprintf('Frekvenca nosilnega signala (f):        %.0f Hz\n', frekvenca);
fprintf('Število vzorcev na bit:       %.0f\n', vzorci);
fprintf('Frekvenca vzorčenja (fs):     %.0f Hz\n', vzorci/T);
fprintf('==========================================\n\n');




%% Funkcija oddajnika
function tx_signal = oddajnik(sporocilo, naslov_odd, naslov_spr, zamik_faze, brezzicno)
    global U_max vzorci encoded T frekvenca
    disp('Pošiljanje sporočila')

    %Ta oddajnik lahko oddaja samo sporočila v tekstu, zato je že vgrajen
    %tip sporočila na tekst.

    %% PRETVORBA BESEDILA V BITE
    sporocilo_ascii = uint8(sporocilo);                % double pretvori znake v ASCII kode in zapise v novo spremenljivko
    sporocilo_bin = dec2bin(sporocilo_ascii,8)';        % vsak znak v pretvori v 8-bitno binarno obliko
    temp = sporocilo_bin(:)';                           % združi vse bite v en niz
    sporocilo_bits = temp - '0';                        % pretvori '0','1' v 0,1. Deluje da od vsakega znaka '0' vrednosti 48 ali '1' vrednosti 49 odšteje znak '0' vrednosti 48, tako da na koncu ostane 0 ali 1
                                                        % Lahko bi nadomestili
                                                        % s sporocilo_bits = (temp == '1')
 
    %% DODAJANJE PREAMBULE
    preambula = [1 0 1 0 1 0 1 0];                      % preambula dolžine 8 bitov (točno definirana)
    
    %% DODAJANJE ZAČETKA OKVIRJA
    SFD = [1 0 1 0 1 0 1 1];          % Start Frame Delimiter (začetek okvirja)

    %% DODAJANJE KONTROLNIH ZNAKOV (tip okvirja, naslov pošiljatelja, naslov sprejemnika)
    tip_okvirja = [0 0 0 1]; %tip okvirja bo v 4 bitih. 1 pomeni tekst. Ta oddajnik omogoča pošiljanje samo teksta.

    % naslov_odd je naslov pošiljatelja
    oddajnik_bin = dec2bin(naslov_odd,4)';        
    temp = oddajnik_bin(:)';                           
    oddajnik_bits = temp - '0';

    % naslov_spr je naslov sprejemnika
    sprejemnik_bin = dec2bin(naslov_spr,4)';        
    temp = sprejemnik_bin(:)';                           
    sprejemnik_bits = temp - '0';

    %% DODAJANJE KONCA OKVIRJA
    EFD = [1 1 1 1 1 1 1 1];        % End frame delimiter (konec okvirja)

    %% ZDRUŽEVANJE V EN OKVIR
    okvir = [preambula SFD tip_okvirja oddajnik_bits sprejemnik_bits sporocilo_bits EFD];


    %% CRC-16 ZAŠČITNO KODIRANJE
    polinom = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];      % standardni polinom za CRC-16-CCITT
    podatki_crc = [okvir zeros(1,16)];                  % dodamo vektor s 16 ničlami za deljenje (ostanek)
    
    temp = podatki_crc;                                 % kopija temp za računanje
    
    for i = 1:length(okvir)                             % zanka dolžine okvirja
        if temp(i) == 1                                 % če je bit 1 izvedemo XOR
            temp(i:i+16) = xor(temp(i:i+16), polinom);  % XOR z generatorjem
        end
    end
    
    crc = temp(end-15:end);                             % ostanek se zapise v zadnjih 16 bitov. Ostanek = CRC.
    okvir_crc = [okvir crc];                            % dodamo CRC na konec.
                                                                                   
    %% HAMMING (7,4) KANALSKO KODIRANJE
    data = okvir_crc;
    n = 7; %Skupno število bitov
    k = 4; %Število podatkovnih bitov
    encoded = encode(data, n, k, 'hamming/binary');     % pridobljeno iz https://www.mathworks.com/help/comm/ref/encode.html

    %% NRZ NIVOJSKA REDSTAVITEV
    NRZ_nivoji = 2*U_max*encoded - U_max;               % Sprememba iz bitov v napetstne nivoje Non Return to zero => 0 -> -u_max, 1-> U_max 
    NRZ_signal = rectpulse(NRZ_nivoji, vzorci);         % Razteg novijev na signal kjer je vsak bit predstavljen z večimi vzorci  

    if brezzicno
        %% MODULACIJA z BPSK
        % Tukaj signal pomnožimo z osnovno sinusoido.
        korak = T/vzorci;
        st_bitov = length(NRZ_nivoji);
        skupni_cas = st_bitov * T;
        t = linspace(korak, skupni_cas, length(NRZ_signal));

        nosilni_signal = sqrt(2/T)*cos(2*pi*frekvenca*t + zamik_faze);
        tx_signal = NRZ_signal.* nosilni_signal;
        disp('Sporočilo poslano brezžično')
    else
        tx_signal = NRZ_signal;
        disp('Sporočilo poslano žično')
    end
end

%% Funkcija za dodajanje šuma. Vhod je signal in SNR [dB]. Izhod je popačen signal
function tx_signal_sum = dodajanje_suma(tx_signal, SNR)
    global vzorci 
    disp('Dodajanje šuma')
    % Uporabljamo AWGN (Additive White Gaussian Noise)
    % BER (Bit Error Rate) = napacni bit / skupno poslanih bitov
    % SNR_lin (Signal to Noise Ratio) = moc signala / moc suma. 
    %  SNR bita [dB]
    SNR_lin = 10^(SNR / 10);
    SNR_vzorca = SNR-10*log10(vzorci);                  % SNR enega vzorca [dB]                                                    
    tx_signal_sum = awgn(tx_signal,SNR_vzorca);
    
    disp('Šum dodan')
    fprintf('SNR: %.1f dB\n', SNR);
    fprintf('Linearen SNR: %f\n\n', SNR_lin)
    
end

%% Funkcija srejemnika. Vhod je signal in spremenljivka, ki pove ali imamo brezžično komunikacijo
function sporocilo = sprejemnik(rx_signal , naslov_spr,brezzicno)
    global vzorci encoded T frekvenca
    disp('Prejemanje sporočila')
    if brezzicno
        %% DEMODULACIJA z BPSK in COSTASOVO zanko
        korak = T/vzorci;
        st_vzorcev = length(rx_signal);
        skupni_cas = korak * st_vzorcev;
        t = linspace(korak, skupni_cas, st_vzorcev);

        mu = 0.0001; %hitrost prilagajanja zanke
        faza = 0; %zacetni ocenjen fazni zamik

        I_out = zeros(size(rx_signal)); %kosinusna veja
        Q_out = zeros(size(rx_signal)); %sinusna veja. Hočemo da je ta 0
        napaka = zeros(size(rx_signal));

        for n = 1:length(rx_signal)
            % Lokalna nosilca s trenutno oceno faze
            cos_nosilec = cos(2*pi*frekvenca*t(n) + faza);
            sin_nosilec = sin(2*pi*frekvenca*t(n) + faza);
            
            % Demodulacija BPSK
            I_out(n) = rx_signal(n) * cos_nosilec;
            Q_out(n) = rx_signal(n) * sin_nosilec;

            napaka(n) = I_out(n) * Q_out(n);

            faza = faza - mu * napaka(n);
        end

        demodulirano=sqrt(2/T)*I_out;   % Uporabimo za normalizacijo energije bita. Eb = 1
        rx_povprecja = [];
        for i=1:vzorci:length(demodulirano)
            rx_povprecja=[rx_povprecja trapz(t(i:i+vzorci-1),demodulirano(i:i+vzorci-1))]; % trapz pogleda ploščino pod integralom
        end
    else
        %% Prejem signala v obliki napetosti. Uporabljamo NRZ:
        rx_matrika = reshape(rx_signal, vzorci, []);        % Pretvorimo signal v matriko ki ima n vrstic. Vsak stolpec tako predstavlja "en bit" ki je vzorčen na n vzorcev. 
        rx_povprecja = mean(rx_matrika);       % To uporabimo da izracunamo povprecje napetosti. Tako se najbolje lahko izognemo temu da sum pokvari sporocilo
    end
    
    rx_bits = rx_povprecja > 0;                            % ce so povprecja vec kot 0V to predstavlja bit 1 (Logicna 1). 
    
    %% BER IZRAČUN
    napake = sum(rx_bits ~= encoded);
    BER = napake/length(rx_bits);
    fprintf('BER: %f\n', BER);
    
    
    %% HAMMING DEKODIRANJE. napako zazna in popravi
    n = 7; %Skupno število bitov
    k = 4; %Število podatkovnih bitov
    decoded = decode(rx_bits, n, k, 'hamming/binary'); % enako kot encode. Pridobljeno iz https://www.mathworks.com/help/comm/ref/decode.html?s_tid=srchtitle_support_results_1_decode
    
    %% PREVERIMO PREAMBULO
    % Ker imamo v istem programu ne sinhroniziramo takta, kot bi ga v
    % realnem sprejemniku.
    %Če je preambula narobe je verjetno obrnjena. Potem obrnemo celotno
    %sporočilo

    preambula = [1 0 1 0 1 0 1 0];

    if preambula == ~decoded(1:length(preambula))
        decoded = ~decoded;
    end
    
    %% CRC PREVERJANJE. Zazna napako
    polinom = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];
    received = decoded;
    
    temp = received;
    
    for i = 1:length(received)-16
        if temp(i) == 1
            temp(i:i+16) = xor(temp(i:i+16), polinom);
        end
    end
    
    crc_check = temp(end-15:end);                       % preverimo ostanek
    
    if sum(crc_check) == 0
        disp('CRC OK - ni napak');
    else
        disp('CRC NAPAKA');
    end
    
    %% ISKANJE ZAČETKA OKVIRJA
    SFD = [1 0 1 0 1 0 1 1];
    
    % Poiščemo, kje se v received pojavi SFD vzorec
    idx = strfind(received, SFD);
    
    if isempty(idx)
        disp('SFD ni bil najden. Sinhronizacija ni uspela.');
        return
    end
    
    start_ptr = idx(1) + length(SFD);

    %% KONTROLNI ZNAKI
    sporocilo = "";
    % Tip okvirja (4 biti)
    tip_bits = received(start_ptr : start_ptr + 3);
    start_ptr = start_ptr + 4;

    if tip_bits ~= [0 0 0 1]
        disp("Sporočilo ni tekst. Ta sprejemnik lahko sprejema le tekst.")
        return
    end
    
    % Naslov pošiljatelja (4 biti)
    posiljatelj_bits = received(start_ptr : start_ptr + 3);
    id_posiljatelja = bin2dec(char(posiljatelj_bits + '0'));
    start_ptr = start_ptr + 4;
    
    % Naslov sprejemnika (4 biti)
    prejemnik_bits = received(start_ptr : start_ptr + 3);
    id_prejemnika = bin2dec(char(prejemnik_bits + '0'));
    start_ptr = start_ptr + 4;
    
    fprintf('Prejet paket od ID: %d, namenjen ID: %d\n', id_posiljatelja, id_prejemnika);

    if id_prejemnika ~= naslov_spr
        fprintf('Paket ni namenjen meni (ID: %d). Prekinjam branje.', naslov_spr);
        return
    end

    %% ISKANJE KONCA OKVIRJA (EFD)
    efd_pattern = [1 1 1 1 1 1 1 1];
    idx_end = strfind(received(start_ptr:end), efd_pattern);
    
    if isempty(idx_end)
        % Če ni EFD, vzamemo do konca brez CRC-ja
        konec_ptr = length(received) - 16; % Odstranimo CRC
    else
        % Če najde EFD, vzame bite do tja
        konec_ptr = start_ptr + idx_end(1) - 1; 
    end
    
    %% IZLUŠČITEV SPOROČILA
    sporocilo_bits = received(start_ptr : konec_ptr); %Vse skupaj odstranimo

    
    %% PRETVORBA NAZAJ V BESEDILO
    sporocilo_bits = sporocilo_bits(1:floor(length(sporocilo_bits)/8)*8); % poravnava
    
    chars = reshape(sporocilo_bits,8,[])';               % razdeli v 8-bitne bloke
    ascii_vals = bin2dec(char(chars + '0'));            % pretvori v ASCII
    sporocilo = char(ascii_vals);               % končni tekst

end

%% ------------------------------------------------------------------------------------------------------------
%% TEST DELOVANJA
%% ------------------------------------------------------------------------------------------------------------

sporocilo = 'Ales Novak E300400500 Slovenija';      % besedilo za pošiljanje
SNR = 11;
brezzicno = true;
zamik_faze = 6*pi/10;
oddajnik_id = 1;
sprejemnik_id = 10;


tx_signal = oddajnik(sporocilo, oddajnik_id, sprejemnik_id, zamik_faze, brezzicno);
tx_signal = dodajanje_suma(tx_signal, SNR);
prejeto_sporocilo = sprejemnik(tx_signal, sprejemnik_id, brezzicno);

fprintf('\n%s\n', "Prejeto sporocilo:");             % nova vrstica in izpis besedila
fprintf('%s\n', prejeto_sporocilo);                 % izpis sporocila
