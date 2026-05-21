# Pregled programa MAIN

Ta projekt je zdruzen robotski program. Glavna aplikacija povezuje:

- CyberGear motorje preko CAN
- DDSM115 pogonske motorje preko RS485
- orientacijske podatke iz BNO086
- UART debug izpis za spremljanje robotskih in stanjskih podatkov

## Struktura Programa

Osrednja logika programa je v `Core/Src/main.c`. Datoteka povezuje posamezne module, vsebuje glavno zanko, zbira merjene podatke in nastavlja ukaze za motorje.

Glavni programski moduli so:

- `main.c`: povezuje vse module in vsebuje glavno zanko
- `CyberGear.c`: funkcije za CAN komunikacijo s CyberGear motorji
- `DDSM115.c`: funkcije za ukaze DDSM115 motorjem preko RS485
- `DDSMove.c`: izracun poti in hitrosti iz DDSM enkoderjev
- `demo_app.c`: komunikacija z BNO086 senzorjem
- `Controllers.c`: zacetne funkcije za regulatorje

Program je razdeljen tako, da so merjenje, priprava podatkov, regulacija in komunikacija z motorji loceni. Za razumevanje delovanja so najpomembnejsi podatkovni tabeli `robot_data` in `state_data` ter izhodne spremenljivke, ki dolocajo cilje motorjev.

## Potek Programa

Po zagonu program pripravi robotske module in motorje nastavi v delovne nacine. CyberGear motorji se pocistijo, nastavijo, postavijo na mehansko niclo in omogocijo. DDSM115 pogonska motorja se nastavita v tokovni nacin in zacneta z nicelnim tokom.

Glavna zanka se nato izvaja neprekinjeno:

1. `BNO_App(robot_data)` prebere nove podatke iz BNO086, ko prekinitev sporoci, da so na voljo novi podatki.
2. `RobotData_Update()` prekopira zadnje podatke motorjev, koles, razdalje, hitrosti in BNO v `robot_data`.
3. `StateData_Update()` ustvari strnjen vektor stanja:
   - pozicija
   - hitrost
   - kot telesa
   - kotna hitrost
4. Vsakih `DEBUG_PRINT_INTERVAL_MS` se vse vrednosti `robot_data` in `state_data` izpisejo preko UART.
5. Ce je bil pritisnjen modri gumb, se spremenijo zeleni cilji motorjev.
6. `DDSM_Service()` poslje en DDSM ukaz naenkrat in pocaka na povratni podatek ali timeout.
7. `CyberGear_Service()` posilja pozicijske ukaze CyberGear motorjem, po en motor naenkrat.

## Naslovi in ID-ji

Program uporablja te logicne naslove in ID-je motorjev:

- BNO086 I2C naslov: `0x4B` (`BNO_ADDR` je definiran kot `0x4B << 1` za uporabo s STM32 HAL)
- DDSM levi pogonski motor: `0x10`
- DDSM desni pogonski motor: `0x30`
- CyberGear host ID: `0xFE`
- CyberGear motor 1: `17`
- CyberGear motor 2: `18`
- CyberGear motor 3: `19`
- CyberGear motor 4: `20`

## Service Funkcije

`DDSM_Service()` skrbi za RS485 pogonska motorja v preprostem neblokirnem zaporedju. Funkcija poslje tokovni ukaz enemu DDSM motorju, ga oznaci kot cakajocega in pocaka, da UART callback potrdi odgovor motorja. Ce povratni podatek ne pride znotraj `DDSM_FEEDBACK_TIMEOUT_MS`, funkcija pobrise cakajoce stanje, poveca stevec timeoutov, ponovno zazene UART sprejem in nadaljuje z naslednjim motorjem. Tako se ukazi ne posiljajo obema DDSM motorjema prehitro, povratni podatek pa je lazje povezati z motorjem, ki je dobil ukaz.

`CyberGear_Service()` skrbi za stiri CyberGear motorje preko CAN. Vsakic, ko se izvede, najprej preveri, ali je prost CAN transmit mailbox. Ce je prost, poslje naslednji pozicijski ukaz z `SetAngle()`. Funkcija krozi cez motor 1, motor 2, motor 3 in motor 4, zato se vsi stirje motorji posodabljajo zaporedno, brez poskusa posiljanja vseh CAN okvirjev hkrati.

## Robot Data

`robot_data` je tabela z 20 vrednostmi za spremljanje in izvoz podatkov. Sama neposredno ne krmili robota. Hrani:

- levo/desno DDSM razdaljo in hitrost
- pozicije CyberGear motorjev
- pozicije DDSM koles
- BNO roll, pitch, yaw
- BNO gyro X/Y/Z
- BNO quaternion W/X/Y/Z

## State Data

`state_data` je manjsi vektor s 4 vrednostmi, namenjen za stanje robota oziroma regulacijo:

- `state_data[0]`: povprecna pozicija robota
- `state_data[1]`: povprecna hitrost robota
- `state_data[2]`: kot telesa robota iz BNO roll
- `state_data[3]`: kotna hitrost robota iz BNO gyro X

## Vrednosti za Regulacijo

Za nadaljnjo regulacijo ali delo s krmilniki so najbolj uporabne vrednosti ze zbrane v `state_data` in `robot_data`.

`state_data` predstavlja strnjen vhodni vektor za regulator:

- pozicija: `state_data[0]`
- hitrost: `state_data[1]`
- kot telesa: `state_data[2]`
- kotna hitrost: `state_data[3]`

`robot_data` vsebuje podrobnejse podatke za primere, kjer regulator potrebuje locene hitrosti levega/desnega kolesa, posamezne pozicije CyberGear motorjev ali celotne BNO orientacijske podatke.

Glavni krmilni izhodi, ki so na voljo v `main.c`, so:

- `desired_angle`, `desired_angle2`, `desired_angle3`, `desired_angle4` za ciljne pozicije CyberGear motorjev
- `ddsm_left_current_cmd` in `ddsm_right_current_cmd` za tokovna ukaza DDSM pogonskih motorjev

Regulator lahko iz `state_data` izracuna ciljne vrednosti in jih zapise v zgornje izhodne spremenljivke. Funkciji `DDSM_Service()` in `CyberGear_Service()` nato poskrbita za dejansko posiljanje ukazov motorjem. V projektu je prisotna tudi koda za regulatorje v `Controllers.c`, med drugim LQR in lateralni regulator.

## Obmocja za Nadaljnji Razvoj

Za razvoj regulacije je smiselno spreminjati predvsem izracun ciljnih vrednosti, medtem ko komunikacijske funkcije ostanejo nespremenjene.

Tipicna obmocja sprememb:

- `desired_angle...` za premik CyberGear motorjev
- `ddsm_left_current_cmd` in `ddsm_right_current_cmd` za pogon DDSM motorjev
- izracun v regulatorju, ki nastavi zgornje vrednosti
- interval debug izpisa `DEBUG_PRINT_INTERVAL_MS`

Obmocja, ki zahtevajo posebno pozornost:

- ID-je motorjev, ker morajo ustrezati dejanskim motorjem
- nacin motorjev, ker mora ukaz ustrezati izbranemu nacinu
- timeout vrednosti, ker vplivajo na zanesljivost komunikacije

Priporocen pristop je, da regulator bere `state_data`, izracuna ukaze in jih zapise v obstojece izhodne spremenljivke. Na ta nacin ostane komunikacija z motorji locena od regulacijskega algoritma.

## Podatkovni Tok

Delovanje programa sledi temu podatkovnemu toku:

1. Senzorji in motorji posljejo povratne podatke.
2. Callback funkcije in BNO funkcija posodobijo zadnje izmerjene vrednosti.
3. `RobotData_Update()` zbere podrobne podatke v `robot_data`.
4. `StateData_Update()` naredi kratek vektor `state_data` za regulacijo.
5. Regulator ali testna logika nastavi ciljne vrednosti motorjev.
6. `DDSM_Service()` in `CyberGear_Service()` posljeta ukaze motorjem.

Ta locitev zmanjsa odvisnost med merjenjem, regulacijo in komunikacijo z motorji ter olajsa nadaljnje razsirjanje programa.

## Obdelava Povratnih Podatkov

DDSM povratni podatki se sprejemajo preko UART callback funkcije. Ko pride veljaven paket levega ali desnega motorja, program posodobi:

- kot kolesa
- RPM kolesa
- tok
- radialno hitrost
- prevozeno razdaljo
- translacijsko hitrost

CyberGear povratne podatke obdeluje CAN receive callback v `CyberGear.c`, kjer dohodna CAN sporocila posodobijo zadnje vrednosti motorjev.

## Obnasanje Gumba

Zunanja prekinitev gumba samo nastavi zastavico za zahtevo. Glavna zanka nato prebere to zastavico in spremeni zelene kote CyberGear motorjev ter tokovne ukaze DDSM motorjev. Tako se ukazi motorjem ne izvajajo neposredno znotraj prekinitve, ampak v glavni zanki.

## Debug Izpis

UART debug vrstica se zacne z `RD` in izpise tako celotne robotske podatke kot tudi strnjene podatke stanja. To je uporabno za logiranje, risanje grafov ali preverjanje, ali se senzorji in motorni povratni podatki pravilno posodabljajo.
