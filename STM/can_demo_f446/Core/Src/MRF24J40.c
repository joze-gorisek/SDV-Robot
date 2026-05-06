/*
 * MRF24J40.c
 *
 *  Created on: Mar 28, 2024
 *      Author: Benjamin Dvorsak
 */

#include"MRF24J40.h"
#include"main.h"
#include"stdbool.h"
#include"MRF24J40_includes.h"
#include "string.h"
#include "stdio.h"

extern SPI_HandleTypeDef hspi2;
extern UART_HandleTypeDef huart2;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// aMaxPHYPacketSize = 127, from the 802.15.4-2006 standard.
 uint8_t rx_buf[127];

// essential for obtaining the data frame only
// bytes_MHR = 2 Frame control + 1 sequence number + 2 panid + 2 shortAddr Destination + 2 shortAddr Source
const int bytes_MHR = 9;
const int bytes_FCS = 2; // FCS length = 2
int bytes_nodata = bytes_MHR + bytes_FCS; // no_data bytes in PHY payload,  header length + FCS
//bytes_nodata = bytes_MHR + bytes_FCS; // no_data bytes in PHY payload,  header length + FCS

int ignoreBytes = 0; // bytes to ignore, some modules behaviour.

bool bufPHY = true; // flag to buffer all bytes in PHY Payload, or not

volatile uint8_t flag_got_rx;
volatile uint8_t flag_got_tx;

rx_info_t rx_info;
tx_info_t tx_info;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Mrf24j_init(void)
{
	uint8_t check;

	// from mrf24j40 datasheet p90 3.2 init
	Mrf24j_write_short(MRF_SOFTRST, 0x7);

	while ((Mrf24j_read_short(MRF_SOFTRST) & 0x7) != 0) {
	; // wait for soft reset to finish
	}
	HAL_Delay(10);
//	Mrf24j_write_short(MRF_PACON2, 0x98); // – Initialize FIFOEN = 1 and TXONTS = 0x6.
	Mrf24j_write_short(MRF_TXSTBL, 0x95); // – Initialize RFSTBL = 0x9.

	// wait for mrf to be in receive mode

//		do {
//			check = Mrf24j_read_long(RFSTATE);
//		} while (check & 0xa0 != 0xa0);

	Mrf24j_address16_write(0x0011);
	Mrf24j_set_pan(0x2222);

	Mrf24j_write_long(MRF_RFCON0, 0x03); // – Initialize RFOPT = 0x03.
	Mrf24j_write_long(MRF_RFCON1, 0x02); // – Initialize VCOOPT = 0x02.
	Mrf24j_write_long(MRF_RFCON2, 0x80); // – Enable PLL (PLLEN = 1).
	Mrf24j_write_long(MRF_RFCON3, 0x00); // – set tx max power
	Mrf24j_write_long(MRF_RFCON6, 0x90); // – Initialize TXFIL = 1 and 20MRECVR = 1.
	Mrf24j_write_long(MRF_RFCON7, 0x80); // – Initialize SLPCLKSEL = 0x2 (100 kHz Internal oscillator).
	Mrf24j_write_long(MRF_RFCON8, 0x10); // – Initialize RFVCO = 1.
//	Mrf24j_read_long(MRF_SLPCON0);
	Mrf24j_write_long(MRF_SLPCON0, 0x02); //interrupt polarity - rising edge
//	Mrf24j_write_long(MRF_SLPCON1, 0x21); // – Initialize CLKOUTEN = 1 and SLPCLKDIV = 0x01.
	HAL_Delay(10);
	//  Configuration for nonbeacon-enabled devices (see Section 3.8 “Beacon-Enabled and
	//  Nonbeacon-Enabled Networks”):
	Mrf24j_write_short(MRF_BBREG2, 0x80); // Set CCA mode to ED
	Mrf24j_write_short(MRF_BBREG6, 0x40); // – Set appended RSSI value to RXFIFO.
	HAL_Delay(10);
	Mrf24j_write_short(MRF_CCAEDTH, 0x60); // – Set CCA ED threshold.


	// wait for mrf to be in receive mode

	do {
		check = Mrf24j_read_long(RFSTATE);
	} while (check & 0xa0 != 0xa0);

	Mrf24j_write_short(MRF_RXMCR, 0x01);

	Mrf24j_set_interrupts();
	//
	Mrf24j_set_channel(20);
	// max power is by default.. just leave it...
	// Set transmitter power - See “REGISTER 2-62: RF CONTROL 3 REGISTER (ADDRESS: 0x203)”.
	Mrf24j_write_short(MRF_RFCTL, 0x04); //  – Reset RF state machine.
	Mrf24j_write_short(MRF_RFCTL, 0x00); // part 2
	HAL_Delay(100);
	Mrf24j_rx_flush();
	//delay(1); // delay at least 192usec

	HAL_Delay(1);
}


void Mrf24j_write_short(uint8_t address, uint8_t data)
{
	  HAL_StatusTypeDef test = 0;
	  HAL_GPIO_WritePin(SPI2_CS_MRF_GPIO_Port, SPI2_CS_MRF_Pin, GPIO_PIN_RESET);

	  address = (address<<1 & 0b01111110) | 0x01;
	  test = HAL_SPI_Transmit(&hspi2, &address, 1,1);
	  if(test != HAL_OK)
	  {
		  HAL_UART_Transmit (&huart2, (uint8_t*)"napaka1\r\n", strlen("napaka1\r\n"), HAL_MAX_DELAY);
	  }
	  else
	  {
		  HAL_UART_Transmit (&huart2, (uint8_t*)"ok1\r\n", strlen("ok1\r\n"), HAL_MAX_DELAY);
	  }
	  test = HAL_SPI_Transmit(&hspi2, &data, 1,1);
	  if(test != HAL_OK)
	  {
		  HAL_UART_Transmit (&huart2, (uint8_t*)"napaka2\r\n", strlen("napaka2\r\n"), HAL_MAX_DELAY);
		  HAL_UART_Transmit (&huart2, &test, 1, HAL_MAX_DELAY);
	  }
	  else
	  {
		  HAL_UART_Transmit (&huart2, (uint8_t*)"ok2\r\n", strlen("ok2\r\n"), HAL_MAX_DELAY);
	  }

	  HAL_GPIO_WritePin(SPI2_CS_MRF_GPIO_Port, SPI2_CS_MRF_Pin, GPIO_PIN_SET);

	}
void Mrf24j_write_long(uint16_t address, uint8_t data)
{
	  HAL_GPIO_WritePin(SPI2_CS_MRF_GPIO_Port, SPI2_CS_MRF_Pin, GPIO_PIN_RESET);

	  uint8_t ahigh = 0x80 | (address >> 3);
	  uint8_t alow = (address << 5) | 0x10;

	  HAL_SPI_Transmit(&hspi2, &ahigh, 1,1);
	  HAL_SPI_Transmit(&hspi2, &alow, 1,1);

	  HAL_SPI_Transmit(&hspi2, &data, 1,1);

	  HAL_GPIO_WritePin(SPI2_CS_MRF_GPIO_Port, SPI2_CS_MRF_Pin, GPIO_PIN_SET);

	}

uint8_t Mrf24j_read_short(uint8_t address)
{

	  HAL_GPIO_WritePin(SPI2_CS_MRF_GPIO_Port, SPI2_CS_MRF_Pin, GPIO_PIN_RESET);

	  address = address<<1 & 0b01111110;
	  HAL_SPI_Transmit(&hspi2, &address, 1,1);

	  uint8_t ret = 0x0;
	  uint8_t bogus = 0x00;
	  HAL_SPI_TransmitReceive(&hspi2, &bogus, &ret, 1, 1);
	//  HAL_SPI_Receive(spi, &ret,1,1);

	  HAL_GPIO_WritePin(SPI2_CS_MRF_GPIO_Port, SPI2_CS_MRF_Pin, GPIO_PIN_SET);

	  return ret;
	}

uint8_t Mrf24j_read_long(uint16_t address)
{
	  HAL_GPIO_WritePin(SPI2_CS_MRF_GPIO_Port, SPI2_CS_MRF_Pin, GPIO_PIN_RESET);

	  uint8_t ahigh = 0x80 | (address >> 3);
	  uint8_t alow = address << 5;

	  HAL_SPI_Transmit(&hspi2, &ahigh, 1,1);
	  HAL_SPI_Transmit(&hspi2, &alow, 1,1);

	  uint8_t ret = 0x0;
	  uint8_t bogus = 0x00;
	  HAL_SPI_TransmitReceive(&hspi2, &bogus, &ret, 1, 1);
	  //HAL_SPI_Receive(spi, &ret,1,1);

	  HAL_GPIO_WritePin(SPI2_CS_MRF_GPIO_Port, SPI2_CS_MRF_Pin, GPIO_PIN_SET);

	  return ret;
	}
void Mrf24j_set_interrupts(void)
{
	uint8_t init_state = 0xFF;
	init_state &= ~0x8;

  // interrupts for rx and tx normal complete
	Mrf24j_write_short(MRF_INTCON, init_state);//~0b11110110);
}
void Mrf24j_set_channel(uint8_t channel)
{
  //  (((channel - 11) << 4) | 0x03));
	Mrf24j_write_long(MRF_RFCON0, (((channel - 11) << 4) | 0x03));
}


void Mrf24j_reset(void)
{
  HAL_GPIO_WritePin(MRF_RESET_GPIO_Port, MRF_RESET_Pin, GPIO_PIN_RESET);
  HAL_Delay(100);
  HAL_GPIO_WritePin(MRF_RESET_GPIO_Port, MRF_RESET_Pin, GPIO_PIN_SET);
  HAL_Delay(200);
}

void Mrf24j_send16(uint16_t dest16, char * data, uint8_t len) {
//  //byte len = strlen(data); // get the length of the char* array
//  int i = 0;
//  write_long(i++, bytes_MHR); // header length
//  // +ignoreBytes is because some module seems to ignore 2 bytes after the header?!.
//  // default: ignoreBytes = 0;
//  write_long(i++, bytes_MHR+ignoreBytes+len);
//
//  // 0 | pan compression | ack | no security | no data pending | data frame[3 bits]
//  write_long(i++, 0b01100001); // first byte of Frame Control
//  // 16 bit source, 802.15.4 (2003), 16 bit dest,
//  write_long(i++, 0b10001000); // second byte of frame control
//  write_long(i++, 1);  // sequence number 1
//
//  word panid = get_pan();
//
//  write_long(i++, panid & 0xff);  // dest panid
//  write_long(i++, panid >> 8);
//  write_long(i++, dest16 & 0xff);  // dest16 low
//  write_long(i++, dest16 >> 8); // dest16 high
//
//  word src16 = address16_read();
//  write_long(i++, src16 & 0xff); // src16 low
//  write_long(i++, src16 >> 8); // src16 high
//
//  // All testing seems to indicate that the next two bytes are ignored.
//  //2 bytes on FCS appended by TXMAC
//  i+=ignoreBytes;
//  for (int q = 0; q < len; q++) {
//    write_long(i++, data[q]);
//  }
//  // ack on, and go!
//  write_short(MRF_TXNCON, (1<<MRF_TXNACKREQ | 1<<MRF_TXNTRIG));

  //byte len = strlen(data); // get the length of the char* array
  int i = 0;
  static int seq_num=1;
  Mrf24j_write_long(i++, bytes_MHR); // header length
  // +ignoreBytes is because some module seems to ignore 2 bytes after the header?!.
  // default: ignoreBytes = 0;


  //bytes_MHR = 2 Frame control + 1 sequence number + 2 panid + 2 shortAddr Destination + 2 shortAddr Source
//  bytes_MHR = 3+2+2+2;
  //  bytes_MHR  = 3+8+8+2; //extended address
  Mrf24j_write_long(i++, bytes_MHR+ignoreBytes+len);

  // 0 | pan compression | ack | no security | no data pending | data frame[3 bits]
  Mrf24j_write_long(i++, 0x61); // first byte of Frame type
  // 16 bit source, 802.15.4 (2003), 16 bit dest,
  //Mrf24j_write_long(i++, 0x8c); // second byte of frame control 0x8c 4 bitni source 8 dest
  Mrf24j_write_long(i++, 0x88); //  0b10001000;	// short dest (10,11) short src (14,15)

  Mrf24j_write_long(i++, seq_num++);  // sequence number 1

  uint16_t panid =0x3322;//= Mrf24j_get_pan();

  Mrf24j_write_long(i++, panid & 0xff);  // dest panid
  Mrf24j_write_long(i++, panid >> 8);

  // Dest
  //dest16=0x5511;
  Mrf24j_write_long(i++, dest16 & 0xff);  // dest16 low
  Mrf24j_write_long(i++, dest16 >> 8); // dest16 high

  // Source
  uint16_t src16 = 0x4455; //Mrf24j_address16_read();

  Mrf24j_write_long(i++, src16 & 0xff); // src16 low
  Mrf24j_write_long(i++, src16 >> 8); // src16 high

  // All testing seems to indicate that the next two bytes are ignored.
  //2 bytes on FCS appended by TXMAC
  i+=ignoreBytes;
  for (int q = 0; q < len; q++) {
	  Mrf24j_write_long(i++, data[q]);
  }
  // ack on, and go!
  Mrf24j_write_short(MRF_TXNCON, (1<<MRF_TXNACKREQ | 1<<MRF_TXNTRIG));
}


uint16_t Mrf24j_get_pan(void) {
  uint8_t panh = Mrf24j_read_short(MRF_PANIDH);
  return (panh << 8) | Mrf24j_read_short(MRF_PANIDL);
}

void Mrf24j_set_pan(uint16_t panid) {
	Mrf24j_write_short(MRF_PANIDH, panid >> 8);
	Mrf24j_write_short(MRF_PANIDL, panid & 0xff);
}

void Mrf24j_address16_write(uint16_t address16) {
	Mrf24j_write_short(MRF_SADRH, address16 >> 8);
	Mrf24j_write_short(MRF_SADRL, address16 & 0xff);
}

uint16_t Mrf24j_address16_read(void) {
  uint8_t a16h = Mrf24j_read_short(MRF_SADRH);
  return a16h << 8 | Mrf24j_read_short(MRF_SADRL);
}

void Mrf24j_rx_flush(void) {
	uint8_t rxflush;

	rxflush = Mrf24j_read_short(RXFLUSH);
	rxflush = (rxflush | 0b00000001);
	rxflush = 0x05;
	Mrf24j_write_short(RXFLUSH, rxflush);
}

void Mrf24j_rx_disable(void){
	Mrf24j_write_short(MRF_BBREG1, 0x04);  // RXDECINV - disable receiver
}

void Mrf24j_rx_enable(void) {
	Mrf24j_write_short(MRF_BBREG1, 0x00);  // RXDECINV - enable receiver
}

/**
 * Call this from within an interrupt handler connected to the MRFs output
 * interrupt pin.  It handles reading in any data from the module, and letting it
 * continue working.
 * Only the most recent data is ever kept.
 */
void Mrf24j_interrupt_handler(void) {
  uint8_t last_interrupt = Mrf24j_read_short(MRF_INTSTAT);
  Mrf24j_rx_flush();
  if (last_interrupt & MRF_I_RXIF) {
    flag_got_rx++;
    // read out the packet data...
    //noInterrupts();
    Mrf24j_rx_disable();
    // read start of rxfifo for, has 2 bytes more added by FCS. frame_length = m + n + 2
    uint8_t frame_length = Mrf24j_read_long(0x300);
    for(int i = 0;i<100;i++);
    Mrf24j_rx_enable();
//    HAL_Delay(100);


    if(frame_length < 128)
    {


    // buffer all bytes in PHY Payload
    if(bufPHY){
      int rb_ptr = 0;
      for (int i = 1; i < frame_length+2; i++) { // from 0x301 to (0x301 + frame_length -1)

        rx_buf[rb_ptr++] = Mrf24j_read_long(0x300 + i);
      }
    }
    } else{

    	Mrf24j_write_short(MRF_RXFLUSH, 0x05);

    }

    // buffer data bytes
    //int rd_ptr = 0;
    // from (0x301 + bytes_MHR) to (0x301 + frame_length - bytes_nodata - 1)
//    for (int i = 0; i < Mrf24j_rx_datalength(); i++) {
//      rx_info.rx_data[rd_ptr++] = Mrf24j_read_long(0x301 + bytes_MHR + i);
//    }

//    rx_info.frame_length = frame_length;
    // same as datasheet 0x301 + (m + n + 2) <-- frame_length
//    rx_info.lqi = Mrf24j_read_long(0x301 + frame_length);
    // same as datasheet 0x301 + (m + n + 3) <-- frame_length + 1
//    rx_info.rssi = Mrf24j_read_long(0x301 + frame_length + 1);

//    Mrf24j_rx_enable();
//    interrupts();
  }
  else
  {
	  last_interrupt = Mrf24j_read_short(MRF_INTSTAT);
//	  Mrf24j_reset();
//	    Mrf24j_init();
//	    Mrf24j_set_channel(20);

  }
  if (last_interrupt & MRF_I_TXNIF) {
    flag_got_tx++;
    uint8_t tmp = Mrf24j_read_short(MRF_TXSTAT);
    // 1 means it failed, we want 1 to mean it worked.
    tx_info.tx_ok = !(tmp & ~(1 << TXNSTAT));
    tx_info.retries = tmp >> 6;
    tx_info.channel_busy = (tmp & (1 << CCAFAIL));
  }
}

int Mrf24j_rx_datalength(void) {
  return rx_info.frame_length - bytes_nodata;
}
