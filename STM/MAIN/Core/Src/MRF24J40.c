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

static SPI_HandleTypeDef* mrfspi;
static UART_HandleTypeDef* duart;
static uint8_t uart_debug=1; // 1 - on, 0 - off;

static GPIO_TypeDef* MRF_CS_Port;
static GPIO_TypeDef* MRF_CE_Port;
static GPIO_TypeDef* MRF_INT_Port;
static uint16_t MRF_CS;
static uint16_t MRF_CE;
static uint16_t MRF_INT;

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
volatile uint8_t mrf_rx_data_len;

rx_info_t rx_info;
tx_info_t tx_info;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

uint8_t last_seq = 0xFF;

static void Mrf24j_print_rx_fifo(void)
{
    uint8_t len = Mrf24j_read_long(0x300);
    uint8_t buf[128];

    if ((len > bytes_nodata) && (len < sizeof(buf))) {
        for (int i = 0; i < len; i++) {
            buf[i] = Mrf24j_read_long(0x301 + i);
        }

        Mrf24j_rx_flush();

        uint8_t *data = &buf[bytes_MHR];
        uint8_t data_len = len - bytes_nodata;

        HAL_UART_Transmit(duart, data, data_len, HAL_MAX_DELAY);
        HAL_UART_Transmit(duart, (uint8_t*)"\r\n", 2, HAL_MAX_DELAY);
    } else {
        Mrf24j_rx_flush();
    }
}

void Mrf24j_poll_uart(void)
{
	uint8_t intstat = Mrf24j_read_short(MRF_INTSTAT);

	if (intstat & MRF_I_RXIF) {
		HAL_UART_Transmit(duart, (uint8_t*)"MRF POLL RX\r\n", 13, 100);
		Mrf24j_print_rx_fifo();
	}
}

void Mrf24j_rx_check_uart(void)
{
	uint8_t intstat = Mrf24j_read_short(MRF_INTSTAT);
	uint8_t intcon = Mrf24j_read_short(MRF_INTCON);
	uint8_t rxmcr = Mrf24j_read_short(MRF_RXMCR);
	uint8_t bbreg1 = Mrf24j_read_short(MRF_BBREG1);
	char msg[80];
	int msg_len = snprintf(msg, sizeof(msg),
						   "RX CHECK INTSTAT=0x%02X INTCON=0x%02X RXMCR=0x%02X BBREG1=0x%02X\r\n",
						   intstat, intcon, rxmcr, bbreg1);
	if (msg_len > 0) {
		HAL_UART_Transmit(duart, (uint8_t*)msg, (uint16_t)msg_len, 100);
	}

	if (intstat & MRF_I_RXIF) {
		Mrf24j_print_rx_fifo();
	}
}

void MRF_InitSPI(SPI_HandleTypeDef *hspi, UART_HandleTypeDef *huart, uint8_t debug){
	mrfspi = hspi;
	duart = huart;
	uart_debug = debug;
}

void MRF_InitGPIO(GPIO_TypeDef* CS_Port, uint16_t CS_Pin, GPIO_TypeDef* CE_Port, uint16_t CE_Pin, GPIO_TypeDef* IRQ_Port, uint16_t IRQ_Pin){
	MRF_CS_Port = CS_Port;
	MRF_CE_Port = CE_Port;
	MRF_INT_Port = IRQ_Port;
	MRF_CS = CS_Pin;
	MRF_CE = CE_Pin;
	MRF_INT = IRQ_Pin;
	HAL_GPIO_WritePin(MRF_CS_Port, MRF_CS, GPIO_PIN_SET);
}

static void Mrf24j_spi_flush(void) {
    while(mrfspi->Instance->SR & SPI_SR_RXNE) {
        volatile uint8_t dummy = mrfspi->Instance->DR;
    }
}

static void Mrf24j_debug_regs(const char *label)
{
	if ((duart != NULL) && (uart_debug == 1))
	{
		char msg[96];
		uint16_t pan = Mrf24j_get_pan();
		uint8_t rfcon0 = Mrf24j_read_long(MRF_RFCON0);
		uint8_t rfcon2 = Mrf24j_read_long(MRF_RFCON2);
		uint8_t rfcon6 = Mrf24j_read_long(MRF_RFCON6);
		uint8_t rfstate = Mrf24j_read_long(RFSTATE);
		int msg_len = snprintf(msg, sizeof(msg),
							   "%s PAN=0x%04X RFCON0=0x%02X RFCON2=0x%02X RFCON6=0x%02X RFSTATE=0x%02X\r\n",
							   label, pan, rfcon0, rfcon2, rfcon6, rfstate);
		if (msg_len > 0)
		{
			HAL_UART_Transmit(duart, (uint8_t*)msg, (uint16_t)msg_len, 100);
		}
	}
}

static void Mrf24j_debug_write_long(const char *label, uint16_t address, uint8_t data)
{
	Mrf24j_write_long(address, data);
	HAL_Delay(2);

	if ((duart != NULL) && (uart_debug == 1))
	{
		char msg[64];
		uint8_t readback = Mrf24j_read_long(address);
		int msg_len = snprintf(msg, sizeof(msg),
							   "%s addr=0x%03X wrote=0x%02X read=0x%02X\r\n",
							   label, address, data, readback);
		if (msg_len > 0)
		{
			HAL_UART_Transmit(duart, (uint8_t*)msg, (uint16_t)msg_len, 100);
		}
	}
}


static void Mrf24j_init_common(uint8_t ch, bool do_reset)
{

	if (do_reset)
	{
		Mrf24j_reset();
	}
	Mrf24j_debug_regs("MRF DBG after reset");
	Mrf24j_debug_write_long("MRF WR RFCON0 after reset", MRF_RFCON0, 0x93);
	Mrf24j_debug_regs("MRF DBG after reset RFCON0 probe");
	uint8_t check;

	/* BASE SETUP */
	// from mrf24j40 datasheet p90 3.2 init
	if (do_reset)
	{
		Mrf24j_write_short(MRF_SOFTRST, 0x7);

		uint8_t softrst;
		uint32_t soft_reset_start = HAL_GetTick();
		do {
			softrst = Mrf24j_read_short(MRF_SOFTRST);
			if ((HAL_GetTick() - soft_reset_start) > 1000U)
			{
				if (duart != NULL)
				{
					char msg[48];
					int msg_len = snprintf(msg, sizeof(msg),
										   "MRF SOFTRST timeout=0x%02X\r\n",
										   softrst);
					if (msg_len > 0)
					{
						HAL_UART_Transmit(duart, (uint8_t*)msg, (uint16_t)msg_len, 100);
					}
				}
				return;
			}
		} while ((softrst & 0x7) != 0);
		HAL_Delay(2);
		Mrf24j_debug_regs("MRF DBG after soft reset");
	}
	else
	{
		if ((duart != NULL) && (uart_debug == 1))
		{
			HAL_UART_Transmit(duart, (uint8_t*)"MRF DBG skip soft reset\r\n", 25, 100);
		}
	}

	/* RF SETUP */
	Mrf24j_write_short(MRF_PACON2, 0x98); // – Initialize FIFOEN = 1 and TXONTS = 0x6.
	Mrf24j_write_short(MRF_TXSTBL, 0x95); // – Initialize RFSTBL = 0x9.

	Mrf24j_debug_write_long("MRF WR RFCON0 base", MRF_RFCON0, 0x03); // – Initialize RFOPT = 0x03.
	Mrf24j_debug_write_long("MRF PROBE RFCON6 before RFCON0", MRF_RFCON6, 0x90);
	HAL_Delay(1);
	uint8_t val = Mrf24j_read_long(MRF_RFCON0);
	Mrf24j_debug_write_long("MRF WR RFCON1", MRF_RFCON1, 0x02); // – Initialize VCOOPT = 0x02.
	Mrf24j_debug_write_long("MRF WR RFCON2", MRF_RFCON2, 0x80); // – Enable PLL (PLLEN = 1).
	Mrf24j_debug_write_long("MRF WR RFCON3", MRF_RFCON3, 0x00); // – set tx max power
	Mrf24j_debug_write_long("MRF WR RFCON6", MRF_RFCON6, 0x90); // – Initialize TXFIL = 1 and 20MRECVR = 1.
	HAL_Delay(1);
	val = Mrf24j_read_long(MRF_RFCON6);
	Mrf24j_debug_write_long("MRF WR RFCON7", MRF_RFCON7, 0x80); // – Initialize SLPCLKSEL = 0x2 (100 kHz Internal oscillator).
	Mrf24j_debug_write_long("MRF WR RFCON8", MRF_RFCON8, 0x10); // – Initialize RFVCO = 1.
	Mrf24j_debug_regs("MRF DBG after RF writes");
//	Mrf24j_read_long(MRF_SLPCON0);
//	Mrf24j_write_long(MRF_SLPCON0, 0x02); //interrupt polarity - rising edge
//	Mrf24j_write_long(MRF_SLPCON1, 0x21); // – Initialize CLKOUTEN = 1 and SLPCLKDIV = 0x01.
	HAL_Delay(10);

	/* BASEBAND SETUP*/
	//  Configuration for nonbeacon-enabled devices (see Section 3.8 “Beacon-Enabled and
	//  Nonbeacon-Enabled Networks”):
	Mrf24j_write_short(MRF_BBREG2, 0x80); // Set CCA mode to ED
	Mrf24j_write_short(MRF_BBREG6, 0x40); // – Set appended RSSI value to RXFIFO.
	Mrf24j_write_short(MRF_CCAEDTH, 0x60); // – Set CCA ED threshold.

	/* NETWORK SETUP */
	Mrf24j_address16_write(0x0011);
	Mrf24j_set_pan(0x2222);
	Mrf24j_set_channel(ch);
	HAL_Delay(2);
	if ((duart != NULL) && (uart_debug == 1))
	{
		char msg[64];
		uint8_t readback = Mrf24j_read_long(MRF_RFCON0);
		int msg_len = snprintf(msg, sizeof(msg),
							   "MRF WR channel ch=%u RFCON0 read=0x%02X\r\n",
							   ch, readback);
		if (msg_len > 0)
		{
			HAL_UART_Transmit(duart, (uint8_t*)msg, (uint16_t)msg_len, 100);
		}
	}
	Mrf24j_debug_regs("MRF DBG after net setup");

	/* MRF RESET */
	// max power is by default.. just leave it...
	// Set transmitter power - See “REGISTER 2-62: RF CONTROL 3 REGISTER (ADDRESS: 0x203)”.
	Mrf24j_write_short(MRF_RFCTL, 0x04); //  – Reset RF state machine.
	Mrf24j_write_short(MRF_RFCTL, 0x00); // part 2
	HAL_Delay(10);
	Mrf24j_debug_regs("MRF DBG after RFCTL");

	/* WAIT FOR RESET */
	Mrf24j_debug_write_long("MRF WR RFCON2 wait", MRF_RFCON2, 0x80); // – Enable PLL (PLLEN = 1).
	HAL_Delay(500);
	val = Mrf24j_read_long(MRF_RFCON2);
	Mrf24j_debug_regs("MRF DBG before RF wait");
	uint32_t rf_wait_start = HAL_GetTick();
	do {
		val = Mrf24j_read_long(MRF_RFCON6); // should return 0x03
		check = Mrf24j_read_long(RFSTATE);
		if ((HAL_GetTick() - rf_wait_start) > 1000U)
		{
			if (duart != NULL)
			{
				char msg[64];
				int msg_len = snprintf(msg, sizeof(msg),
									   "MRF RF timeout RFCON6=0x%02X RFSTATE=0x%02X\r\n",
									   val, check);
				if (msg_len > 0)
				{
					HAL_UART_Transmit(duart, (uint8_t*)msg, (uint16_t)msg_len, 100);
				}
			}
			return;
		}
	} while ((check & 0xa0) != 0xa0);

	/* ENABLE RECEIVER */
	Mrf24j_write_short(MRF_RXMCR, 0x01);
	Mrf24j_rx_flush();

	/* ENABLE RX INTERRUPTS ONLY */
	Mrf24j_set_interrupts();

	HAL_Delay(2);
}

void Mrf24j_init(uint8_t ch)
{
	Mrf24j_init_common(ch, true);
}

void Mrf24j_init_no_reset(uint8_t ch)
{
	Mrf24j_init_common(ch, false);
}


void Mrf24j_write_short(uint8_t address, uint8_t data)
{
	  HAL_StatusTypeDef test = 0;
	  HAL_GPIO_WritePin(MRF_CS_Port, MRF_CS, GPIO_PIN_RESET);

	  address = (address<<1 & 0b01111110) | 0x01;
	  test = HAL_SPI_Transmit(mrfspi, &address, 1,100);
	  if(test != HAL_OK && uart_debug == 1)
	  {
		  HAL_UART_Transmit (duart,(uint8_t*) "napaka1\r\n", strlen("napaka1\r\n"), HAL_MAX_DELAY);
	  }
	  else if( uart_debug == 1)
	  {
		  HAL_UART_Transmit (duart, (uint8_t*)"ok1\r\n", strlen("ok1\r\n"), HAL_MAX_DELAY);
	  }
	  test = HAL_SPI_Transmit(mrfspi, &data, 1,100);
	  if(test != HAL_OK && uart_debug == 1)
	  {
		  HAL_UART_Transmit (duart, (uint8_t*)"napaka2\r\n", strlen("napaka2\r\n"), HAL_MAX_DELAY);
		  HAL_UART_Transmit (duart, &test, 1, HAL_MAX_DELAY);
	  }
	  else if(uart_debug == 1)
	  {
		  HAL_UART_Transmit (duart, (uint8_t*)"ok2\r\n", strlen("ok2\r\n"), HAL_MAX_DELAY);
	  }

	  HAL_Delay(1);
	  HAL_GPIO_WritePin(MRF_CS_Port, MRF_CS, GPIO_PIN_SET);
	  HAL_Delay(1);

	}
void Mrf24j_write_long(uint16_t address, uint8_t data)
{
    HAL_GPIO_WritePin(MRF_CS_Port, MRF_CS, GPIO_PIN_RESET);

    uint8_t ahigh = 0x80 | ((address >> 3) & 0x7F);  // long address flag
    uint8_t alow  = ((address & 0x07) << 5) | 0x10;  // addr bits + WRITE bit (bit4=1)

    HAL_SPI_Transmit(mrfspi, &ahigh, 1, 10);
    HAL_SPI_Transmit(mrfspi, &alow,  1, 10);
    HAL_SPI_Transmit(mrfspi, &data,  1, 10);

    HAL_GPIO_WritePin(MRF_CS_Port, MRF_CS, GPIO_PIN_SET);
}


uint8_t Mrf24j_read_short(uint8_t address)
{

	Mrf24j_spi_flush();

	  HAL_GPIO_WritePin(MRF_CS_Port, MRF_CS, GPIO_PIN_RESET);

	  address = address<<1 & 0b01111110;
	  HAL_SPI_Transmit(mrfspi, &address, 1,100);

	  uint8_t ret = 0x0;
	  uint8_t bogus = 0x00;
	  HAL_SPI_TransmitReceive(mrfspi, &bogus, &ret, 1, 100);
	//  HAL_SPI_Receive(spi, &ret,1,1);

	  HAL_GPIO_WritePin(MRF_CS_Port, MRF_CS, GPIO_PIN_SET);

	  return ret;
	}

uint8_t Mrf24j_read_long(uint16_t address)
{
	Mrf24j_spi_flush();

	  HAL_GPIO_WritePin(MRF_CS_Port, MRF_CS, GPIO_PIN_RESET);

	  uint8_t ahigh = 0x80 | ((address >> 3) & 0x7F); // R/W=1, long address
	  uint8_t alow  = (uint8_t) ((address & 0x07) << 5);          // lower 3 bits in bits 7–5

	  HAL_SPI_Transmit(mrfspi, &ahigh, 1,100);
	  HAL_SPI_Transmit(mrfspi, &alow, 1,100);

	  uint8_t ret = 0x0;
	  uint8_t bogus = 0x00;
	  HAL_SPI_TransmitReceive(mrfspi, &bogus, &ret, 1, 100);
	  //HAL_SPI_Receive(spi, &ret,1,1);

	  HAL_GPIO_WritePin(MRF_CS_Port, MRF_CS, GPIO_PIN_SET);

	  return ret;
	}
void Mrf24j_set_interrupts(void)
{
	uint8_t intcon = 0xFF;

	// enable ONLY RX interrupt
	intcon &= ~(1 << 3); // RXIE = 0
	Mrf24j_write_short(MRF_INTCON, intcon);
	HAL_Delay(1);
	uint8_t val = Mrf24j_read_short(MRF_INTCON);
	HAL_Delay(1);
}
void Mrf24j_set_channel(uint8_t channel)
{
  //  (((channel - 11) << 4) | 0x03));
	Mrf24j_write_long(MRF_RFCON0, (((channel - 11) << 4) | 0x03));
}


void Mrf24j_reset(void)
{
  HAL_GPIO_WritePin(MRF_CE_Port, MRF_CE, GPIO_PIN_RESET);
  HAL_Delay(100);
  HAL_GPIO_WritePin(MRF_CE_Port, MRF_CE, GPIO_PIN_SET);
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

  uint16_t panid = Mrf24j_get_pan();

  Mrf24j_write_long(i++, panid & 0xff);  // dest panid
  Mrf24j_write_long(i++, panid >> 8);

  // Dest
  Mrf24j_write_long(i++, dest16 & 0xff);  // dest16 low
  Mrf24j_write_long(i++, dest16 >> 8); // dest16 high
  // Source
  uint16_t src16 = Mrf24j_address16_read();

  Mrf24j_write_long(i++, src16 & 0xff); // src16 low
  Mrf24j_write_long(i++, src16 >> 8); // src16 high

  // All testing seems to indicate that the next two bytes are ignored.
  //2 bytes on FCS appended by TXMAC
  i+=ignoreBytes;
  for (int q = 0; q < len; q++) {
	  Mrf24j_write_long(i++, data[q]);
  }
  Mrf24j_write_short(MRF_TXNCON, (1<<MRF_TXNTRIG));
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

  if (last_interrupt & MRF_I_RXIF) {
    uint8_t frame_length = Mrf24j_read_long(0x300);
    rx_info.frame_length = frame_length;

    if ((frame_length > bytes_nodata) && (frame_length < sizeof(rx_buf)))
    {
      int rb_ptr = 0;
      for (int i = 1; i < frame_length + 2; i++) {
        rx_buf[rb_ptr++] = Mrf24j_read_long(0x300 + i);
      }

      int data_len = frame_length - bytes_nodata;
      if (data_len > 0 && data_len <= (int)sizeof(rx_info.rx_data))
      {
        memcpy(rx_info.rx_data, &rx_buf[bytes_MHR], (size_t)data_len);
        mrf_rx_data_len = (uint8_t)data_len;
        flag_got_rx++;
        HAL_GPIO_TogglePin(LD1_GPIO_Port, LD1_Pin);
      }
      Mrf24j_rx_flush();
    }
    else
    {
    	Mrf24j_write_short(MRF_RXFLUSH, 0x05);
    }
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
