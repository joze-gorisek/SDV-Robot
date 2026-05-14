/*
   Library: DDSM115.c

    Created on: Mar 7, 2025
       Author: Andy
    https://www.waveshare.com/wiki/DDSM115

 ******* Message structure ****************************************************************************************************************************

  After sending a command, DDSM115 always sends the status of the motor.
    DATA[0]  DATA[1]      DATA[2]        DATA[3]          DATA[4]         DATA[5]           DATA[6]         DATA[7]        DATA[8]    DATA[9]
      [ID]     [MODE]   [Torque MSB]     [Torque LSB]  [Velocity MSB]   [Velocity LSB]   [Position MSB]   [Position LSB]    [ERROR]    [CRC8]

*******************************************************************************************************************************************************/
/* SET number of motors in DDSM115.h */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "main.h"
#include "DDSM115.h"


extern UART_HandleTypeDef huart5; //RS485
extern UART_HandleTypeDef huart2; //VCP


float current = 0.0f; //Desired values
float RPM     = 0.0f;
float angle   = 0.0f;

uint8_t motor_status_flag=0;

int toggle = 0;
uint8_t HEX_Buffer[15];
uint16_t angle_u;
uint8_t MotorID[MOTOR_NUMBER]={0x87}; //SET Motors ID !! SET numbers of motors in DDSM115.h!
float   MotorData[MOTOR_NUMBER][4];   // Current [-8 8]A [－32767 32767], Velocity [-300 300], Temperature, Position [0 360]deg [0 32767]
uint8_t ERRORDDSM[MOTOR_NUMBER][5];   // Troubleshooting, Stall error, Phase overcurrent error, Overcurrent error, Sensor error


uint8_t RS485_RxBuffer[RS485_BUFFER_SIZE];
uint8_t AllocBuffer[RS485_BUFFER_SIZE];
uint8_t modeCmd[10] = {0};


uint8_t MotorMode=3;      //FOR real-time motor parameters 1-current 2-velocity 3-position
uint16_t LastCommand=0;

uint8_t command[10] = {
    0xC8,  // Reserved
    0x64,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0xDE,   // CRC8 checksum
};


//  Current loop command（－32767～32767 corresponds to －8A～8A）
//  Velocity loop commands（－330～330 rpm）
//  Position loop commands（0～32767 corresponds to 0～360°）
uint8_t mode[10] = {
    0x01,  // Motor ID
    0xA0,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x03,  // CRC8 checksum 0x01-current 0x02-velocity[RPM] 0x03-position[deg]
};


uint8_t  ID_query[10] = {
    0xC8,  // Reserved
    0x64,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0x00,  // Reserved
    0xDE,   // CRC8 checksum
};

/******  ID query function **********************************************************/

uint8_t GetMotorID() //Only one motor on the bus!!!
{

	uint8_t ID;

	HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_SET);
	HAL_UART_Transmit(&huart5, ID_query, 10 , 100);
	HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);
	/*Feedback in IRQ_UART4 -> RS485_RxBuffer[10]*/

	/*Wait*/
	HAL_Delay(100);

	/*VCP message*/
	//HAL_UART_Transmit(&huart2, HEX_Buffer, strlen(HEX_Buffer),1000);


   return  RS485_RxBuffer[0];
}

/**** Motor modes *********************************************************************/
void CurrentMode(uint8_t motorID)
{
	  mode[0]=motorID;
	  mode[9]=0x01;
	  HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_SET);
	  HAL_UART_Transmit(&huart5, mode, 10 , 1000);
	  HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);
	  MotorMode=1;

}

void  VelocityMode(uint8_t motorID)
{
	mode[0]=motorID;
	mode[9]=0x02;
	HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_SET);
    HAL_UART_Transmit(&huart5, mode, 10 , 1000);
    HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);
    MotorMode=2;
}

void  PositionMode(uint8_t motorID)
{
	mode[0]=motorID;
	mode[9]=0x03;
	HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_SET);
	HAL_UART_Transmit(&huart5, mode, 10 , 1000);
	HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);
	MotorMode=3;

}

/******   Current functions ***************************************************************/

	static int16_t CurrentToValue(float current) {
		// Clamp the angle between -8A and 8A
		if(current < -8.0f) {
			current = -8.0f;
		}
		if(current > 8.0f) {
			current = 8.0f;
		}
		// Map [-8A <-> 8A] to [－32767 <-> 32767]
		return (int16_t)((current / 8.0f) * 32767.0f);
	}

	void sendCurrentCommand(uint8_t motorID, float current) {
		uint8_t command[10] = {0};
		uint16_t target_value =CurrentToValue(current);
		LastCommand = target_value; //For real-time data acquisition

		// Fill command packet according to documentation:
		// Byte 0: Motor ID, Byte 1: Command code (0x64 for drive command)
		command[0] = motorID;
		command[1] = 0x64;

		// Bytes 2-3: 16-bit target position (big-endian)
		command[2] = (uint8_t)(target_value >> 8) & 0xFF;   // High byte
		command[3] = (uint8_t)(target_value & 0xFF);   // Low byte

		// Bytes 4-8: Reserved (set to 0)
		command[4] = 0x00;
		command[5] = 0x00;
		command[6] = 0x00;
		command[7] = 0x00;
		command[8] = 0x00;

		// Byte 9: CRC8 checksum over bytes 0 to 8
		command[9] = compute_crc8(command, 9);

		// Set RS485 transceiver to transmit mode
		HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_SET);
		// Transmit the command
		HAL_UART_Transmit(&huart5, command, 10, HAL_MAX_DELAY);
		// Return RS485 transceiver to receive mode
		HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);
		// After sending a command, DDSM115 always sends the status of the motor.

	}



/******   Velocity functions *****************************************************************/

	 static int16_t VelocityToValue(float velocity) {
		// Clamp the angle between -300RPM and 300RPM
		if(velocity < -300.0f) {
			velocity = -300.0f;
		}
		if(velocity > 300.0f) {
			velocity = 300.0f;
		}
		// Map [-300RPM <-> 300RPM]
		return (int16_t)velocity ;
	}

	void sendVelocityCommand(uint8_t motorID, float velocity) {
		uint8_t command[10] = {0};
		uint16_t target_value = VelocityToValue(velocity);
		LastCommand = target_value; //For real-time data acquisition

		// Fill command packet according to documentation:
		// Byte 0: Motor ID, Byte 1: Command code (0x64 for drive command)
		command[0] = motorID;
		command[1] = 0x64;

		// Bytes 2-3: 16-bit target position (big-endian)
		command[2] = (uint8_t)(target_value >> 8) & 0xFF;   // High byte
		command[3] = (uint8_t)(target_value & 0xFF);   // Low byte

		// Bytes 4-8: Reserved (set to 0)
		command[4] = 0x00;
		command[5] = 0x00;
		command[6] = 0x00;  //MAX acceleration time 1RPM - 0.1ms
		command[7] = 0x00;
		command[8] = 0x00;

		// Byte 9: CRC8 checksum over bytes 0 to 8
		command[9] = compute_crc8(command, 9);

		// Set RS485 transceiver to transmit mode
		HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_SET);
		// Transmit the command
		HAL_UART_Transmit(&huart5, command, 10, HAL_MAX_DELAY);
		// Return RS485 transceiver to receive mode
		HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);
		// After sending a command, DDSM115 always sends the status of the motor.


		                //Parse motor data

		            //  MotorData [0][0] = (float) ( (int16_t) ( (RS485_RxBuffer[2]<<8) | RS485_RxBuffer[3] ) ) * 8.0f/32767;   //Current
					 //	MotorData [0][1] = (float) ( (int16_t) ( (RS485_RxBuffer[4]<<8) | RS485_RxBuffer[5] ) );   //Velocity
					//	MotorData [0][2] = 22.0f;                                              //Temperature;
					//	MotorData [0][3] = ( (float)  ((RS485_RxBuffer[6]<<8) | RS485_RxBuffer[7])) * 360.0f/32767.0f;   //Position



	}

	void MotorStop(uint8_t motorID) {  // Works only in velocity mode
		uint8_t command[10] = {0};


		// Fill command packet according to documentation:
		// Byte 0: Motor ID, Byte 1: Command code (0x64 for drive command)
		command[0] = motorID;
		command[1] = 0x64;

		// Bytes 2-3: 16-bit target position (big-endian)
		command[2] = 0x00;   // High byte
		command[3] = 0x00;   // Low byte

		// Bytes 4-8: Reserved (set to 0)
		command[4] = 0x00;
		command[5] = 0x00;
		command[6] = 0x00;
		command[7] = 0xFF; //Brake
		command[8] = 0x00;

		// Byte 9: CRC8 checksum over bytes 0 to 8
		command[9] = compute_crc8(command, 9);

		// Set RS485 transceiver to transmit mode
		HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_SET);
		// Transmit the command
		HAL_UART_Transmit(&huart5, command, 10, HAL_MAX_DELAY);
		// Return RS485 transceiver to receive mode
		HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);
		// After sending a command, DDSM115 always sends the status of the motor.

	}


/******   Position functions *************************************************************************/

    // Function to convert an angle in degrees to a 16-bit value (0 to 32767) ENCODER resolution  4096
	static uint16_t AngleToValue(float angle_deg) {
		// Clamp the angle between 0 and 360
		if(angle_deg < 0.0f) {
			angle_deg = 0.0f;
		}
		if(angle_deg > 360.0f) {
			angle_deg = 360.0f;
		}
		// Map 0-360° to 0-32767 (note: 32767 is the maximum unsigned 16-bit value used)
		return (uint16_t)((angle_deg / 360.0f) * 32767.0f);
	}

	void sendPositionCommand(uint8_t motorID, float angle_deg) {
		uint8_t command[10] = {0};
		uint16_t target_value = AngleToValue(angle_deg);
		LastCommand = target_value; //For real-time data acquisition


		// Fill command packet according to documentation:
		// Byte 0: Motor ID, Byte 1: Command code (0x64 for drive command)
		command[0] = motorID;
		command[1] = 0x64;

		// Bytes 2-3: 16-bit target position (big-endian)
		command[2] = (uint8_t)(target_value >> 8) & 0xFF;   // High byte
		command[3] = (uint8_t)(target_value & 0xFF);   // Low byte

		// Bytes 4-8: Reserved (set to 0)
		command[4] = 0x00;
		command[5] = 0x00;
		command[6] = 0x00;
		command[7] = 0x00;
		command[8] = 0x00;

		// Byte 9: CRC8 checksum over bytes 0 to 8
		command[9] = compute_crc8(command, 9);

		// Set RS485 transceiver to transmit mode
		HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_SET);
		// Transmit the command
		HAL_UART_Transmit(&huart5, command, 10, HAL_MAX_DELAY);
		// Return RS485 transceiver to receive mode
		HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);
		// After sending a command, DDSM115 always sends the status of the motor.

	}

/******   Switch mode *************************************************************************/
/*
 * mode = 1 current
 * mode = 2 velocity
 * mode = 3 position
 * **********************************************************************************************/
	void MotorSwitchMode(uint8_t motorID,uint8_t mode) {
             /*motorCnt start with 1 for first motor*/
		   uint8_t command[10] = {0};

             if(mode<1 || mode > 3)
            	 mode=3;  //Default mode


			command[0] = motorID;
			command[1] = 0xA0; // Status

			// Bytes 2-8: Reserved (set to 0)
			command[2] = 0x00;
			command[3] = 0x00;
			command[4] = 0x00;
			command[5] = 0x00;
			command[6] = 0x00;
			command[7] = 0x00;
			command[8] = 0x00;


			command[9] = mode;

			// Set RS485 transceiver to transmit mode
			HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_SET);
			// Transmit the command
			HAL_UART_Transmit(&huart5, command, 10, HAL_MAX_DELAY);
			// Return RS485 transceiver to receive mode
			HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);


		}


/***********  Motor Status    ***********************************************************/

	//	     DATA[0]  DATA[1]      DATA[2]        DATA[3]          DATA[4]         DATA[5]           DATA[6]         DATA[7]        DATA[8]    DATA[9]
	//	      [ID]     [MODE]   [Torque MSB]     [Torque LSB]  [Velocity MSB]   [Velocity LSB]    [Temperature]   [Position U8]     [ERROR]    [CRC8]
	void MotorStatus(uint8_t motorID, uint8_t motorCnt) { //NOT WROKING

			    uint8_t command[10] = {0};

			    		// Fill command packet according to documentation:
			    		// Byte 0: Motor ID, Byte 1: Command code (0x64 for drive command)
			    		command[0] = motorID;
			    		command[1] = 0x64;

			    		// Bytes 2-3: 16-bit target position (big-endian)
			    		command[2] = (uint8_t)(LastCommand >> 8) & 0xFF;   // High byte
			    		command[3] = (uint8_t)(LastCommand & 0xFF);   // Low byte

			    		// Bytes 4-8: Reserved (set to 0)
			    		command[4] = 0x00;
			    		command[5] = 0x00;
			    		command[6] = 0x00;
			    		command[7] = 0x00;
			    		command[8] = 0x00;

						// Byte 9: CRC8 checksum over bytes 0 to 8
						command[9] = compute_crc8(command, 9);

				// Set RS485 transceiver to transmit mode
				//motor_status_flag=1;
				HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_SET);
				// Transmit the command
				HAL_UART_Transmit(&huart5, command, 10, HAL_MAX_DELAY);
				// Return RS485 transceiver to receive mode
				HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);


				//Parse data from UART4 interrupt in RS485_RxBuffer[]

		         MotorData [0][0] = (float) ( (int16_t) ( (RS485_RxBuffer[2]<<8) | RS485_RxBuffer[3] ) ) * 8.0f/32767;   //Current
				 MotorData [0][1] = (float) ( (int16_t) ( (RS485_RxBuffer[4]<<8) | RS485_RxBuffer[5] ) );                //Velocity
			     MotorData [0][2] = 0;                                                //Temperature;
			     MotorData [0][3] = ( (float)  ((RS485_RxBuffer[4]<<8) |RS485_RxBuffer[7])) * 360.0f/255.0f;   //Position


	           //ERRORs
				ERRORDDSM [motorCnt-1][0]= RS485_RxBuffer[8] & 0x10; //Troubleshooting
				ERRORDDSM [motorCnt-1][1]= RS485_RxBuffer[8] & 0x08; //Stall error
				ERRORDDSM [motorCnt-1][2]= RS485_RxBuffer[8] & 0x04; //Phase over current error
				ERRORDDSM [motorCnt-1][3]= RS485_RxBuffer[8] & 0x02; //Over current error
				ERRORDDSM [motorCnt-1][4]= RS485_RxBuffer[8] & 0x01; //Sensor error


			}

//	     DATA[0]  DATA[1]      DATA[2]        DATA[3]          DATA[4]         DATA[5]           DATA[6]         DATA[7]        DATA[8]    DATA[9]
//	      [ID]     [MODE]   [Torque MSB]     [Torque LSB]  [Velocity MSB]   [Velocity LSB]    [Temperature]   [Position U8]     [ERROR]    [CRC8]
	void MotorStatusOLD(uint8_t motorID, uint8_t motorCnt) { //NOT WROKING!!!!
             /*motorCnt start with 1 for first motor*/
		    uint8_t command[10] = {0};

			// Fill command packet according to documentation:
			// Byte 0: Motor ID, Byte 1: Command code (0x64 for drive command)

			command[0] = motorID;
			command[1] = 0x74; // Status

			// Bytes 2-8: Reserved (set to 0)
			command[2] = 0x00;
			command[3] = 0x00;
			command[4] = 0x00;
			command[5] = 0x00;
			command[6] = 0x00;
			command[7] = 0x00;
			command[8] = 0x00;

			// Byte 9: CRC8 checksum over bytes 0 to 8
			command[9] = compute_crc8(command, 9);

			// Set RS485 transceiver to transmit mode
			//motor_status_flag=1;
			HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_SET);
			// Transmit the command
			HAL_UART_Transmit(&huart5, command, 10, HAL_MAX_DELAY);
			// Return RS485 transceiver to receive mode
			HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);



			//while(motor_status_flag==1);
            //HAL_Delay(10);
			//Parse data from UART4 interrupt in RS485_RxBuffer[]
	         MotorData [0][0] = (float) ( (int16_t) ( (RS485_RxBuffer[2]<<8) | RS485_RxBuffer[3] ) ) * 8.0f/32767;   //Current
		     MotorData [0][1] = (float) ( (int16_t) ( (RS485_RxBuffer[4]<<8) | RS485_RxBuffer[5] ) );                //Velocity
			 MotorData [0][2] = (float) ((int16_t) RS485_RxBuffer[6]);                                                //Temperature;
			 MotorData [0][3] = ( (float)  (RS485_RxBuffer[7])) * 360.0f/255.0f;   //Position


           //ERRORs
			ERRORDDSM [motorCnt-1][0]= RS485_RxBuffer[8] & 0x10; //Troubleshooting
			ERRORDDSM [motorCnt-1][1]= RS485_RxBuffer[8] & 0x08; //Stall error
			ERRORDDSM [motorCnt-1][2]= RS485_RxBuffer[8] & 0x04; //Phase over current error
			ERRORDDSM [motorCnt-1][3]= RS485_RxBuffer[8] & 0x02; //Over current error
			ERRORDDSM [motorCnt-1][4]= RS485_RxBuffer[8] & 0x01; //Sensor error
			char RX_Buff[40];

			sprintf(RX_Buff,"MOT%d  %.1f  %.1f  %.1f    %.1f\n\r", motorID,
					                                               MotorData [motorCnt-1][3],
																   MotorData [motorCnt-1][1],
																   MotorData [motorCnt-1][0],
																   MotorData [motorCnt-1][2]);


		     HAL_UART_Transmit(&huart2, (uint8_t *)RX_Buff, strlen(RX_Buff),1000);


		}

/******  Change Motor ID  **********************************************************/

//USE  WHEN JUST ONE MOTOR IS ON NETWORK!!!
	void ChangeMotorID(uint8_t NewMotorID)
	{

	        if(NewMotorID>255) NewMotorID=255;
	        if(NewMotorID<0)   NewMotorID=0x02;

			command[0] = 0XAA;
			command[1] = 0x55; // Status
			command[2] = 0x53;
			command[3] = NewMotorID;
			command[4] = 0x00;
			command[5] = 0x00;
			command[6] = 0x00;
			command[7] = 0x00;
			command[8] = 0x00;
			command[9] = 0x00;


			for(int i=0;i<5;i++)
			{
				HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_SET);
				// Transmit the command
				HAL_UART_Transmit(&huart5, command, 10, HAL_MAX_DELAY);
				// Return RS485 transceiver to receive mode
				HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);
				HAL_Delay(10);
			}


	}


/*********** CRC8 *************************************************************************/

		// CRC8 calculation using CRC-8/MAXIM (often uses a reversed polynomial of 0x8C)
		uint8_t compute_crc8(uint8_t *data, uint8_t len) {
			uint8_t crc = 0;
			for (uint8_t i = 0; i < len; i++) {
				crc ^= data[i];
				for (uint8_t j = 0; j < 8; j++) {
					if (crc & 0x01)
						crc = (crc >> 1) ^ 0x8C;
					else
						crc >>= 1;
				}
			}
			return crc;
		}

/********************************************************************/



void ParseData(uint8_t *Buffer,uint8_t motorNUM)
{
	   //motorNUM - motor number

        MotorData [motorNUM-1][0] = (float) ( (int16_t) ( (Buffer[2]<<8) | Buffer[3] ) ) * 8.0f/32767;   //Current
	    MotorData [motorNUM-1][1] = (float) ( (int16_t) ( (Buffer[4]<<8) | Buffer[5] ) );                //Velocity
	    MotorData [motorNUM-1][2] = 22.0f;  //NO TEMPERATURE MEASUREMENT in mode 0x64
	    MotorData [motorNUM-1][3] = ( (float)  ((Buffer[6]<<8) | Buffer[7])) * 360.0f/32767.0f;         //Position

}


/* EXAMPLE: RS485 interrupt in  UART4_IRQHandler(void) */

/*

 HAL_UART_IRQHandler(&huart4);

 if((UART4->SR & 0x0525U)!= RESET)
	{

	   RS485_RxBuffer[cnt]=UART4->DR;

	   cnt++;

	   if(cnt==10)
		   {
		      cnt=0;  //Rest counter
		      for(int i=0; i<RS485_BUFFER_SIZE; i++)
		      {
		    	  sprintf(HEX_Buffer,"0x%x ",RS485_RxBuffer[i]);
		    	  HAL_UART_Transmit(&huart2, HEX_Buffer, strlen(HEX_Buffer),1000);

		      }
		   }
	}
*/

