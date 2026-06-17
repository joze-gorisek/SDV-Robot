/*
   Library: DDSM115.c

    Created on: Mar 7, 2025
       Author: Andy
    https://www.waveshare.com/wiki/DDSM115

 ******* Message structuchangere ****************************************************************************************************************************

  After sending a command, DDSM115 always sends the status of the motor.
    DATA[0]  DATA[1]      DATA[2]        DATA[3]          DATA[4]         DATA[5]           DATA[6]         DATA[7]        DATA[8]    DATA[9]
      [ID]     [MODE]   [Torque MSB]     [Torque LSB]  [Velocity MSB]   [Velocity LSB]   [Position MSB]   [Position LSB]    [ERROR]    [CRC8]

*******************************************************************************************************************************************************/

#include "main.h"
#include "DDSM115.h"


extern UART_HandleTypeDef huart5; //RS485
extern UART_HandleTypeDef huart2; //VCP




float current = 1.0f; //Desired values
float RPM     = 1.0f;
float angle   = 1.0f;

uint8_t motor_status_flag=0;

int toggle = 0;
char HEX_Buffer[15];
uint16_t angle_u;
uint8_t MotorID[MOTOR_NUMBER]={0x77}; //SET Motors ID !!
float   MotorData[MOTOR_NUMBER][4];   // Current, Velocity, Temp, Position
uint8_t ERRORDDSM[MOTOR_NUMBER][5];   // Troubleshooting, Stall error, Phase overcurrent error, Overcurrent error, Sensor error


uint8_t RS485_RxBuffer[RS485_BUFFER_SIZE];
uint8_t modeCmd[10] = {0};


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
	HAL_Delay(10);

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

}

void  VelocityMode(uint8_t motorID)
{
	mode[0]=motorID;
	mode[9]=0x02;
	HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_SET);
    HAL_UART_Transmit(&huart5, mode, 10 , 1000);
    HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);
}

void  PositionMode(uint8_t motorID)
{
	mode[0]=motorID;
	mode[9]=0x03;
	HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_SET);
	HAL_UART_Transmit(&huart5, mode, 10 , 1000);
	HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);

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

		// Fill command packet according to documentation:
		// Byte 0: Motor ID, Byte 1: Command code (0x64 for drive command)
		command[0] = motorID;
		command[1] = 0x64;

		// Bytes 2-3: 16-bit target position (big-endian)
		command[2] = (uint8_t)(target_value >> 8);   // High byte
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
		HAL_UART_Transmit(&huart5, command, 10, 5);
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

		// Fill command packet according to documentation:
		// Byte 0: Motor ID, Byte 1: Command code (0x64 for drive command)
		command[0] = motorID;
		command[1] = 0x64;

		// Bytes 2-3: 16-bit target position (big-endian)
		command[2] = (uint8_t)(target_value >> 8);   // High byte
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

	void BrakeVelocity(uint8_t motorID) {  //Works only in veloctiy mode
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
		command[7] = 0xFF; //Brake command
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

		// Fill command packet according to documentation:
		// Byte 0: Motor ID, Byte 1: Command code (0x64 for drive command)
		command[0] = motorID;
		command[1] = 0x64;

		// Bytes 2-3: 16-bit target position (big-endian)
		command[2] = (uint8_t)(target_value >> 8);   // High byte
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


/***********  Motor Status    ***********************************************************/


//	     DATA[0]  DATA[1]      DATA[2]        DATA[3]          DATA[4]         DATA[5]           DATA[6]         DATA[7]        DATA[8]    DATA[9]
//	      [ID]     [MODE]   [Torque MSB]     [Torque LSB]  [Velocity MSB]   [Velocity LSB]    [Temperature]   [Position U8]     [ERROR]    [CRC8]
	void MotorStatus(uint8_t motorID,uint8_t motorCnt) {

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
			motor_status_flag=1;
			HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_SET);
			// Transmit the command
			HAL_UART_Transmit(&huart5, command, 10, HAL_MAX_DELAY);
			// Return RS485 transceiver to receive mode
			HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);

			while(motor_status_flag==1);

			//Parse data from UART4 interrupt in RS485_RxBuffer[]
			MotorData [motorCnt][0] = (float)(256*RS485_RxBuffer[2]+RS485_RxBuffer[3]);   //Current
			MotorData [motorCnt][1] = (float)(256*RS485_RxBuffer[4]+RS485_RxBuffer[5]);   //Velocity
			MotorData [motorCnt][2] = (float)(RS485_RxBuffer[6]);                         //Temperature;
			MotorData [motorCnt][3] = (float)(RS485_RxBuffer[7])*360.0f/255.0f;           //Position only 8bit data

           //ERRORs
			ERRORDDSM [motorCnt][0]= RS485_RxBuffer[8] & 0x10; //Troubleshooting
			ERRORDDSM [motorCnt][1]= RS485_RxBuffer[8] & 0x08; //Stall error
			ERRORDDSM [motorCnt][2]= RS485_RxBuffer[8] & 0x04; //Phase overcurrent error
			ERRORDDSM [motorCnt][3]= RS485_RxBuffer[8] & 0x02; //Overcurrent error
			ERRORDDSM [motorCnt][4]= RS485_RxBuffer[8] & 0x01; //Sensor error


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



/********************************************************************/
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

