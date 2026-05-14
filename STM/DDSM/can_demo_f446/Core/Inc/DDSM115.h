/*
 * DDSM115.h
 *
 *  Created on: Mar 7, 2025
 *      Author: Andy
 *      https://www.waveshare.com/wiki/DDSM115
 */

#ifndef INC_DDSM115_H_
#define INC_DDSM115_H_

#define MOTOR_NUMBER 1   //Motors number in the network
#define RS485_BUFFER_SIZE 10
#define PACKET_SIZE 10



/*ID*/
uint8_t GetMotorID();
void MotorStatus(uint8_t motorID,uint8_t motorCnt);
void ChangeMotorID(uint8_t NewMotorID);


/*Current control functions */
	void  CurrentMode(uint8_t motorID);
	static int16_t CurrentToValue(float current);
	void sendCurrentCommand(uint8_t motorID, float current);

/*Velocity control functions*/
	void  VelocityMode(uint8_t motorID);
	static int16_t VelocityToValue(float velocity);
	void sendVelocityCommand(uint8_t motorID, float velocity);
	void MotorStop(uint8_t motorID);

/*Position control functions*/
	void  PositionMode(uint8_t motorID);
	static uint16_t AngleToValue(float angle_deg);
	void sendPositionCommand(uint8_t motorID, float angle_deg);


	void MotorSwitchMode(uint8_t motorID,uint8_t mode);

	void MotorStatus(uint8_t motorID, uint8_t motorCnt); //For real-time motor status
	void MotorStatusOLD(uint8_t motorID, uint8_t motorCnt); // NOT WORKING!!
/*CRC8  MAXIM-DOW*/
	uint8_t compute_crc8(uint8_t *data, uint8_t len);


#endif /* INC_DDSM115_H_ */
