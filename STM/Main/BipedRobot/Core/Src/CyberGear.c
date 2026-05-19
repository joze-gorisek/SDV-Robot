/*
 * CyberGear.c
 *
 *  Created on: Mar 21, 2025
 *      Author: Andy
 */

#include "main.h"
#include "CyberGear.h"

#include "string.h"
#include "stdbool.h"
#include "stdio.h"


/*Extern periferial unit*/
extern CAN_HandleTypeDef hcan1;
extern UART_HandleTypeDef huart2;
//void HAL_CAN_RxFifo0MsgPendingCallback(CAN_HandleTypeDef *hcan);


uint8_t received_data_cg[8]; //CAN data buffer
uint16_t angle_u_cg;      //Current motor angle
uint32_t CANReadSingleValue=0;

float MOTangle=0;      //Current angle    [0-65535] -> [-4pi  4pi]
float MOTvelocity=0;   //Current velocity [0-65535] -> [-30rad/s  30rad/s]
float MOTtorque=0;     //Current Torque   [0-65535] -> [-12Nm   12Nm]
float MOTtemp=0;       //Temperature  10*Celsius
float MOTrpm =0;       //Current RPM

extern float MOTangle100;      //Current angle    [0-65535] -> [-4pi  4pi]
extern float MOTvelocity100;   //Current velocity [0-65535] -> [-30rad/s  30rad/s]
extern float MOTtorque100;     //Current Torque   [0-65535] -> [-12Nm   12Nm]
extern float MOTtemp100;       //Temperature  10*Celsius
extern float MOTrpm100;       //Current RPM

extern float MOTangle10;      //Current angle    [0-65535] -> [-4pi  4pi]
extern float MOTvelocity10;   //Current velocity [0-65535] -> [-30rad/s  30rad/s]
extern float MOTtorque10;     //Current Torque   [0-65535] -> [-12Nm   12Nm]
extern float MOTtemp10;       //Temperature  10*Celsius
extern float MOTrpm10 ;       //Current RPM

extern float MOTangle9;      //Current angle    [0-65535] -> [-4pi  4pi]
extern float MOTvelocity9;   //Current velocity [0-65535] -> [-30rad/s  30rad/s]
extern float MOTtorque9;     //Current Torque   [0-65535] -> [-12Nm   12Nm]
extern float MOTtemp9;       //Temperature  10*Celsius
extern float MOTrpm9;       //Current RPM


extern float MOTangle11;      //Current angle    [0-65535] -> [-4pi  4pi]
extern float MOTvelocity11;   //Current velocity [0-65535] -> [-30rad/s  30rad/s]
extern float MOTtorque11;     //Current Torque   [0-65535] -> [-12Nm   12Nm]
extern float MOTtemp11;       //Temperature  10*Celsius
extern float MOTrpm11;       //Current RPM

#define P_MIN -12.5f
#define P_MAX 12.5f
#define V_MIN -30.0f
#define V_MAX 30.0f
#define KP_MIN 0.0f
#define KP_MAX 500.0f
#define KD_MIN 0.0f
#define KD_MAX 5.0f
#define T_MIN -12.0f
#define T_MAX 12.0f

int motorID_Prekinitev = 0;
/**
 * Ask the motor for its device ID => type=0
 * The motor (if it hears this) should respond with type=0,
 * data=64-bit unique MCU ID, or an info frame. CANB
 */
void getMotorDeviceID(uint8_t hostID, uint8_t motorID)
{
    CAN_TxHeaderTypeDef txHeader;
    uint32_t txMailbox;
    uint8_t txData[8] = {0};

    // 0 in bits28..24 => get device ID
    uint32_t extId = ((uint32_t)0 << 24)
                   | ((uint32_t)hostID << 8)
                   | (uint32_t)motorID;

    txHeader.ExtId = extId;
    txHeader.IDE   = CAN_ID_EXT;     // extended frame
    txHeader.RTR   = CAN_RTR_DATA;
    txHeader.DLC   = 8;
    txHeader.TransmitGlobalTime = DISABLE;

    // Typically data can be all 0
    memset(txData , 0, 8);

    HAL_CAN_AddTxMessage(&hcan1, &txHeader, txData, &txMailbox);
}



/**
 * Ask the motor for its device ID => type=0
 * The motor (if it hears this) should respond with type=0,
 * data=64-bit unique MCU ID, or an info frame. CANB
 */

void clearMotorFault(uint8_t hostID, uint8_t motorID)
{
    CAN_TxHeaderTypeDef txHeader;
    uint32_t txMailbox;
    uint8_t txData[8] = {0};

    // 4 in bits28..24 => Stop command
    uint32_t extId = ((uint32_t)4 << 24)
                   | ((uint32_t)hostID << 8)
                   | (uint32_t)motorID;

    txHeader.ExtId = extId;
    txHeader.IDE   = CAN_ID_EXT;     // extended frame
    txHeader.RTR   = CAN_RTR_DATA;
    txHeader.DLC   = 8;
    txHeader.TransmitGlobalTime = DISABLE;

    // data[0] = 1 => clear fault
    txData[0] = 1;

    HAL_CAN_AddTxMessage(&hcan1, &txHeader, txData, &txMailbox);
}



/**
 * Write a 4-byte parameter at index paramIndex (e.g. 0x7005 for run_mode).
 *   paramValue must point to 4 bytes (e.g. float).
 *   type=18 (0x12) in bits28..24
 *
 *
 *   Reply frame: Reply motor feedback frame (see communication type 2) WITH ALL motor parameters
 */
void writeParameter(uint16_t paramIndex, const volatile void* paramValue,
                                 uint8_t hostID, uint8_t motorID)

{
    CAN_TxHeaderTypeDef txHeader;
    uint32_t txMailbox;
    uint8_t txData[8] = {0};

    // Build extended ID => (type=18)
    uint32_t extId = ((uint32_t)0x12 << 24)
                   | ((uint32_t)hostID << 8)
                   | (uint32_t)motorID;

    txHeader.ExtId = extId;
    txHeader.IDE   = CAN_ID_EXT;
    txHeader.RTR   = CAN_RTR_DATA;
    txHeader.DLC   = 8;
    txHeader.TransmitGlobalTime = DISABLE;

    // Byte0..1 = paramIndex (little-endian)
    txData[0] = (uint8_t)(paramIndex & 0xFF);
    txData[1] = (uint8_t)(paramIndex >> 8);
    // Byte2..3 = 0
    // Byte4..7 = paramValue
    memcpy(&txData[4], paramValue, 4);

    HAL_CAN_AddTxMessage(&hcan1, &txHeader, txData, &txMailbox);
}
/**
 * Write a 4-byte parameter at index paramIndex (e.g. 0x7005 for run_mode).
 *   paramValue must point to 4 bytes (e.g. float).
 *   type=18 (0x12) in bits28..24
 */
void readParameter(uint16_t paramIndex, uint8_t hostID, uint8_t motorID)
{
    CAN_TxHeaderTypeDef txHeader;
    uint32_t txMailbox;
    uint8_t txData[8] = {0};

    //Single parameter read  ID  (communication type 17  -> hex 0x11)
    uint32_t extId = ((uint32_t)0x11 << 24) //Single parameter read
                   | ((uint32_t)hostID << 8)
                   | (uint32_t)motorID;

    txHeader.ExtId = extId;
    txHeader.IDE   = CAN_ID_EXT;
    txHeader.RTR   = CAN_RTR_DATA;
    txHeader.DLC   = 8;
    txHeader.TransmitGlobalTime = DISABLE;


    // paramValue
    // 0x7019 - mechPos (angle - rad)
    // 0x701B - mechVel (speed -30..30 rad/s)
    // 0x701A - iq      (Current -23..23 A)-Torque


    // Byte0..1 = paramIndex (little-endian)
    txData[0] = (uint8_t)(paramIndex & 0xFF);
    txData[1] = (uint8_t)(paramIndex >> 8);
    txData[2] = 0;
    txData[3] = 0;
    txData[4] = 0;
    txData[5] = 0;
    txData[6] = 0;
    txData[7] = 0;
    // Byte2..3 = 0
    // Byte4..7 = paramValue
    //memcpy(&txData[2], 0, 6);

     HAL_CAN_AddTxMessage(&hcan1, &txHeader, txData, &txMailbox);
}


/*
 * Motor control mode 1-positon 2-velocity 3-current
 *
 * Control_mode -  (0 - Operation control mode) (1 - Position mode) (2 - Speed mode) (3 - Current mode)
 * */

void MotorControlMode(uint8_t Control_mode, uint8_t hostID, uint8_t motorID)
{

	writeParameter(0x7005, &Control_mode, hostID, motorID);

}

/*
 *
 * POSITON MODE FUNCTIONS
 *
 * */

void PositionSpeedLimit(float SpeedLimit, uint8_t hostID, uint8_t motorID)
{
	// 0x7017 - Position mode speed limit: SpeedLimit - 0 ~ 30rad/s   (Communication type=18)
	writeParameter(0x7017, &SpeedLimit, hostID, motorID);

}

void SetAngle(float Angle, uint8_t hostID, uint8_t motorID) //Works only in POSITION MODE!
{
	// 0x7016 - Position mode angle command:  Angle[rad]   (Communication type=18)
	 writeParameter(0x7016, &Angle, hostID, motorID); //Position
}

void ReadAngle(uint8_t hostID, uint8_t motorID)
{
	//0x7019  Load end lap counting mechanical angle [RAD]
	readParameter(0x7019,hostID, motorID); //Data in CAN interrupr
}

void ReadRotation(uint8_t hostID, uint8_t motorID)
{
	//0x701D  Number of turns [int16]
	readParameter(0x701D,hostID, motorID); //Data in CAN interrupr
}

/*
 *
 * VELOCITY MODE FUNCTIONS
 *
 * */

//Current limit for Velocity control
void VelocityCurrentLimit(float CurrentLimit, uint8_t hostID, uint8_t motorID) //Current limiter for Velocity mode
{
	//0x7018 Speed Position Mode Current Limit:  0 ~ 23A
	writeParameter(0x7018, &CurrentLimit, hostID, motorID);
}

void SetVelocity(float Velocity, uint8_t hostID, uint8_t motorID) //Works only in VELOCITY MODE!
{
	// 0x700A - PSpeed mode speed command Velocity: -30 ~ 30rad/s  (Communication type=18)
	 writeParameter(0x700A, &Velocity, hostID, motorID); //Position
}

void ReadVelocity(uint8_t hostID, uint8_t motorID)
{
	//0x701B (mechVel) Load end speed 	-30 ~ 30rad/s
	readParameter(0x701B,hostID, motorID); //Data in CAN interrupr
}


/*
 *
 * CURRENT-TORQUE MODE FUNCTIONS
 *
 * */

void TorqueLimit(float TorqueLimit, uint8_t hostID, uint8_t motorID) //Current limiter for Velocity mode
{
	//0x700B Torque limit:  0~12Nm
	writeParameter(0x700B, &TorqueLimit, hostID, motorID);
}

//Set TOURQE - Iq current
void SetIq(float Iq, uint8_t hostID, uint8_t motorID) //Works only in VELOCITY MODE!
{
	// 0x7006 - iq_ref Current Mode Iq Command Iqy: -23 ~ 23A (Communication type=18)
	 writeParameter(0x7006, &Iq, hostID, motorID); //Position
}

void ReadIq(uint8_t hostID, uint8_t motorID)
{
	//0x701A (iq) filter value 	-23 ~ 23A
	readParameter(0x701A,hostID, motorID); //Data in CAN interrupr
}

void ReadVBUS(uint8_t hostID, uint8_t motorID)
{
	//0x701C VBUS	bus voltage V
	readParameter(0x701C,hostID, motorID); //Data in CAN interrupr
}
/**
 * Enable motor => type=3
 */
void motorEnable(uint8_t hostID, uint8_t motorID)
{
    CAN_TxHeaderTypeDef txHeader;
    uint8_t txData[8] = {0};
    uint32_t txMailbox;

    uint32_t extId = ((uint32_t)3 << 24) |
                     ((uint32_t)hostID << 8) |
                     motorID;

    txHeader.ExtId = extId;
    txHeader.IDE   = CAN_ID_EXT;
    txHeader.RTR   = CAN_RTR_DATA;
    txHeader.DLC   = 8;
    txHeader.TransmitGlobalTime = DISABLE;

    HAL_CAN_AddTxMessage(&hcan1, &txHeader, txData, &txMailbox);
}

/**
 * Stop motor => type=4
 */
void motorStop(uint8_t hostID, uint8_t motorID)
{
    CAN_TxHeaderTypeDef txHeader;
    uint8_t txData[8] = {0};
    uint32_t txMailbox;

    uint32_t extId = ((uint32_t)4 << 24) |
                     ((uint32_t)hostID << 8) |
                     motorID;

    txHeader.ExtId = extId;
    txHeader.IDE   = CAN_ID_EXT;
    txHeader.RTR   = CAN_RTR_DATA;
    txHeader.DLC   = 8;
    txHeader.TransmitGlobalTime = DISABLE;

    // data[0] = 1 => if you want to clear faults, else 0
    // txData[0] = 1;

     HAL_CAN_AddTxMessage(&hcan1, &txHeader, txData, &txMailbox);
}


/*Set mechanical zero - Encoder data*/
void setMechanicalZero(uint8_t hostID, uint8_t motorID) {
    CAN_TxHeaderTypeDef txHeader;
    uint8_t txData[8] = {0};
    uint32_t txMailbox;

    // Communication type 6 is used to set mechanical zero.
    // Build the extended ID: type (6) is in bits28..24,
    // hostID is in bits15..8 and motorID in bits7..0.
    uint32_t extId = ((uint32_t)6 << 24) | ((uint32_t)hostID << 8) | motorID;

    txHeader.ExtId = extId;
    txHeader.IDE   = CAN_ID_EXT;
    txHeader.RTR   = CAN_RTR_DATA;
    txHeader.DLC   = 8;
    txHeader.TransmitGlobalTime = DISABLE;

    // Set data byte 0 to 1 to indicate zeroing the angle.
    txData[0] = 1;

 HAL_CAN_AddTxMessage(&hcan1, &txHeader, txData, &txMailbox);
}


/*
 * SET mew motor ID - take effect immediately
 * */
void setNewID(uint8_t hostID,uint8_t old_motorID, uint8_t new_motorID)
{

    CAN_TxHeaderTypeDef txHeader;
    uint8_t txData[8] = {0};
    uint32_t txMailbox;

    // Communication type 6 is used to set mechanical zero.
    // Build the extended ID: type (6) is in bits28..24,
    // hostID is in bits15..8 and motorID in bits7..0.
    uint32_t extId = ((uint32_t)7 << 24)          |
    		         ((uint32_t)hostID << 16)     |
					 ((uint32_t)new_motorID << 8) |
					 old_motorID;

    txHeader.ExtId = extId;
    txHeader.IDE   = CAN_ID_EXT;
    txHeader.RTR   = CAN_RTR_DATA;
    txHeader.DLC   = 8;
    txHeader.TransmitGlobalTime = DISABLE;

    // Set data byte 0 to 1 to indicate zeroing the angle.
    txData[0] = 1;

     HAL_CAN_AddTxMessage(&hcan1, &txHeader, txData, &txMailbox);


}

void ReadAllMotorData(uint8_t hostID, uint8_t motorID)
{
	//ReadAngle( hostID, motorID);
	//MOTangle= (float)(CANReadSingleValue);     //- 32767 //Current angle    [0-65535] -> [-4pi  4pi]

	//HAL_Delay(10);
	ReadVelocity(hostID, motorID);
    //MOTvelocity=(float)(CANReadSingleValue);   //Current velocity [0-65535] -> [-30rad/s  30rad/s]

    //HAL_Delay(10);
    ReadIq(hostID, motorID);
    //MOTtorque=(float)(CANReadSingleValue);     //Current Torque   [0-65535] -> [-12Nm   12Nm]
    ReadAngle( hostID, motorID);
}



void Motor_MITrun(uint8_t hostID, uint8_t motorID, float torque, float MechPosition, float speed, float kp, float kd)
{

    CAN_TxHeaderTypeDef txHeader;
    uint8_t txData[8] = {0};
    uint32_t txMailbox;

    uint32_t extId = ((uint32_t)1 << 24) |
                     ((uint32_t)hostID << 8) |
                     motorID;

    txHeader.ExtId = extId;
    txHeader.IDE   = CAN_ID_EXT;
    txHeader.RTR   = CAN_RTR_DATA;
    txHeader.DLC   = 8;
    txHeader.TransmitGlobalTime = DISABLE;


    // Byte0..1 = paramIndex (little-endian)
    txData[0] = float_to_uint(MechPosition,P_MIN,P_MAX)>>8;
    txData[1] = float_to_uint(MechPosition,P_MIN,P_MAX);
    txData[2] = float_to_uint(speed,V_MIN,V_MAX)>>8;
    txData[3] = float_to_uint(speed,V_MIN,V_MAX);
    txData[4] = float_to_uint(kp,KP_MIN,KP_MAX)>>8;
    txData[5] = float_to_uint(kp,KP_MIN,KP_MAX);
    txData[6] = float_to_uint(kd,KD_MIN,KD_MAX)>>8;
    txData[7] = float_to_uint(kd,KD_MIN,KD_MAX);
    // Byte2..3 = 0
    // Byte4..7 = paramValue
    //memcpy(&txData[2], 0, 6);

     HAL_CAN_AddTxMessage(&hcan1, &txHeader, txData, &txMailbox);





}




/*
 *
 *
 * Motor response over CAN interrupt
 *
 * */
void HAL_CAN_RxFifo0MsgPendingCallback(CAN_HandleTypeDef *hcan)
{
    HAL_StatusTypeDef status;
    CAN_RxHeaderTypeDef pRxHeader;
    status = HAL_CAN_GetRxMessage(hcan, CAN_RX_FIFO0, &pRxHeader, received_data_cg);
    if (status != HAL_OK)
    {
        // Error reading the message
        return;
    }

    uint32_t extId = pRxHeader.ExtId;
    uint8_t type   = (extId >> 24) & 0x1F;       // bits28..24
    uint8_t motorID= (extId>>8) & 0xFF;          // OFFSET 0xFE  (bits7..0)  MOTOR ID (bits23..8)
    motorID_Prekinitev = motorID;
    // For debug, print the entire extended ID
   // char dbg[128];
   // sprintf(dbg, "RX: ExtId=0x%08lX (type=%u, motorID=%u), DLC=%u\n\r",
    //        (unsigned long)extId, (unsigned)type, (unsigned)motorID, pRxHeader.DLC);
    //serialWrite(dbg);

    // Check what type we received
    if (type == 2)  //RESPONSE WHEM WRITE SINGLE PARAMETER%
    {
        // This is typically the motor feedback frame
        // Byte0..1 might be current angle, Byte2..3 velocity, Byte4..5 torque, Byte6..7 temperature, etc.
        //sprintf(dbg, "Received type=2 (feedback) from motorID=%u\n\r", (unsigned)motorID);
      //  serialWrite(dbg);

        // Optionally parse the data
        // e.g. parse angle:
        angle_u_cg = (received_data_cg[0] << 8) | received_data_cg[1];

        MOTangle  =((float)(((received_data_cg[0] << 8) | received_data_cg[1])) * 8 * PI/65535.0f) - 4 * PI;    //Current angle    [0-65535] -> [-4pi  4pi]
        MOTvelocity=((float)(((received_data_cg[2] << 8) | received_data_cg[3])) * 2 * V_MAX/65535.0f) - V_MAX;   //Current velocity [0-65535] -> [-30rad/s  30rad/s]
        MOTtorque  =((float)(((received_data_cg[4] << 8) | received_data_cg[5])) * 2 * T_MAX/65535.0f) - T_MAX;   //Current Torque   [0-65535] -> [-12Nm   12Nm]
        MOTtemp    =(float)(((received_data_cg[6] << 8) | received_data_cg[7])) / 10.0f ;        //Temperature  10*Celsius


        MOTrpm= MOTvelocity/(2*PI)*60.0; //Convert to RPM



  	  switch(motorID_Prekinitev) {
  	  case 100:
  		  MOTangle100 = MOTangle;
  		  MOTvelocity100=MOTvelocity;
  		  MOTtorque100=MOTtorque;
  		  MOTtemp100=MOTtemp;
  		  MOTrpm100 =MOTrpm;
  		  break;
  	  case 9:
  		  MOTangle9 = MOTangle;
  		  MOTvelocity9=MOTvelocity;
  		  MOTtorque9=MOTtorque;
  		  MOTtemp9=MOTtemp;
  		  MOTrpm9 =MOTrpm;
  		  break;
  	  case 11:
  		  MOTangle11 = MOTangle;
  		  MOTvelocity11=MOTvelocity;
  		  MOTtorque11=MOTtorque;
  		  MOTtemp11=MOTtemp;
  		  MOTrpm11 =MOTrpm;
  		  break;
  	  case 10:
  		  MOTangle10 = MOTangle;
  		  MOTvelocity10=MOTvelocity;
  		  MOTtorque10=MOTtorque;
  		  MOTtemp10=MOTtemp;
  		  MOTrpm10 =MOTrpm;
  		  break;

  	  }
        // etc. Then print or store
    }
    else if (type == 21)
    {
        // This is typically the fault/error frame
        // Byte0..3 might be a fault code
        //sprintf(dbg, "Received type=21 (fault) from motorID=%u\n\r", (unsigned)motorID);
        //serialWrite(dbg);

        // Optionally parse the fault code in received_data[0..3]
        // For example:
        uint32_t fault = (received_data_cg[0])
                       | (received_data_cg[1] << 8)
                       | (received_data_cg[2] << 16)
                       | (received_data_cg[3] << 24);
        //sprintf(dbg, "Fault code = 0x%08lX\n", (unsigned long)fault);
        //serialWrite(dbg);
    }
    else if (type == 0)
    {
        // Possibly a "Get Device ID" response
        //sprintf(dbg, "Received type=0 (device info) from motorID=%u\n\r", (unsigned)motorID);
        //serialWrite(dbg);
        // The data bytes may contain the 64-bit unique ID
    }else if (type == 17)
    {
        //Response from single parameter read

    	  CANReadSingleValue= (received_data_cg[4])
                      | (received_data_cg[5] << 8)
                      | (received_data_cg[6] << 16)
                      | (received_data_cg[7] << 24) & 0xFFFF;


    	  float data_test;
    	 memcpy (&data_test,  received_data_cg+4, sizeof (float));

    	  if( (received_data_cg[0] | (received_data_cg[1] << 8) )== 0x701B)//Velocity
    		  MOTvelocity= uint_to_float( CANReadSingleValue, -30.0f, 30.0f);
    	  else if( (received_data_cg[0] | (received_data_cg[1] << 8)) == 0x7019)//Position
    	       MOTangle = uint_to_float( CANReadSingleValue, -12.56f, 12.56f);
    	  else if((received_data_cg[0] | (received_data_cg[1] << 8)) == 0x701A)//Current
    	      MOTtorque= uint_to_float( CANReadSingleValue, -23.00f, 23.00f);



    }
    else
    {
        // Some other type
        //sprintf(dbg, "Received unknown type=%u from motorID=%u\n\r",
         //       (unsigned)type, (unsigned)motorID);
        //serialWrite(dbg);
    }
}

float uint_to_float(uint32_t x, float x_min, float x_max){
    uint32_t type_max = 0xFFFF;
    float span = x_max - x_min;
    return (float) (x / (type_max * span)) + x_min;
}
/*
 * Function used for MIT control.
 */
static uint16_t float_to_uint(float x, float x_min, float x_max)
{
    // Clamp x into [x_min, x_max]
    if(x < x_min) x = x_min;
    if(x > x_max) x = x_max;

    float span = (x_max - x_min);
    float offset = (x - x_min);

    // Scale into [0..65535]
    // The Xiaomi snippet does:
    //   (x - offset)*(65535.0f)/span
    // but we can do something like:
    return (uint16_t)((offset * 65535.0f) / span + 0.5f);
}


/*
 * UART functions
 * */

void serialWrite(char data[]){
	HAL_UART_Transmit(&huart2, (uint8_t *) data, strlen(data), 10);
	HAL_UART_Transmit(&huart2,(uint8_t *)"\n",1,10);
}

void serialProcessRxData(){

}

