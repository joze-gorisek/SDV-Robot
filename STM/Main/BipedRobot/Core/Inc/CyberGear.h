/*
 * CyberGear.h
 *
 *  Created on: Mar 21, 2025
 *      Author: Andy
 */

#ifndef INC_CYBERGEAR_H_
#define INC_CYBERGEAR_H_

#include "stm32f4xx_hal.h"  // This includes the HAL definitions including HAL_StatusTypeDef.
#include "stdbool.h"
#include "stdint.h"
#include "stdlib.h"

#define P_MIN  -12.5f     // min angle (or -4π or some doc example)
#define P_MAX  12.5f      // max angle
#define V_MIN  -30.0f
#define V_MAX   30.0f
#define KP_MIN  0.0f
#define KP_MAX  500.0f
#define KD_MIN  0.0f
#define KD_MAX  5.0f
#define T_MIN  -12.0f
#define T_MAX   12.0f
#define PI      3.14159f


/* MOTOR ID*/
void getMotorDeviceID(uint8_t hostID, uint8_t motorID);

void motorEnable(uint8_t hostID, uint8_t motorID);
void motorStop(uint8_t hostID, uint8_t motorID);

/* MOTOR FEEDBACK CONTROL*/
void MotorControlMode(uint8_t Control_mode,
		              uint8_t hostID, uint8_t motorID);  // Control_mode 1-positon 2-velocity 3-current



/*MIT mode*/
void Motor_MITrun(uint8_t hostID, uint8_t motorID, float torque, float MechPosition, float speed, float kp, float kd);

/*Position mode*/
void PositionSpeedLimit(float SpeedLimit, uint8_t hostID, uint8_t motorID);
void SetAngle(float Angle, uint8_t hostID, uint8_t motorID);
void ReadAngle(uint8_t hostID, uint8_t motorID);

/*Velocity mode*/
void VelocityCurrentLimit(float CurrentLimit, uint8_t hostID, uint8_t motorID); //Current limiter for Velocity mode
void SetVelocity(float Velocity, uint8_t hostID, uint8_t motorID); //Works only in Velocity/Speed MODE!
void ReadVelocity(uint8_t hostID, uint8_t motorID);

/*Current mode*/
void TorqueLimit(float TorqueLimit, uint8_t hostID, uint8_t motorID); //Current limiter for Velocity mode
void SetIq(float Iq, uint8_t hostID, uint8_t motorID); //Works only in VELOCITY MODE!
void ReadIq(uint8_t hostID, uint8_t motorID);



/*Read write data*/
void readParameter(uint16_t paramIndex, uint8_t hostID, uint8_t motorID);
void writeParameter(uint16_t paramIndex, const volatile void* paramValue,
                                 uint8_t hostID, uint8_t motorID);

void clearMotorFault(uint8_t hostID, uint8_t motorID);
void ReadAllMotorData(uint8_t hostID, uint8_t motorID); //Read Iq, Velocity and Angle

void setMechanicalZero(uint8_t hostID, uint8_t motorID);
static uint16_t float_to_uint(float x, float x_min, float x_max);
void setNewID(uint8_t hostID,uint8_t old_motorID, uint8_t new_motorID);



/*Motor movements*/
void CyberFrontAngle(void);
void CyberBackAngle(void);
void CyberUpDown(void);
void SpeedLimit(void);
void CyberRollLeft(void);
void CyberRollRight(void);
void CyberDummyRequest(void);


//STM32 Communication functions
void serialWrite(char data[]);
void serialProcessRxData();
//void HAL_CAN_RxFifo0MsgPendingCallback(CAN_HandleTypeDef *hcan);
//void HAL_UART_RxCpltCallback(UART_HandleTypeDef *huart);

float uint_to_float(uint32_t x, float x_min, float x_max);
#endif /* INC_CYBERGEAR_H_ */
