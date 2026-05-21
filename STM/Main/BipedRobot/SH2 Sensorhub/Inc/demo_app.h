/*
 * demo_app.h
 *
 *  Created on: Apr 7, 2026
 *      Author: blpongrac
 */

#ifndef INC_DEMO_APP_H_
#define INC_DEMO_APP_H_

#include "main.h"
#include "sh2.h"
#include "sh2_SensorValue.h"
#include "sh2_err.h"

#define ROBOT_DATA_SIZE                     20U

void BNO_App(float robot_data[ROBOT_DATA_SIZE]);
void BNO_EXTI_Callback(uint16_t GPIO_Pin);

// BNO086 SparkFun breakout: ADDR pin low = 0x4A, high = 0x4B
#define BNO_ADDR  (0x4B << 1)

void BNO_Init(I2C_HandleTypeDef *hi2c, UART_HandleTypeDef *huart, uint16_t intPin);

// Rotation results (written by sensor callback, read by application)
extern volatile float bno_roll, bno_pitch, bno_yaw;   // degrees, from quaternion
extern volatile float bno_qw, bno_qx, bno_qy, bno_qz; // raw quaternion
extern volatile float bno_gx, bno_gy, bno_gz;

#endif /* INC_DEMO_APP_H_ */
