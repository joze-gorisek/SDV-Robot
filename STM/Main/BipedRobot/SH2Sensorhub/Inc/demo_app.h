/*
 * demo_app.h
 *
 *  Created on: Apr 7, 2026
 *      Author: blpongrac
 */

#ifndef INC_DEMO_APP_H_
#define INC_DEMO_APP_H_

#include "main.h"

#include "../../SH2Sensorhub/Inc/sh2.h"
#include "../../SH2Sensorhub/Inc/sh2_err.h"
#include "../../SH2Sensorhub/Inc/sh2_SensorValue.h"

// dodano
#ifndef BNO_APP_H
#define BNO_APP_H
void BNO_App(float vector[4]);
#endif
// dodano

// BNO086 SparkFun breakout: ADDR pin low = 0x4A, high = 0x4B
#define BNO_ADDR  (0x4B << 1)

void BNO_Init(I2C_HandleTypeDef *hi2c, UART_HandleTypeDef *huart, uint16_t intPin);

// Rotation results (written by sensor callback, read by application)
extern volatile float bno_roll, bno_pitch, bno_yaw;   // degrees, from quaternion
extern volatile float bno_qw, bno_qx, bno_qy, bno_qz; // raw quaternion
extern volatile float bno_gx, bno_gy, bno_gz;

#endif /* INC_DEMO_APP_H_ */
