/*
 * Controllers.h
 *
 *  Created on: Apr 17, 2025
 *      Author: Aljaž Pfifer, Karik Meijer
 */

#ifndef INC_CONTROLLERS_H_
#define INC_CONTROLLERS_H_


#include "stdint.h"

void LQR_controller(uint8_t MOTOR_ID, uint8_t ControllerStatus);
void LQR_controllerLR(uint8_t MOTOR_ID, uint8_t ControllerStatus);
void LateralController(uint8_t MOTOR_ID, uint8_t ControllerStatus);

#endif /* INC_CONTROLLERS_H_ */
