/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2022 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include "string.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <ctype.h>
#include "DDSM115.h"
#include "mpu6050.h"
#include "CyberGear.h"
#include "MRF24J40.h"
#include "demo_app.h"
#include <math.h>
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */

#define DDSM_LEFT_ID      0x10U
#define DDSM_RIGHT_ID     0x30U
#define DDSM_PENDING_NONE 0x00U
#define DDSM_FEEDBACK_TIMEOUT_MS 10U
#define DEBUG_PRINT_INTERVAL_MS 20U
#define BUTTON_DEBOUNCE_MS 300U
#define DEG_TO_RAD (PI / 180.0f)
#define MRF_COMMAND_MAX_LEN 32U
#define MRF_UP_TARGET_STEP 0.7f
#define MRF_DOWN_TARGET_STEP 0.0f
#define MRF_POSITION_RAMP_SPEED_MPS 0.50f
#define MRF_MAX_COMMAND_DISTANCE_M 1.00f

/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
CAN_HandleTypeDef hcan1;

I2C_HandleTypeDef hi2c3;

SPI_HandleTypeDef hspi2;

TIM_HandleTypeDef htim6;

UART_HandleTypeDef huart5;
UART_HandleTypeDef huart2;
DMA_HandleTypeDef hdma_uart5_rx;

/* USER CODE BEGIN PV */
uint8_t i = 0;

CAN_TxHeaderTypeDef pTxHeader;
CAN_RxHeaderTypeDef pRxHeader;
CAN_FilterTypeDef sFilterConfig;
uint32_t pTxMailbox;

//extern void FrontAngle(void);
//extern void BackAngle(void);


extern float MOTangle;      //Current angle    [0-65535] -> [-4pi  4pi]
extern float MOTvelocity;   //Current velocity [0-65535] -> [-30rad/s  30rad/s]
extern float MOTtorque;     //Current Torque   [0-65535] -> [-12Nm   12Nm]
extern float MOTtemp;       //Temperature  10*Celsius
extern float MOTrpm;       //Current RPMe

float MOTangle100=0;      //Current angle    [0-65535] -> [-4pi  4pi]
float MOTvelocity100=0;   //Current velocity [0-65535] -> [-30rad/s  30rad/s]
float MOTtorque100=0;     //Current Torque   [0-65535] -> [-12Nm   12Nm]
float MOTtemp100=0;       //Temperature  10*Celsius
float MOTrpm100 =0;       //Current RPM

float MOTangle10=0;      //Current angle    [0-65535] -> [-4pi  4pi]
float MOTvelocity10=0;   //Current velocity [0-65535] -> [-30rad/s  30rad/s]
float MOTtorque10=0;     //Current Torque   [0-65535] -> [-12Nm   12Nm]
float MOTtemp10=0;       //Temperature  10*Celsius
float MOTrpm10 =0;       //Current RPM

float MOTangle9=0;      //Current angle    [0-65535] -> [-4pi  4pi]
float MOTvelocity9=0;   //Current velocity [0-65535] -> [-30rad/s  30rad/s]
float MOTtorque9=0;     //Current Torque   [0-65535] -> [-12Nm   12Nm]
float MOTtemp9=0;       //Temperature  10*Celsius
float MOTrpm9 =0;       //Current RPM

volatile uint8_t button_step_requested = 0;
float robot_data[ROBOT_DATA_SIZE] = {0};
float state_data_L[4] = {0};
float state_data_R[4] = {0};

float MOTangle11=0;      //Current angle    [0-65535] -> [-4pi  4pi]
float MOTvelocity11=0;   //Current velocity [0-65535] -> [-30rad/s  30rad/s]
float MOTtorque11=0;     //Current Torque   [0-65535] -> [-12Nm   12Nm]
float MOTtemp11=0;       //Temperature  10*Celsius
float MOTrpm11 =0;       //Current RPM

/*Drive Wheels*/
extern float DDSangle01;       //Current angle    [0-32767] -> [0 360]   address= 0x10 over RS485
extern float DDSrpm01;         //Current velocity [0-     ] -> [-330RPM  330RPM]
extern float DDScurrent01;     //Current Torque   [-32767   32767] -> [-8A   8A]
extern float DDSvelocityRadial01;    //Velocity
extern float DDSvelocity01;          //Translational velocity [m/s]
extern float distnacex01;            //Translational distance [m]
extern float revolutionx01;
extern float zerocrossx01up;
extern float zerocrossx01down;

extern float DDSangle30;       //Current angle    [0-32767] -> [0 360]   address= 0x30 over RS485
extern float DDSrpm30;         //Current velocity [0-     ] -> [-330RPM  330RPM]
extern float DDScurrent30;     //Current Torque   [-32767   32767] -> [-8A   8A]
extern float DDSvelocityRadial30;    //Velocity
extern float DDSvelocity30;          //Translational velocity [m/s]
extern float distnacex30;            //Translational distance [m]
extern float revolutionx30;
extern float zerocrossx30up;
extern float zerocrossx30down;

float desired_angle ;
float desired_angle2 ;
float desired_angle3 ;
float desired_angle4 ;
float ddsm_left_current_cmd = 0.0f;
float ddsm_right_current_cmd = 0.0f;
volatile uint8_t ddsm_pending_motor_id = DDSM_PENDING_NONE;
volatile uint32_t ddsm_left_feedback_count = 0U;
volatile uint32_t ddsm_right_feedback_count = 0U;
volatile uint32_t ddsm_feedback_timeout_count = 0U;
static uint8_t ddsm_next_motor_id = DDSM_LEFT_ID;
static uint32_t ddsm_pending_tick = 0U;
static uint32_t debug_last_tick = 0U;
static uint8_t cyber_next_motor_index = 0U;

extern float front_angle;
extern float back_angle;

float targetValue_100=0.0;
float targetValue_10=0.0;
float targetValue_9=0.0;
float targetValue_11=0.0;

float speedLimit_100 = 1.0f;
//float speedLimit_10  = 1.0f;
//float speedLimit_9   = 1.0f;
//float speedLimit_11  = 1.0f;

// Tuning parameters MIT (Adjust these carefully!)

static float target = 0.0f;
static int direction = 1;

const float amplitude = 0.0f;
float target_step = 0.1f;

float torque = 0.0f;
float velocity = 0.0f;
const float kp = 3.0f;
const float kd = 0.5f;

extern Kalman_t KalmanX;

extern int motorID_Prekinitev;

extern uint8_t HEX_Buffer[20];

int8_t Buffer[40];

//USART2 DMA RX
#define RX_BUFFER_SIZE 10 // Adjust based on your needs
#define DELIMITER '%'      // Define the delimiter character

uint8_t rxBuffer[RX_BUFFER_SIZE]; // DMA receive buffer
uint8_t rxData[RX_BUFFER_SIZE];   // Processed message buffer

// initial target angle. Value should be 6.28... for a full rotation
float targetValue = 0.0f;
float iq_test=0.0;
extern uint8_t VCPBuffer[10];


uint8_t received_data[8];
// UART
char RxUARTBuffer[256]="";
uint8_t RxUARTLength=0;
uint8_t RxSingleByte;
uint8_t onFlag = 0;
uint8_t MOTOR_ID = 0x00;
char RxData[10];
// Global buffer for RS485 data reception
volatile uint8_t RS485_RxIndex = 0;

static float mrf_active_xd_L = 0.0f;
static float mrf_active_xd_R = 0.0f;
static float mrf_final_xd_L = 0.0f;
static float mrf_final_xd_R = 0.0f;
static volatile uint32_t mrf_exti_count = 0U;
static volatile uint8_t mrf_irq_pending = 0U;
static uint8_t mrf_regulation_active = 0U;
static uint8_t mrf_ramp_active = 0U;
static uint32_t mrf_ramp_last_tick = 0U;


extern double roll_kalman, pitch_kalman;

// We'll do a simple function that sets position mode and moves the motor



/* USER CODE END PV */

// Retarget printf to UART
int __io_putchar(int ch)
{
    HAL_UART_Transmit(&huart2, (uint8_t *)&ch, 1, HAL_MAX_DELAY);
    return ch;
}

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_DMA_Init(void);
static void MX_USART2_UART_Init(void);
static void MX_CAN1_Init(void);
static void MX_TIM6_Init(void);
static void MX_I2C3_Init(void);
static void MX_SPI2_Init(void);
static void MX_UART5_Init(void);
/* USER CODE BEGIN PFP */
//void HAL_CAN_RxFifo0MsgPendingCallback(CAN_HandleTypeDef *hcan);
void serialWrite(char data[]);
void serialProcessRxData();
//void HAL_UART_RxCpltCallback(UART_HandleTypeDef *huart);
static uint16_t float_to_uint(float x, float x_min, float x_max);
static HAL_StatusTypeDef UART5_StartReceiveDmaIfReady(void);
static void DDSM_Service(void);
static void DDSM_HandleFeedback(uint8_t motor_id);
static void CyberGear_Service(void);
static void RobotData_Update(void);
static void StateData_Update(void);
static void MRF_ServiceCommand(void);
static bool MRF_ParseMotionCommand(const uint8_t *data, uint8_t len, float *distance_m, uint8_t *up_bit);
static void MRF_ApplyCommand(float distance_m, uint8_t up_bit);
static void MRF_UpdatePositionRamp(uint32_t now);
static void MRF_PrintInitDiag(void);
static void MRF_BitBangPrimeRfRegs(void);
void HAL_UART_RxCpltCallback(UART_HandleTypeDef *huart);
//void HAL_UARTEx_RxEventCallback(UART_HandleTypeDef *huart, uint16_t Size);

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */


  //Variables from mpu6050.c
	MPU6050_t MPU6050;
	char buffer[100];
	extern double roll_kalman, pitch_kalman;


	//Variables from DDSM115.c
	extern int ControllerONLaterl;
    extern float current;
    extern float RPM;
	extern float angle;
	extern uint8_t RS485_RxBuffer[RS485_BUFFER_SIZE];
	extern uint8_t modeCmd[10];
	extern uint8_t position_mode[10];
	extern uint8_t command[10];
	extern uint8_t  ID_query[10];
	extern float Up_angle;
	extern float speedLimit; // CyberGear speed limit e.g. 1 rad/s  max 30rad/s
	//Variable from USART2 interrupt
	extern float VPCnum;

	//Variables MRF
	uint8_t read_intstat, read_intcon, read_RXMCR, read_BBREG1,read_RXFLUSH, frame_length;
	char buffer485[10];
	char Allocbuffer485[10];
	extern float anglex01_offset;
	extern float anglex30_offset;
	extern float delatx01;
	extern float delatx30;
	extern float DDSangle01_old ;
	extern float DDSangle30_old ;
	extern float DDSangle01m_old;
	extern float DDSangle30m_old;

static HAL_StatusTypeDef UART5_StartReceiveDmaIfReady(void)
{
	if (huart5.RxState == HAL_UART_STATE_READY)
	{
		return HAL_UART_Receive_DMA(&huart5, (uint8_t *)buffer485, sizeof(buffer485));
	}

	return HAL_BUSY;
}

static void RobotData_Update(void)
{
	robot_data[0] = distnacex01;             // DDSM left distance [m]
	robot_data[1] = distnacex30;             // DDSM right distance [m]
	robot_data[2] = DDSvelocity01;           // DDSM left speed [m/s]
	robot_data[3] = DDSvelocity30;           // DDSM right speed [m/s]

	robot_data[4] = MOTangle100;             // CyberGear motor 1 position [rad]
	robot_data[5] = MOTangle9;               // CyberGear motor 2 position [rad]
	robot_data[6] = MOTangle10;              // CyberGear motor 3 position [rad]
	robot_data[7] = MOTangle11;              // CyberGear motor 4 position [rad]
	robot_data[8] = DDSangle01 * DEG_TO_RAD; // DDSM left position [rad]
	robot_data[9] = DDSangle30 * DEG_TO_RAD; // DDSM right position [rad]
}

static void StateData_Update(void)
{
	state_data_L[0] = distnacex01;       // Robot position x [m]
	state_data_L[1] = DDSvelocity01;     // Robot speed v [m/s]
	state_data_L[2] = robot_data[10] * (3.141592653589793f / 180.0f);    // Robot angle theta, BNO roll [deg]
	state_data_L[3] = robot_data[13];    // Robot angular speed omega, BNO gyro X [rad/s]

	state_data_R[0] = distnacex30;        // Robot position x [m]
	state_data_R[1] = DDSvelocity30;      // Robot speed v [m/s]
	state_data_R[2] = state_data_L[2];    // Robot angle theta, BNO roll [deg]
	state_data_R[3] = state_data_L[3];    // Robot angular speed omega, BNO gyro X [rad/s]


}

static void DDSM_Service(void)
{
	uint8_t pending_motor = ddsm_pending_motor_id;

	if (pending_motor != DDSM_PENDING_NONE)
	{
		if ((HAL_GetTick() - ddsm_pending_tick) >= DDSM_FEEDBACK_TIMEOUT_MS)
		{
			ddsm_feedback_timeout_count++;
			ddsm_pending_motor_id = DDSM_PENDING_NONE;
			ddsm_next_motor_id = (pending_motor == DDSM_LEFT_ID) ? DDSM_RIGHT_ID : DDSM_LEFT_ID;
			HAL_UART_AbortReceive(&huart5);
			(void)UART5_StartReceiveDmaIfReady();
		}

		return;
	}

	uint8_t motor_id = ddsm_next_motor_id;
	float current_cmd = (motor_id == DDSM_LEFT_ID) ? ddsm_left_current_cmd : ddsm_right_current_cmd;

	ddsm_pending_motor_id = motor_id;
	ddsm_pending_tick = HAL_GetTick();
	sendCurrentCommand(motor_id, current_cmd);
	ddsm_next_motor_id = (motor_id == DDSM_LEFT_ID) ? DDSM_RIGHT_ID : DDSM_LEFT_ID;
}

static void DDSM_HandleFeedback(uint8_t motor_id)
{
	if (motor_id == DDSM_LEFT_ID)
	{
		ddsm_left_feedback_count++;
	}
	else if (motor_id == DDSM_RIGHT_ID)
	{
		ddsm_right_feedback_count++;
	}

	if (ddsm_pending_motor_id == motor_id)
	{
		ddsm_pending_motor_id = DDSM_PENDING_NONE;
	}
}

static void CyberGear_Service(void)
{
	if (HAL_CAN_GetTxMailboxesFreeLevel(&hcan1) == 0U)
	{
		return;
	}

	switch (cyber_next_motor_index)
	{
		case 0U:
			//SetAngle(desired_angle, CYBER_HOST_ID, CYBER_MOTOR_1_ID);
			Motor_MITrun(CYBER_HOST_ID, CYBER_MOTOR_1_ID, torque,  target_step, velocity, kp, kd);  //17 ID 24 -- plus dol
			break;
		case 1U:
			//SetAngle(desired_angle2, CYBER_HOST_ID, CYBER_MOTOR_2_ID);
			Motor_MITrun(CYBER_HOST_ID, CYBER_MOTOR_2_ID, torque,  -target_step, velocity, kp, kd);  //18 ID 22 -- pus gor
			break;
		case 2U:
			//SetAngle(desired_angle3, CYBER_HOST_ID, CYBER_MOTOR_3_ID);
			Motor_MITrun(CYBER_HOST_ID, CYBER_MOTOR_3_ID, torque,  -target_step, velocity, kp, kd);  //19 ID 21 -- pus gor
			break;
		default:
			//SetAngle(desired_angle4, CYBER_HOST_ID, CYBER_MOTOR_4_ID);
			Motor_MITrun(CYBER_HOST_ID, CYBER_MOTOR_4_ID, torque,  target_step, velocity, kp, kd);  //20 ID 23 -- NIKAMOR NE GRE!
			break;
	}

	cyber_next_motor_index++;
	if (cyber_next_motor_index >= 4U)
	{
		cyber_next_motor_index = 0U;
	}
}

static bool MRF_ParseMotionCommand(const uint8_t *data, uint8_t len, float *distance_m, uint8_t *up_bit)
{
	char cmd[MRF_COMMAND_MAX_LEN];

	if ((data == NULL) || (len == 0U) || (len >= MRF_COMMAND_MAX_LEN))
	{
		return false;
	}

	memcpy(cmd, data, len);
	cmd[len] = '\0';

	char *p = cmd;
	while (isspace((unsigned char)*p))
	{
		p++;
	}

	char *end = p;
	long distance_mm = strtol(p, &end, 10);
	if (end == p)
	{
		return false;
	}

	while (isspace((unsigned char)*end))
	{
		end++;
	}

	if (((end[0] != 'm') && (end[0] != 'M')) || ((end[1] != 'm') && (end[1] != 'M')))
	{
		return false;
	}
	end += 2;

	while (isspace((unsigned char)*end))
	{
		end++;
	}

	if ((end[0] != '0') && (end[0] != '1'))
	{
		return false;
	}
	*up_bit = (uint8_t)(end[0] - '0');
	end++;

	while (isspace((unsigned char)*end))
	{
		end++;
	}

	if (*end != '\0')
	{
		return false;
	}

	*distance_m = ((float)distance_mm) / 1000.0f;
	return true;
}

static void MRF_ServiceCommand(void)
{
	static uint8_t last_rx_count = 0U;

	if ((mrf_irq_pending != 0U) ||
		(HAL_GPIO_ReadPin(MRF_INT_GPIO_Port, MRF_INT_Pin) == GPIO_PIN_RESET))
	{
		mrf_irq_pending = 0U;
		Mrf24j_interrupt_handler();
	}

	if (flag_got_rx != last_rx_count)
	{
		last_rx_count = flag_got_rx;

		float distance_m = 0.0f;
		uint8_t up_bit = 0U;
		if (MRF_ParseMotionCommand(rx_info.rx_data, mrf_rx_data_len, &distance_m, &up_bit))
		{
			printf("MRF CMD distance=%.3f m up=%u\r\n", distance_m, up_bit);
			if (mrf_regulation_active != 0U)
			{
				MRF_ApplyCommand(distance_m, up_bit);
			}
			else
			{
				printf("MRF CMD ignored: regulation not started\r\n");
			}
		}
	}
}

static void MRF_ApplyCommand(float distance_m, uint8_t up_bit)
{
	if (fabsf(distance_m) > MRF_MAX_COMMAND_DISTANCE_M)
	{
		printf("MRF CMD ignored: distance too large\r\n");
		return;
	}

	mrf_final_xd_L = state_data_L[0] + distance_m;
	mrf_final_xd_R = state_data_R[0] + distance_m;
	target_step = (up_bit == 1U) ? MRF_UP_TARGET_STEP : MRF_DOWN_TARGET_STEP;
	mrf_ramp_active = 1U;
	mrf_ramp_last_tick = HAL_GetTick();

	printf("MRF TARGET final_L=%.3f final_R=%.3f target_step=%.2f\r\n",
		   mrf_final_xd_L,
		   mrf_final_xd_R,
		   target_step);
}

static void MRF_UpdatePositionRamp(uint32_t now)
{
	if (mrf_ramp_active == 0U)
	{
		mrf_ramp_last_tick = now;
		return;
	}

	uint32_t elapsed_ms = now - mrf_ramp_last_tick;
	mrf_ramp_last_tick = now;

	if (elapsed_ms == 0U)
	{
		return;
	}

	float max_step = MRF_POSITION_RAMP_SPEED_MPS * ((float)elapsed_ms / 1000.0f);
	float err_L = mrf_final_xd_L - mrf_active_xd_L;
	float err_R = mrf_final_xd_R - mrf_active_xd_R;

	if (fabsf(err_L) <= max_step)
	{
		mrf_active_xd_L = mrf_final_xd_L;
	}
	else
	{
		mrf_active_xd_L += (err_L > 0.0f) ? max_step : -max_step;
	}

	if (fabsf(err_R) <= max_step)
	{
		mrf_active_xd_R = mrf_final_xd_R;
	}
	else
	{
		mrf_active_xd_R += (err_R > 0.0f) ? max_step : -max_step;
	}

	if ((mrf_active_xd_L == mrf_final_xd_L) &&
		(mrf_active_xd_R == mrf_final_xd_R))
	{
		mrf_ramp_active = 0U;
	}
}

static void MRF_PrintInitDiag(void)
{
	printf("MAIN MRF PAN=0x%04X RFCON0=0x%02X RFCON6=0x%02X RXMCR=0x%02X INTCON=0x%02X BBREG1=0x%02X INTpin=%u\r\n",
		   Mrf24j_get_pan(),
		   Mrf24j_read_long(MRF_RFCON0),
		   Mrf24j_read_long(MRF_RFCON6),
		   Mrf24j_read_short(MRF_RXMCR),
		   Mrf24j_read_short(MRF_INTCON),
		   Mrf24j_read_short(MRF_BBREG1),
		   (unsigned int)HAL_GPIO_ReadPin(MRF_INT_GPIO_Port, MRF_INT_Pin));
}

static void MRF_BitBangDelay(void)
{
	for (volatile uint32_t n = 0U; n < 100U; n++)
	{
	}
}

static uint8_t MRF_BitBangTransfer(uint8_t tx)
{
	uint8_t rx = 0U;

	for (uint8_t mask = 0x80U; mask != 0U; mask >>= 1U)
	{
		HAL_GPIO_WritePin(SPI2_MOSI_MRF_GPIO_Port, SPI2_MOSI_MRF_Pin,
						  (tx & mask) ? GPIO_PIN_SET : GPIO_PIN_RESET);
		MRF_BitBangDelay();
		HAL_GPIO_WritePin(SPI2_SCK_MRF_GPIO_Port, SPI2_SCK_MRF_Pin, GPIO_PIN_SET);
		MRF_BitBangDelay();
		if (HAL_GPIO_ReadPin(SPI2_MISO_MRF_GPIO_Port, SPI2_MISO_MRF_Pin) == GPIO_PIN_SET)
		{
			rx |= mask;
		}
		HAL_GPIO_WritePin(SPI2_SCK_MRF_GPIO_Port, SPI2_SCK_MRF_Pin, GPIO_PIN_RESET);
		MRF_BitBangDelay();
	}

	return rx;
}

static void MRF_BitBangWriteLong(uint16_t address, uint8_t data)
{
	uint8_t ahigh = 0x80U | ((address >> 3) & 0x7FU);
	uint8_t alow = (uint8_t)(((address & 0x07U) << 5) | 0x10U);

	HAL_GPIO_WritePin(SPI2_CS_MRF_GPIO_Port, SPI2_CS_MRF_Pin, GPIO_PIN_RESET);
	MRF_BitBangDelay();
	(void)MRF_BitBangTransfer(ahigh);
	(void)MRF_BitBangTransfer(alow);
	(void)MRF_BitBangTransfer(data);
	MRF_BitBangDelay();
	HAL_GPIO_WritePin(SPI2_CS_MRF_GPIO_Port, SPI2_CS_MRF_Pin, GPIO_PIN_SET);
}

static void MRF_BitBangPrimeRfRegs(void)
{
	GPIO_InitTypeDef GPIO_InitStruct = {0};

	GPIO_InitStruct.Pin = SPI2_MOSI_MRF_Pin|SPI2_CS_MRF_Pin;
	GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
	GPIO_InitStruct.Pull = GPIO_NOPULL;
	GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
	HAL_GPIO_Init(GPIOC, &GPIO_InitStruct);

	GPIO_InitStruct.Pin = SPI2_SCK_MRF_Pin;
	GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
	GPIO_InitStruct.Pull = GPIO_NOPULL;
	GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
	HAL_GPIO_Init(SPI2_SCK_MRF_GPIO_Port, &GPIO_InitStruct);

	GPIO_InitStruct.Pin = SPI2_MISO_MRF_Pin;
	GPIO_InitStruct.Mode = GPIO_MODE_INPUT;
	GPIO_InitStruct.Pull = GPIO_NOPULL;
	HAL_GPIO_Init(SPI2_MISO_MRF_GPIO_Port, &GPIO_InitStruct);

	HAL_GPIO_WritePin(SPI2_CS_MRF_GPIO_Port, SPI2_CS_MRF_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(SPI2_SCK_MRF_GPIO_Port, SPI2_SCK_MRF_Pin, GPIO_PIN_RESET);
	HAL_GPIO_WritePin(SPI2_MOSI_MRF_GPIO_Port, SPI2_MOSI_MRF_Pin, GPIO_PIN_RESET);
	HAL_Delay(5);

	MRF_BitBangWriteLong(MRF_RFCON0, 0x93U);
	MRF_BitBangWriteLong(MRF_RFCON1, 0x02U);
	MRF_BitBangWriteLong(MRF_RFCON2, 0x80U);
	MRF_BitBangWriteLong(MRF_RFCON6, 0x90U);
	MRF_BitBangWriteLong(MRF_RFCON7, 0x80U);
	MRF_BitBangWriteLong(MRF_RFCON8, 0x10U);
}


/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{

  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();


  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */
  HAL_Delay(500); //this helps with the problem of self moving motors. This lets the physical power spike dissipate completely before the STM32 touches the communication lines.
  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */

  MX_GPIO_Init();
  MX_DMA_Init();
  MX_USART2_UART_Init();
  HAL_GPIO_WritePin(SPI2_CS_MRF_GPIO_Port, SPI2_CS_MRF_Pin, GPIO_PIN_SET);
  HAL_GPIO_WritePin(MRF_RESET_GPIO_Port, MRF_RESET_Pin, GPIO_PIN_RESET);
  HAL_Delay(100);
  HAL_GPIO_WritePin(MRF_RESET_GPIO_Port, MRF_RESET_Pin, GPIO_PIN_SET);
  HAL_Delay(200);
  MRF_BitBangPrimeRfRegs();
  MX_SPI2_Init();
  MRF_InitSPI(&hspi2, &huart2, 0);
  MRF_InitGPIO(SPI2_CS_MRF_GPIO_Port, SPI2_CS_MRF_Pin,
		  	   MRF_RESET_GPIO_Port, MRF_RESET_Pin,
			   MRF_INT_GPIO_Port, MRF_INT_Pin);
  Mrf24j_init_no_reset(20u);
  MRF_PrintInitDiag();
  MX_CAN1_Init();
  //MX_TIM6_Init();
  MX_I2C3_Init();
  MX_UART5_Init();
  /* USER CODE BEGIN 2 */

  // UART
  //__HAL_UART_ENABLE_IT(&huart2, UART_IT_TC);
  __HAL_UART_ENABLE_IT(&huart2, UART_IT_RXNE); //
 //__HAL_UART_ENABLE_IT(&huart5, UART_IT_TC);
 // __HAL_UART_ENABLE_IT(&huart5, UART_IT_RXNE);
 // HAL_UART_Receive_IT(&huart5, RxData, 10);

  /*MPU 6050*/
  //MPU6050_Init(&hi2c3);

  BNO_Init(&hi2c3, &huart2, INT_Pin);

  /*  Start CAN */
  HAL_CAN_Start(&hcan1);
  HAL_CAN_ActivateNotification(&hcan1, CAN_IT_RX_FIFO0_MSG_PENDING);
  HAL_NVIC_SetPriority(CAN1_RX0_IRQn, 1, 1);
  HAL_NVIC_EnableIRQ(CAN1_RX0_IRQn);

  //HAL_Delay(10);

  /****************     CyberGear Settings   ****************************/
	  int motor_mode_flag=1;

	  if(motor_mode_flag==0) //Position mode
	  {
		  // ADD code
		  /* Clear fault */
		    clearMotorFault(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_1_ID);
		    HAL_Delay(10);
		    clearMotorFault(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_2_ID);
		    HAL_Delay(10);
		    clearMotorFault(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_3_ID);
		    HAL_Delay(10);
		    clearMotorFault(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_4_ID);
		    HAL_Delay(10);

		    /*Position Mode - 0x01*/
		    MotorControlMode(/*Mode*/0x01,/*hostID=*/ CYBER_HOST_ID, /*motorID=*/ CYBER_MOTOR_1_ID);
		    HAL_Delay(50);
		    MotorControlMode(0x01,CYBER_HOST_ID,CYBER_MOTOR_2_ID);
		    HAL_Delay(50);
		    MotorControlMode(0x01,CYBER_HOST_ID,CYBER_MOTOR_3_ID);
		    HAL_Delay(50);
		    MotorControlMode(0x01,CYBER_HOST_ID,CYBER_MOTOR_4_ID);
		    HAL_Delay(50);

		    PositionSpeedLimit(speedLimit, CYBER_HOST_ID, CYBER_MOTOR_1_ID);
		    HAL_Delay(50);
		    PositionSpeedLimit(speedLimit, CYBER_HOST_ID, CYBER_MOTOR_2_ID);
		    HAL_Delay(50);
		    PositionSpeedLimit(speedLimit, CYBER_HOST_ID, CYBER_MOTOR_3_ID);
		    HAL_Delay(50);
		    PositionSpeedLimit(speedLimit, CYBER_HOST_ID, CYBER_MOTOR_4_ID);
		    HAL_Delay(50);

		    /* Mechanical ZERO */
		    setMechanicalZero(CYBER_HOST_ID,CYBER_MOTOR_1_ID);
		    HAL_Delay(50);
		    setMechanicalZero(CYBER_HOST_ID,CYBER_MOTOR_2_ID);
		    HAL_Delay(50);
		    setMechanicalZero(CYBER_HOST_ID,CYBER_MOTOR_3_ID);
		    HAL_Delay(50);
		    setMechanicalZero(CYBER_HOST_ID,CYBER_MOTOR_4_ID);
		    HAL_Delay(100);

		    /* Enable motor */
		    motorEnable(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_1_ID);
		    HAL_Delay(50);
		    motorEnable(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_2_ID);
		    HAL_Delay(50);
		    motorEnable(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_3_ID);
		    HAL_Delay(50);
		    motorEnable(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_4_ID);
		    HAL_Delay(50);

	  }
	  else //MIT mode
	  {
		  /* Clear fault */
		    clearMotorFault(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_1_ID);
		    HAL_Delay(10);
		    clearMotorFault(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_2_ID);
		    HAL_Delay(10);
		    clearMotorFault(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_3_ID);
		    HAL_Delay(10);
		    clearMotorFault(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_4_ID);
		    HAL_Delay(10);

		    /* Mechanical ZERO */
		    setMechanicalZero(CYBER_HOST_ID,CYBER_MOTOR_1_ID);
		    HAL_Delay(50);
		    setMechanicalZero(CYBER_HOST_ID,CYBER_MOTOR_2_ID);
		    HAL_Delay(50);
		    setMechanicalZero(CYBER_HOST_ID,CYBER_MOTOR_3_ID);
		    HAL_Delay(50);
		    setMechanicalZero(CYBER_HOST_ID,CYBER_MOTOR_4_ID);
		    HAL_Delay(100);


		    /* MIT Mode- 0x00 */
		    MotorControlMode(/*Mode*/0x00,/*hostID=*/ CYBER_HOST_ID, /*motorID=*/ CYBER_MOTOR_1_ID);
		    HAL_Delay(50);
		    MotorControlMode(0x00,CYBER_HOST_ID,CYBER_MOTOR_2_ID);
		    HAL_Delay(50);
		    MotorControlMode(0x00,CYBER_HOST_ID,CYBER_MOTOR_3_ID);
		    HAL_Delay(50);
		    MotorControlMode(0x00,CYBER_HOST_ID,CYBER_MOTOR_4_ID);
		    HAL_Delay(50);

		    /* Enable motor */
		    motorEnable(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_1_ID);
		    HAL_Delay(50);
		    motorEnable(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_2_ID);
		    HAL_Delay(50);
		    motorEnable(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_3_ID);
		    HAL_Delay(50);
		    motorEnable(/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_4_ID);
		    HAL_Delay(50);

		  // ADD code

	  }

   /* LEVEL CyberGear MOtors*/

  /****************     END of CyberGear Settings   *******************/

  /****************   DDSM115 Motor settings  ************************/

//int IDsetFlag = 1;
//while(IDsetFlag ==1){
//	ChangeMotorID(0x30);
//	GetMotorID();
//	HAL_Delay(100);
//	if (GetMotorID() == DDSM_LEFT_ID) break;
//}

	  (void)UART5_StartReceiveDmaIfReady();

	  CurrentMode(DDSM_LEFT_ID);
	  HAL_Delay(5);
	  sendCurrentCommand(DDSM_LEFT_ID, 0.00f); // LEVI
	  HAL_Delay(5);
	  CurrentMode(DDSM_RIGHT_ID);
	  HAL_Delay(5);
	  sendCurrentCommand(DDSM_RIGHT_ID, 0.00f); // DESNI
	  HAL_Delay(5);

  /****************  END of  DDSM115 Motor settings  *****************/




   /*Start timer interrupt for MPU6050 and Controller execution*/
	 //HAL_TIM_Base_Start_IT(&htim6);
	 //HAL_UART_Transmit(&huart2,"\n\rRun :",7,1000);
  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */


	 float xd_L, x_L, vd_L, v_L, fid_L, fi_L, wd_L, w_L;
	 float xd_R, x_R, vd_R, v_R, fid_R, fi_R, wd_R, w_R;
	 float u_L, u_R;
	 float k = 0.2;
	 float k1 = 23;   // 23.8188;					// linearna hitrost robota
	 float k2 = 16;    // 16.7196;		// 5		// linearni pospesek robota
	 float k3 = 77;   // 77.3774;					// naklon robota
	 float k4 = 3;    // 4.5646;		// 3		// kotna hitrost naklona robota
	 float senzor_offset = 0.03f;       // 0.03

	 //Preverjanje ce je v zraku
	 // Varnostne spremenljivke za detekcijo dviga (Pick-up detection)
	 uint32_t air_detect_start_time = 0;
	 uint8_t robot_is_in_air = 0;

	 float CG1_air_position = 0;
	 float CG2_air_position = 0;
	 float CG3_air_position = 0;
	 float CG4_air_position = 0;

	 // Nastavi mejo glede na tvoje enote (npr. če je v_L v m/s, je 2.0 m/s verjetno prehitro za normalno stanje)
	 const float MAX_SAFE_SPEED = 1.47f;
	 // Čas (v milisekundah), ki mora miniti nad to hitrostjo, da prožimo napako
	 const uint32_t AIR_TIME_THRESHOLD_MS = 250;
	 const uint32_t FALLING_TIME_THRESHOLD_MS = 5;


  while (1)
  {
	  MRF_ServiceCommand();

	  if (button_step_requested)
	  {
		  button_step_requested = 0;


		  HAL_Delay(2000);

		  BNO_App(robot_data);
		  RobotData_Update();
		  StateData_Update();

		  mrf_active_xd_L = state_data_L[0];
		  mrf_active_xd_R = state_data_R[0];
		  mrf_final_xd_L = mrf_active_xd_L;
		  mrf_final_xd_R = mrf_active_xd_R;
		  mrf_ramp_active = 0U;
		  mrf_ramp_last_tick = HAL_GetTick();
		  mrf_regulation_active = 1U;
		  printf("MRF regulation started\r\n");

		  i = 0;

		  while (1)
		  {
			  /* robot_data is only for monitoring/export; it does not change control logic.
			   * BNO_App updates BNO fields when the BNO086 interrupt has new data.
			   * RobotData_Update copies the latest DDSM and CyberGear feedback values.
			   * StateData_Update creates x = [position, speed, angle, angular speed].
			   */
			  BNO_App(robot_data);
			  RobotData_Update();
			  StateData_Update();
			  MRF_ServiceCommand();
//			  printf("X: %6.3f | Y: %6.3f | Z: %6.3f\r\n",
//			          robot_data[30],
//			          robot_data[31],
//			          robot_data[32]);

			  uint32_t now = HAL_GetTick();
			  if ((now - debug_last_tick) >= DEBUG_PRINT_INTERVAL_MS)
			  {
				  debug_last_tick = now;
				  MRF_UpdatePositionRamp(now);

//				  printf("x_L: %7.3f\r\n"
//						 "v_L: %6.3f\r\n"
//						 "theta_L: %6.2f\r\n"
//						 "omega_L: %6.3f\r\n"
//						 "U_L = : %6.2f\r\n\n\n\n",
//						 state_data_L[0],  // Robot position x [m]
//						 state_data_L[1],  // Robot speed v [m/s]
//						 state_data_L[2],  // Robot angle theta [rad]
//						 state_data_L[3],  // Robot angular speed omega [rad/s]
//						 u_L);
//				  printf("Kot1: %6.3f\r\n"
//						 "Kot2: %6.2f\r\n"
//						 "Kot3: %6.3f\r\n"
//						 "Kot4 = : %6.2f\r\n\n\n\n",
//						 robot_data[4],  // Robot position x [m]
//						 robot_data[5],  // Robot speed v [m/s]
//						 robot_data[6],  // Robot angle theta [rad]
//						 robot_data[7]  // Robot angular speed omega [rad/s]
//						 );

				  //////////    REGULACIJA     /////////

				  xd_L = mrf_active_xd_L;
				  x_L = state_data_L[0];
				  vd_L = 0;
				  v_L = state_data_L[1];
				  fid_L = 0;
				  fi_L = state_data_L[2] + senzor_offset;
				  wd_L = 0;
				  w_L = state_data_L[3];

				  xd_R = mrf_active_xd_R;
				  x_R = state_data_R[0];
				  vd_R = 0;
				  v_R = state_data_R[1];
				  fid_R = 0;
				  fi_R = state_data_R[2] + senzor_offset;
				  wd_R = 0;
				  w_R = state_data_R[3];



				  if (fabs(v_L) > MAX_SAFE_SPEED || fabs(v_R) > MAX_SAFE_SPEED)
				  {
					  if (air_detect_start_time == 0)
					  {
						  air_detect_start_time = now; // Začni meriti čas
					  }
					  else if ((now - air_detect_start_time) > AIR_TIME_THRESHOLD_MS)
					  {
						  robot_is_in_air = 1; // Robot je defenitivno v zraku

						  CG1_air_position = robot_data[4]; //shrani pozicije, da lahko zaznamo kdaj je spet na tleh
						  CG2_air_position = robot_data[5];
						  CG3_air_position = robot_data[6];
						  CG4_air_position = robot_data[7];

					  }
				  }
				  else
				  {
					  air_detect_start_time = 0; // Resetiraj časovnik, če hitrost pade nazaj na varno mejo
				  }

				  //////////    UPRAVLJANJE MOTORJEV    /////////
				  if (robot_is_in_air)
				  {
					  // Varnostni izklop DDSM motorjev

					  u_L = 0.0f;
					  u_R = 0.0f;


					  // Opcijsko izpiši opozorilo
					  printf("NAPAKA: Robot je v zraku! Izklapljam pogon.\r\n");

					  if (robot_data[4] + 0.1 < CG1_air_position && robot_data[5] + 0.1 > CG2_air_position && robot_data[6] + 0.1 > CG3_air_position && robot_data[7] + 0.1 < CG4_air_position)
					  {
						  robot_is_in_air = 0;

						  revolutionx01 = 0;
						  //delatx01 = DDSangle01 * Rwheel * PI / 180.0f;
						  zerocrossx01up = 0;
						  zerocrossx01down = 0;

						  revolutionx30 = 0;
						  //delatx30 = DDSangle01 * Rwheel * PI / 180.0f;
						  zerocrossx30up = 0;
						  zerocrossx30down = 0;


						  distnacex01 = 0.0f;
						  distnacex30 = 0.0f;
						  x_L = 0;
						  x_R = 0;
						  v_L = 0;
						  v_R = 0;



						  u_L = (k1*(xd_L-x_L) + k2*(vd_L-v_L) + k3*(fid_L - fi_L) + k4*(wd_L - w_L))*k;
						  u_R = (k1*(xd_R-x_R) + k2*(vd_R-v_R) + k3*(fid_R - fi_R) + k4*(wd_R - w_R))*k;

					  }
				  }
				  else
				  {
					  // Normalno delovanje

					  //////////    DDSM     /////////
					  u_L = (k1*(xd_L-x_L) + k2*(vd_L-v_L) + k3*(fid_L - fi_L) + k4*(wd_L - w_L))*k;
					  u_R = (k1*(xd_R-x_R) + k2*(vd_R-v_R) + k3*(fid_R - fi_R) + k4*(wd_R - w_R))*k;

				  }
				  ddsm_left_current_cmd = -(u_R/2)/0.75;		// upostevanje konstante motorja
				  ddsm_right_current_cmd = (u_L/2)/0.75;
			  }
			  /*
			 * Simple up/down MIT motion.
			 * Speed is controlled by target_step.
			 */

			  //////////    CYBERGEAR     /////////
			  //Uporabli bomo kaskadni pristop. Enačba za DDSM bo ostala enaka.
			  //Dodali bomo še posebno regulacijo za Cybergear kot vzmet dušilec.

			(void)UART5_StartReceiveDmaIfReady();


			  if (button_step_requested)
			  {
			      button_step_requested = 0;

			      if (i == 0)
			      {

			    	  //MIT
					  torque = 0;
					  target_step = MRF_UP_TARGET_STEP;
					  velocity = 0;
					  i = 1;
			      }
			      else
			      {
			          desired_angle  += 1;
			          desired_angle2 += -1;
			          desired_angle3 += -1;
			          desired_angle4 += 1;
			          ddsm_left_current_cmd = 0.0f;
			          ddsm_right_current_cmd = 0.0f;

			          //MIT
			    	  torque = 0;
			    	  target_step = MRF_DOWN_TARGET_STEP;
			    	  velocity = 0;

			    	  i = 0;
			      }
			  }

			  DDSM_Service();

			  CyberGear_Service();

		  }
		  button_step_requested = 0;
	  }
	}
}

/************************************
   END MAIN
/**************************************



/*DDS motor Callback*/
void HAL_UART_RxCpltCallback(UART_HandleTypeDef *huart)
 {


	if (huart->Instance == UART5)
	{

			 memcpy( Allocbuffer485,buffer485,10);
			 (void)UART5_StartReceiveDmaIfReady();

			 /*Process data*/
			 if(Allocbuffer485[0]== DDSM_LEFT_ID) //Motor Address
			 {
				  DDSangle01    = (float)(Allocbuffer485[6]<<8 | Allocbuffer485[7]) * 360.0f/32767.0f  ;    //Angle   [0-32767] -> [0 360]   address= 0x10 over RS485
				  /*MAX - MIN filter*/
				  if(abs(DDSangle01-DDSangle01_old)<0.2) //filter noise 0.2deg
				  {
					  DDSangle01=DDSangle01_old;
				  }
				  DDSangle01_old=DDSangle01;

				  //DDSangle01    = (float)(buffer485[6]<<8 | buffer485[7]) * 360.0f/32767.0f + (delatx01) ;    //Angle   [0-32767] -> [0 360]   address= 0x10 over RS485
				  DDSrpm01            = (float)((int16_t)(Allocbuffer485[4]<<8 | Allocbuffer485[5])) * 1.0f;              //RPM [-330RPM  330RPM]
				  DDScurrent01        = (float)((int16_t)(Allocbuffer485[2]<<8 | Allocbuffer485[3])) * 8.0f/32767.0f;              //Current Torque   [-32767   32767] -> [-8A   8A]
				  DDSvelocityRadial01 = DDSrpm01 * 0.10472 ;    //Velocity rad/s
				  Distnacex01();
				  DDSM_HandleFeedback(DDSM_LEFT_ID);

			 }
			 else if(Allocbuffer485[0]== DDSM_RIGHT_ID)
			 {
				   DDSangle30    = (float)(Allocbuffer485[6]<<8 | Allocbuffer485[7]) * 360.0f/32767.0f ;    //Angle   [0-32767] -> [0 360]   address= 0x30 over RS485

				   /*MAX - MIN filter*/
					   if(abs(DDSangle30-DDSangle30_old)<0.2) //filter noise 0.2deg
					   {
						  DDSangle30=DDSangle30_old;
					   }
					  DDSangle30_old=DDSangle30;

				   //DDSangle30    = (float)(buffer485[6]<<8 | buffer485[7]) * 360.0f/32767.0f + (delatx30);    //Angle   [0-32767] -> [0 360]   address= 0x30 over RS485
				  DDSrpm30            = (float)((int16_t)(Allocbuffer485[4]<<8 | Allocbuffer485[5])) * 1.0f;        //RPM [-330RPM  330RPM]
				  DDScurrent30        = (float)((int16_t)(Allocbuffer485[2]<<8 | Allocbuffer485[3])) * 8.0f/32767.0f;       //Current Torque   [-32767   32767] -> [-8A   8A]
				  DDSvelocityRadial30 = DDSrpm30 * 0.10472 ;    //Velocity rad/s
				  Distnacex30();
				  DDSM_HandleFeedback(DDSM_RIGHT_ID);
			 }

	}/*END of Instance UART5*/

 }

void HAL_GPIO_EXTI_Callback(uint16_t GPIO_Pin)
{
    static uint32_t last_button_tick = 0U;

    if(GPIO_Pin == BlueButton_Pin)
    {
        uint32_t now = HAL_GetTick();

        if((now - last_button_tick) >= BUTTON_DEBOUNCE_MS)
        {
            last_button_tick = now;
            button_step_requested = 1;
        }
    }

    if(GPIO_Pin == MRF_INT_Pin)
    {
        mrf_exti_count++;
        mrf_irq_pending = 1U;
    }

    BNO_EXTI_Callback(GPIO_Pin);
}


/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Configure the main internal regulator output voltage
  */
  __HAL_RCC_PWR_CLK_ENABLE();
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE3);

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI;
  RCC_OscInitStruct.HSIState = RCC_HSI_ON;
  RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSI;
  RCC_OscInitStruct.PLL.PLLM = 16;
  RCC_OscInitStruct.PLL.PLLN = 336;
  RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV4;
  RCC_OscInitStruct.PLL.PLLQ = 2;
  RCC_OscInitStruct.PLL.PLLR = 2;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV2;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_2) != HAL_OK)
  {
    Error_Handler();
  }
}

/**
  * @brief CAN1 Initialization Function
  * @param None
  * @retval None
  */
static void MX_CAN1_Init(void)
{

  /* USER CODE BEGIN CAN1_Init 0 */

  /* USER CODE END CAN1_Init 0 */

  /* USER CODE BEGIN CAN1_Init 1 */

  /* USER CODE END CAN1_Init 1 */
  hcan1.Instance = CAN1;
  hcan1.Init.Prescaler = 2;
  hcan1.Init.Mode = CAN_MODE_NORMAL;
  hcan1.Init.SyncJumpWidth = CAN_SJW_1TQ;
  hcan1.Init.TimeSeg1 = CAN_BS1_12TQ;
  hcan1.Init.TimeSeg2 = CAN_BS2_8TQ;
  hcan1.Init.TimeTriggeredMode = DISABLE;
  hcan1.Init.AutoBusOff = DISABLE;
  hcan1.Init.AutoWakeUp = DISABLE;
  hcan1.Init.AutoRetransmission = DISABLE;
  hcan1.Init.ReceiveFifoLocked = DISABLE;
  hcan1.Init.TransmitFifoPriority = DISABLE;
  if (HAL_CAN_Init(&hcan1) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN CAN1_Init 2 */
  CAN_FilterTypeDef sFilterConfig;
  sFilterConfig.FilterBank = 0;
  sFilterConfig.FilterMode = CAN_FILTERMODE_IDMASK;
  sFilterConfig.FilterScale = CAN_FILTERSCALE_32BIT;
  sFilterConfig.FilterIdHigh = 0x0000;
  sFilterConfig.FilterIdLow  = 0x0000;
  sFilterConfig.FilterMaskIdHigh = 0x0000;
  sFilterConfig.FilterMaskIdLow  = 0x0000;
  sFilterConfig.FilterFIFOAssignment = CAN_FILTER_FIFO0;
  sFilterConfig.FilterActivation = ENABLE;
  sFilterConfig.SlaveStartFilterBank = 14;
  HAL_CAN_ConfigFilter(&hcan1, &sFilterConfig);
  /* USER CODE END CAN1_Init 2 */

}

/**
  * @brief I2C3 Initialization Function
  * @param None
  * @retval None
  */
static void MX_I2C3_Init(void)
{

  /* USER CODE BEGIN I2C3_Init 0 */

  /* USER CODE END I2C3_Init 0 */

  /* USER CODE BEGIN I2C3_Init 1 */

  /* USER CODE END I2C3_Init 1 */
  hi2c3.Instance = I2C3;
  hi2c3.Init.ClockSpeed = 400000;
  hi2c3.Init.DutyCycle = I2C_DUTYCYCLE_2;
  hi2c3.Init.OwnAddress1 = 0;
  hi2c3.Init.AddressingMode = I2C_ADDRESSINGMODE_7BIT;
  hi2c3.Init.DualAddressMode = I2C_DUALADDRESS_DISABLE;
  hi2c3.Init.OwnAddress2 = 0;
  hi2c3.Init.GeneralCallMode = I2C_GENERALCALL_DISABLE;
  hi2c3.Init.NoStretchMode = I2C_NOSTRETCH_DISABLE;
  if (HAL_I2C_Init(&hi2c3) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN I2C3_Init 2 */

  /* USER CODE END I2C3_Init 2 */

}

/**
  * @brief SPI2 Initialization Function
  * @param None
  * @retval None
  */
static void MX_SPI2_Init(void)
{

  /* USER CODE BEGIN SPI2_Init 0 */

  /* USER CODE END SPI2_Init 0 */

  /* USER CODE BEGIN SPI2_Init 1 */

  /* USER CODE END SPI2_Init 1 */
  hspi2.Instance = SPI2;
  hspi2.Init.Mode = SPI_MODE_MASTER;
  hspi2.Init.Direction = SPI_DIRECTION_2LINES;
  hspi2.Init.DataSize = SPI_DATASIZE_8BIT;
  hspi2.Init.CLKPolarity = SPI_POLARITY_LOW;
  hspi2.Init.CLKPhase = SPI_PHASE_1EDGE;
  hspi2.Init.NSS = SPI_NSS_SOFT;
  hspi2.Init.BaudRatePrescaler = SPI_BAUDRATEPRESCALER_256;
  hspi2.Init.FirstBit = SPI_FIRSTBIT_MSB;
  hspi2.Init.TIMode = SPI_TIMODE_DISABLE;
  hspi2.Init.CRCCalculation = SPI_CRCCALCULATION_DISABLE;
  hspi2.Init.CRCPolynomial = 10;
  if (HAL_SPI_Init(&hspi2) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN SPI2_Init 2 */

  /* USER CODE END SPI2_Init 2 */

}

/**
  * @brief TIM6 Initialization Function
  * @param None
  * @retval None
  */
static void MX_TIM6_Init(void)
{

  /* USER CODE BEGIN TIM6_Init 0 */

  /* USER CODE END TIM6_Init 0 */

  TIM_MasterConfigTypeDef sMasterConfig = {0};

  /* USER CODE BEGIN TIM6_Init 1 */

  /* USER CODE END TIM6_Init 1 */
  htim6.Instance = TIM6;
  htim6.Init.Prescaler = 84;
  htim6.Init.CounterMode = TIM_COUNTERMODE_UP;
  htim6.Init.Period = 10000;
  htim6.Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
  if (HAL_TIM_Base_Init(&htim6) != HAL_OK)
  {
    Error_Handler();
  }
  sMasterConfig.MasterOutputTrigger = TIM_TRGO_RESET;
  sMasterConfig.MasterSlaveMode = TIM_MASTERSLAVEMODE_DISABLE;
  if (HAL_TIMEx_MasterConfigSynchronization(&htim6, &sMasterConfig) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN TIM6_Init 2 */

  /* USER CODE END TIM6_Init 2 */

}

/**
  * @brief UART5 Initialization Function
  * @param None
  * @retval None
  */
static void MX_UART5_Init(void)
{

  /* USER CODE BEGIN UART5_Init 0 */

  /* USER CODE END UART5_Init 0 */

  /* USER CODE BEGIN UART5_Init 1 */

  /* USER CODE END UART5_Init 1 */
  huart5.Instance = UART5;
  huart5.Init.BaudRate = 115200;
  huart5.Init.WordLength = UART_WORDLENGTH_8B;
  huart5.Init.StopBits = UART_STOPBITS_1;
  huart5.Init.Parity = UART_PARITY_NONE;
  huart5.Init.Mode = UART_MODE_TX_RX;
  huart5.Init.HwFlowCtl = UART_HWCONTROL_NONE;
  huart5.Init.OverSampling = UART_OVERSAMPLING_16;
  if (HAL_UART_Init(&huart5) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN UART5_Init 2 */

  /* USER CODE END UART5_Init 2 */

}

/**
  * @brief USART2 Initialization Function
  * @param None
  * @retval None
  */
static void MX_USART2_UART_Init(void)
{

  /* USER CODE BEGIN USART2_Init 0 */

  /* USER CODE END USART2_Init 0 */

  /* USER CODE BEGIN USART2_Init 1 */

  /* USER CODE END USART2_Init 1 */
  huart2.Instance = USART2;
  huart2.Init.BaudRate = 115200;
  huart2.Init.WordLength = UART_WORDLENGTH_8B;
  huart2.Init.StopBits = UART_STOPBITS_1;
  huart2.Init.Parity = UART_PARITY_NONE;
  huart2.Init.Mode = UART_MODE_TX_RX;
  huart2.Init.HwFlowCtl = UART_HWCONTROL_NONE;
  huart2.Init.OverSampling = UART_OVERSAMPLING_16;
  if (HAL_UART_Init(&huart2) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN USART2_Init 2 */

  /* USER CODE END USART2_Init 2 */

}

/**
  * Enable DMA controller clock
  */
static void MX_DMA_Init(void)
{

  /* DMA controller clock enable */
  __HAL_RCC_DMA1_CLK_ENABLE();

  /* DMA interrupt init */
  /* DMA1_Stream0_IRQn interrupt configuration */
  HAL_NVIC_SetPriority(DMA1_Stream0_IRQn, 0, 1);
  HAL_NVIC_EnableIRQ(DMA1_Stream0_IRQn);

}

/**
  * @brief GPIO Initialization Function
  * @param None
  * @retval None
  */
static void MX_GPIO_Init(void)
{
  GPIO_InitTypeDef GPIO_InitStruct = {0};
  /* USER CODE BEGIN MX_GPIO_Init_1 */
  /* USER CODE END MX_GPIO_Init_1 */

  /* GPIO Ports Clock Enable */
  __HAL_RCC_GPIOC_CLK_ENABLE();
  __HAL_RCC_GPIOH_CLK_ENABLE();
  __HAL_RCC_GPIOA_CLK_ENABLE();
  __HAL_RCC_GPIOD_CLK_ENABLE();
  __HAL_RCC_GPIOB_CLK_ENABLE();

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(GPIOC, MRF_RESET_Pin|SPI2_CS_MRF_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(GPIOA, RS485_DIR_Pin|LD1_Pin|LD2_Pin|LD3_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(RST_GPIO_Port, RST_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin : BlueButton_Pin */
  GPIO_InitStruct.Pin = BlueButton_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_FALLING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(BlueButton_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : INT_Pin */
  GPIO_InitStruct.Pin = INT_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_FALLING;
  GPIO_InitStruct.Pull = GPIO_PULLUP;
  HAL_GPIO_Init(INT_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : MRF_INT_Pin */
  GPIO_InitStruct.Pin = MRF_INT_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_FALLING;
  GPIO_InitStruct.Pull = GPIO_PULLUP;
  HAL_GPIO_Init(MRF_INT_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pins : MRF_RESET_Pin SPI2_CS_MRF_Pin */
  GPIO_InitStruct.Pin = MRF_RESET_Pin|SPI2_CS_MRF_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOC, &GPIO_InitStruct);

  /*Configure GPIO pins : RS485_DIR_Pin LD1_Pin LD2_Pin LD3_Pin */
  GPIO_InitStruct.Pin = RS485_DIR_Pin|LD1_Pin|LD2_Pin|LD3_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOA, &GPIO_InitStruct);

  /*Configure GPIO pin : RST_Pin */
  GPIO_InitStruct.Pin = RST_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(RST_GPIO_Port, &GPIO_InitStruct);

  /* EXTI interrupt init*/
  HAL_NVIC_SetPriority(EXTI1_IRQn, 0, 1);
  HAL_NVIC_EnableIRQ(EXTI1_IRQn);

  HAL_NVIC_SetPriority(EXTI9_5_IRQn, 0, 1);
  HAL_NVIC_EnableIRQ(EXTI9_5_IRQn);

  HAL_NVIC_SetPriority(EXTI15_10_IRQn, 2, 0);
  HAL_NVIC_EnableIRQ(EXTI15_10_IRQn);

  /* USER CODE BEGIN MX_GPIO_Init_2 */
  /* USER CODE END MX_GPIO_Init_2 */
}

/* USER CODE BEGIN 4 */


//void HAL_UARTEx_RxEventCallback(UART_HandleTypeDef *huart, uint16_t Size)
//{
//	  HAL_UARTEx_ReceiveToIdle_IT(&huart4, RS485_RxBuffer, RS485_BUFFER_SIZE);
//}



// Function to send a position command over RS485 (using UART4)
// motorID: The motor ID to address (for example, 0xC8 or 0x01, depending on your configuration)
// angle_deg: The desired target angle in degrees




/*
 * Function used for MIT control.
 */


/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
