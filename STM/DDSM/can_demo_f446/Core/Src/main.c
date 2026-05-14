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
  * STUDENTJE EDIT
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include "stdio.h"
#include "string.h"
#include "mpu6050.h"
#include "DDSM115.h"
#include "MRF24J40.h"
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */

#ifndef XIAOMI_CYBERGEAR_DEFS_H
#define XIAOMI_CYBERGEAR_DEFS_H

#define CMD_POSITION                  0x1
#define CMD_REQUEST                   0x2
#define CMD_ENABLE                    0x3
#define CMD_STOP                      0x4
#define CMD_SET_MECH_POSITION_TO_ZERO 0x6
#define CMD_SET_CAN_ID                0x7
#define CMD_RAM_WRITE                0x12
#define CMD_GET_STATUS               0x15

#define ADDR_SPEED_KP              0x2014
#define ADDR_SPEED_KI              0x2015
#define ADDR_POSITION_KP           0x2016
#define ADDR_RUN_MODE              0x7005
#define ADDR_I_REF                 0x7006
#define ADDR_SPEED_REF             0x700A
#define ADDR_LIMIT_TORQUE          0x700B
#define ADDR_CURRENT_KP            0x7010
#define ADDR_CURRENT_KI            0x7011
#define ADDR_CURRENT_FILTER_GAIN   0x7014
#define ADDR_POSITION_REF          0x7016
#define ADDR_LIMIT_SPEED           0x7017
#define ADDR_LIMIT_CURRENT         0x7018

#define MODE_MOTION                  0x00
#define MODE_POSITION                0x01
#define MODE_SPEED                   0x02
#define MODE_CURRENT                 0x03

#define POS_MIN                   -12.5f
#define POS_MAX                    12.5f
#define V_MIN                     -30.0f
#define V_MAX                      30.0f
#define KP_MIN                      0.0f
#define KP_MAX                    500.0f
#define KI_MIN                      0.0f
#define KI_MAX                     10.0f
#define KD_MIN                      0.0f
#define KD_MAX                      5.0f
#define T_MIN                     -12.0f
#define T_MAX                      12.0f
#define I_MIN                     -27.0f
#define I_MAX                      27.0f
#define CYBERGEAR_4PI             12.566371f
#define CYBERGEAR_8PI             25.132741f
#define CURRENT_FILTER_GAIN_MIN     0.0f
#define CURRENT_FILTER_GAIN_MAX     1.0f

#define RET_CYBERGEAR_OK              0x00
#define RET_CYBERGEAR_MSG_NOT_AVAIL   0x01
#define RET_CYBERGEAR_INVALID_CAN_ID  0x02
#define RET_CYBERGEAR_INVALID_PACKET  0x03


#endif // !XIAOMI_CYBERGEAR_DEFS_H

/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
CAN_HandleTypeDef hcan1;

I2C_HandleTypeDef hi2c3;

SPI_HandleTypeDef hspi1;
SPI_HandleTypeDef hspi2;

UART_HandleTypeDef huart5;
UART_HandleTypeDef huart2;
DMA_HandleTypeDef hdma_uart5_rx;

/* USER CODE BEGIN PV */
CAN_TxHeaderTypeDef pTxHeader;
CAN_RxHeaderTypeDef pRxHeader;
CAN_FilterTypeDef sFilterConfig;
uint32_t pTxMailbox;

uint8_t received_data[8];

// UART
char RxUARTBuffer[256]="";
uint8_t RxUARTLength=0;
uint8_t RxSingleByte;
uint8_t onFlag = 0;
// uint8_t MOTOR_ID = 18;

extern rx_info_t rx_info;//MRF
extern tx_info_t tx_info;//MRF
volatile uint8_t button_step_requested = 0;
uint8_t motorIDs[4] = {17, 18, 19, 20};
volatile float motor_positions[4] = {0.0f, 0.0f, 0.0f, 0.0f};
volatile uint8_t motor_position_valid[4] = {0, 0, 0, 0};
/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_DMA_Init(void);
static void MX_USART2_UART_Init(void);
static void MX_CAN1_Init(void);
static void MX_I2C3_Init(void);
static void MX_UART5_Init(void);
static void MX_SPI1_Init(void);
static void MX_SPI2_Init(void);
/* USER CODE BEGIN PFP */
void HAL_CAN_RxFifo0MsgPendingCallback(CAN_HandleTypeDef *hcan);
void serialWrite(char data[]);
void serialProcessRxData();
void HAL_UART_RxCpltCallback(UART_HandleTypeDef *huart);
static HAL_StatusTypeDef CAN_AddTxMessageWait(CAN_TxHeaderTypeDef *txHeader, uint8_t *txData);
static void request_motor_status(uint8_t motor_ID);
static HAL_StatusTypeDef set_mechanical_zero(uint8_t motor_ID);
static int8_t motor_index_from_id(uint8_t motor_ID);

void stop_motor(uint8_t motor_ID);
void enable_motor(uint8_t motor_ID);
void init_motor(uint8_t mode,uint8_t motor_ID);
void set_run_mode (uint8_t mode,  uint8_t motor_ID);
void send_can_float_package(uint8_t motor_ID, uint16_t addr, float value, float min, float max);
void set_limit_torque (float value, uint8_t motor_ID,  float min, float max);
void set_limit_current (float value, uint8_t motor_ID,  float min, float max);
void set_limit_speed (float value, uint8_t motor_ID,  float min, float max);
void set_position_ref(float position, uint8_t motorID);
HAL_StatusTypeDef writeParameter(uint16_t paramIndex, const volatile void* paramValue,
                                 uint8_t hostID, uint8_t motorID);
HAL_StatusTypeDef getMotorDeviceID(uint8_t hostID, uint8_t motorID);


extern uint8_t RS485_RxBuffer[RS485_BUFFER_SIZE];
extern uint8_t AllocBuffer[RS485_BUFFER_SIZE];
/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */
MPU6050_t MPU6050;
char buffer[220];
extern double roll_kalman, pitch_kalman;
extern float VPCnum;
uint8_t read_intstat, read_intcon, read_RXMCR, read_BBREG1,read_RXFLUSH, frame_length;
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

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
 // MX_DMA_Init();
  MX_USART2_UART_Init();
  MX_CAN1_Init();
  MX_I2C3_Init();
  MX_UART5_Init();
  MX_SPI1_Init();
  MX_SPI2_Init();
  /* USER CODE BEGIN 2 */
  // UART
  __HAL_UART_ENABLE_IT(&huart2, UART_IT_TC);
  __HAL_UART_ENABLE_IT(&huart2, UART_IT_RXNE);
  HAL_UART_Receive_IT(&huart2, &RxSingleByte, 1);
  //char buffer[11] = "HelloWorld\n";
  //HAL_UART_Transmit_IT(&huart2,(uint8_t *) buffer, 11);
  HAL_UART_Transmit(&huart2, buffer, strlen(buffer),1000);

  /*  UART 4 - RS485 */
  HAL_UART_Receive_DMA(&huart5,RS485_RxBuffer,10); //DMA interrupt RECEIVE  //__HAL_UART_ENABLE_IT(&huart4, UART_IT_RXNE);


  // CAN
  // Konfiguracija zaglavja sporocila (posiljanje sporocila)


	// Konfiguracija CAN filtra (prejemanje sporocila)
  sFilterConfig.FilterBank = 0;
  sFilterConfig.FilterMode = CAN_FILTERMODE_IDMASK;
  sFilterConfig.FilterScale = CAN_FILTERSCALE_32BIT;
  sFilterConfig.FilterIdHigh = 0x00;//0x12<<5; // ID naprave je 0x11
  sFilterConfig.FilterIdLow = 0;
  sFilterConfig.FilterMaskIdHigh = 0;
  sFilterConfig.FilterMaskIdLow = 0;
  sFilterConfig.FilterFIFOAssignment = CAN_FILTER_FIFO0;
  sFilterConfig.FilterActivation = ENABLE;
  sFilterConfig.SlaveStartFilterBank = 14;

 HAL_CAN_ConfigFilter(&hcan1, &sFilterConfig);

	 // Zagon CAN 1 vmesnika
 HAL_CAN_Start(&hcan1);
	 // Zagon prekinitev na prejemni strani
 HAL_CAN_ActivateNotification(&hcan1, CAN_IT_RX_FIFO0_MSG_PENDING); // potrebujemo za prejemanje

	  //I2C INIT

	  // MPU6050_Init(&hi2c3);
	  //uint8_t motor_ID = 18;

	  //init_motor(MODE_POSITION,  motor_ID);

  for(int i = 0; i < 4; i++)
  {
	  getMotorDeviceID(0x01, motorIDs[i]);
	  HAL_Delay(50);
  }

//	getMotorDeviceID(0x01, motor_ID);
//
//	getMotorDeviceID(0x01, motor_ID2);
//	getMotorDeviceID(0x01, motor_ID3);
//	getMotorDeviceID(0x01, motor_ID4);

//  Mrf24j_reset();
//  Mrf24j_init();
//
//  Mrf24j_set_channel(20);
//  Mrf24j_rx_flush();
//  uint8_t rxflush = 15;
//  //  Mrf24j_rx_enable();
//  rxflush = Mrf24j_read_short(MRF_RXFLUSH);
//
//  char one[] = "Hello world";
//  int da_size = sizeof(one) / sizeof(char);
//
//		// This is _our_ address
//
//  Mrf24j_send16(0xAAff, one, da_size);
//
//
//  while(1){
//	  read_intstat = Mrf24j_read_short(MRF_INTSTAT);
//	  if(read_intstat & 0b00001000) //interrupt occured
//	  {
//		  frame_length = Mrf24j_read_long(0x300);
//		  //HAL_UART_Transmit (&huart3, "interrupt main\r\n", strlen("interrupt main\r\n"), HAL_MAX_DELAY);
//
//	  }
//  }

  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */

  float curr_pos = 0;
  uint32_t last_position_read_tick = HAL_GetTick();

  for(int i = 0; i < 4; i++)
  {
	  set_limit_speed(3.0f, motorIDs[i], 0.0f, V_MAX);
	  set_limit_current(1.0f, motorIDs[i], 0.0f, I_MAX);

	  stop_motor(motorIDs[i]);
	  HAL_Delay(50);

	  set_mechanical_zero(motorIDs[i]);
	  HAL_Delay(50);

	  motor_positions[i] = 0.0f;
	  motor_position_valid[i] = 0;

	  enable_motor(motorIDs[i]);
	  HAL_Delay(100);

	  uint8_t runMode = 1; 			// position mode
	  writeParameter(0x7005, &runMode, 0xFE, motorIDs[i]);

	  HAL_Delay(20);
	  //set_position_ref(0.0f, motorIDs[i]);
  }



//  set_limit_speed(10.0f, motor_ID, 0.0f, V_MAX); //set the maximum speed of the motor
//  set_limit_current(5.0,motor_ID, 0.0f, I_MAX); //current limit allows faster operation
//  stop_motor(motor_ID); //spremenil
//  HAL_Delay(100);		//spremenil
//  enable_motor(motor_ID); // turn on the
//  HAL_Delay(200); 		//spremenil
//  set_run_mode(1, motor_ID);
//  uint8_t runMode = 1; // 1 => position mode
//
//  writeParameter(0x7005, &runMode, /*hostID=*/0xFE, /*motorID=*/motor_ID);
//  set_position_ref(0.0,motor_ID);




  /*DDSM_115*/

  /* SET motor ID */
  int  IDset_fag = 1;
  uint8_t motor_ID = 01;
  uint8_t curr_ID = GetMotorID();
  char msg[30];
  int len;

  while(IDset_fag == 1) {
	  ChangeMotorID(motor_ID);
	  HAL_Delay(100);


	  if (GetMotorID() == motor_ID) {
		  len = sprintf(msg, "Motor set for ID: %02X\r\n", curr_ID);
		  break;
	  } else {
		  len = sprintf(msg, "Motor ID not set\n");
	  }
	  HAL_UART_Transmit(&huart2, (uint8_t *)msg, len, 1000);

  }




  if(IDset_fag == 1){
	  return IDset_fag;
  }


          /* GET motor ID */
//          	  GetMotorID();
//          	  MOTOR_ID= RS485_RxBuffer[0]; //GET ID VALUE -EXAMPLE
//
//          	  /* Motor mode*/
//          		  //CurrentMode(MOTOR_ID);
//          		   VelocityMode(MOTOR_ID);
//          	  	  //PositionMode(MOTOR_ID);
//          		 VPCnum=50;


  while (0) {/*
	  if(received_data[0] == 1){
		  HAL_GPIO_WritePin(GPIOA, LD1_Pin, GPIO_PIN_SET);
	  }else{
		  HAL_GPIO_WritePin(GPIOA, LD1_Pin, GPIO_PIN_RESET);
	  }*/
	   //init_motor(MODE_POSITION);

	  if(button_step_requested)
	  {
		  button_step_requested = 0;
		  if(curr_pos < 10)
		  {
			  curr_pos = curr_pos + 3.0f;
		  }
		  else
		  {
			  curr_pos = curr_pos - 3.0f;
		  }


		  set_position_ref(curr_pos, motorIDs[0]);
		  HAL_Delay(5);
		  set_position_ref(curr_pos, motorIDs[1]);
		  HAL_Delay(5);
		  set_position_ref(curr_pos, motorIDs[2]);
		  HAL_Delay(5);
		  set_position_ref(curr_pos, motorIDs[3]);
		  HAL_Delay(5);

		  sprintf(buffer, "B1 pressed: position %.1f rad\r\n", curr_pos);
		  HAL_UART_Transmit(&huart2, (uint8_t *)buffer, strlen(buffer), 1000);
	  }
	  else
	  {
		  HAL_Delay(10);
	  }

	  if((HAL_GetTick() - last_position_read_tick) >= 1000U)
	  {
		  last_position_read_tick += 1000U;

		  for(int i = 0; i < 4; i++)
		  {
			  request_motor_status(motorIDs[i]);
			  HAL_Delay(5);
		  }

		  float pos0 = motor_positions[0];
		  float pos1 = motor_positions[1];
		  float pos2 = motor_positions[2];
		  float pos3 = motor_positions[3];
		  uint8_t valid0 = motor_position_valid[0];
		  uint8_t valid1 = motor_position_valid[1];
		  uint8_t valid2 = motor_position_valid[2];
		  uint8_t valid3 = motor_position_valid[3];

		  sprintf(buffer, "pos: %u=%s%.2f %u=%s%.2f %u=%s%.2f %u=%s%.2f\r\n",
				  (unsigned)motorIDs[0], valid0 ? "" : "?", pos0,
				  (unsigned)motorIDs[1], valid1 ? "" : "?", pos1,
				  (unsigned)motorIDs[2], valid2 ? "" : "?", pos2,
				  (unsigned)motorIDs[3], valid3 ? "" : "?", pos3);
		  HAL_UART_Transmit(&huart2, (uint8_t *)buffer, strlen(buffer), 1000);
	  }



	  /*DDSM_115 Control*/
	      //sendCurrentCommand(MOTOR_ID, VPCnum);
	  	 //  sendVelocityCommand(MOTOR_ID, VPCnum);//-200 200
	  	  //sendPositionCommand(MOTOR_ID, VPCnum);
	//  HAL_Delay(5000);
	 // Mrf24j_send16(0xAAff, one, da_size);
	 // HAL_CAN_AddTxMessage(&hcan1, &pTxHeader, data, &pTxMailbox);

	  //HAL_Delay(1000);


    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
		 /*Read*/

		// MPU6050_Read_All(&hi2c3, &MPU6050);

		 //Calculate Pitch and Roll angle from ACC and call Kalman filter function
		// Acc_ptich_roll_Kalman(&MPU6050);

		//sprintf(buffer,"Roll: %.2f  Pitch: %.2f \n\r",roll_kalman, pitch_kalman);
		//sprintf(buffer,"Roll: %.2f  Pitch: %.2f \n\r",3, 1);
	//	HAL_UART_Transmit(&huart2, buffer, strlen(buffer), 1000);
  }
  /* USER CODE END 3 */
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
  hi2c3.Init.ClockSpeed = 100000;
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
  * @brief SPI1 Initialization Function
  * @param None
  * @retval None
  */
static void MX_SPI1_Init(void)
{

  /* USER CODE BEGIN SPI1_Init 0 */

  /* USER CODE END SPI1_Init 0 */

  /* USER CODE BEGIN SPI1_Init 1 */

  /* USER CODE END SPI1_Init 1 */
  /* SPI1 parameter configuration*/
  hspi1.Instance = SPI1;
  hspi1.Init.Mode = SPI_MODE_MASTER;
  hspi1.Init.Direction = SPI_DIRECTION_2LINES;
  hspi1.Init.DataSize = SPI_DATASIZE_8BIT;
  hspi1.Init.CLKPolarity = SPI_POLARITY_LOW;
  hspi1.Init.CLKPhase = SPI_PHASE_1EDGE;
  hspi1.Init.NSS = SPI_NSS_SOFT;
  hspi1.Init.BaudRatePrescaler = SPI_BAUDRATEPRESCALER_2;
  hspi1.Init.FirstBit = SPI_FIRSTBIT_MSB;
  hspi1.Init.TIMode = SPI_TIMODE_DISABLE;
  hspi1.Init.CRCCalculation = SPI_CRCCALCULATION_DISABLE;
  hspi1.Init.CRCPolynomial = 10;
  if (HAL_SPI_Init(&hspi1) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN SPI1_Init 2 */

  /* USER CODE END SPI1_Init 2 */

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
  /* SPI2 parameter configuration*/
  hspi2.Instance = SPI2;
  hspi2.Init.Mode = SPI_MODE_MASTER;
  hspi2.Init.Direction = SPI_DIRECTION_2LINES;
  hspi2.Init.DataSize = SPI_DATASIZE_8BIT;
  hspi2.Init.CLKPolarity = SPI_POLARITY_LOW;
  hspi2.Init.CLKPhase = SPI_PHASE_1EDGE;
  hspi2.Init.NSS = SPI_NSS_SOFT;
  hspi2.Init.BaudRatePrescaler = SPI_BAUDRATEPRESCALER_2;
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
  HAL_NVIC_SetPriority(DMA1_Stream0_IRQn, 0, 0);
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
  __HAL_RCC_GPIOB_CLK_ENABLE();
  __HAL_RCC_GPIOD_CLK_ENABLE();

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(GPIOC, MRF_RESET_Pin|SPI2_CS_MRF_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(RS485_DIR_GPIO_Port, RS485_DIR_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(SPI1_CS_HILSCHER_GPIO_Port, SPI1_CS_HILSCHER_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin : B1_Pin */
  GPIO_InitStruct.Pin = B1_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_FALLING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(B1_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pins : MRF_RESET_Pin SPI2_CS_MRF_Pin */
  GPIO_InitStruct.Pin = MRF_RESET_Pin|SPI2_CS_MRF_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOC, &GPIO_InitStruct);

  /*Configure GPIO pins : PA0 PA1 */
  GPIO_InitStruct.Pin = GPIO_PIN_0|GPIO_PIN_1;
  GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_VERY_HIGH;
  GPIO_InitStruct.Alternate = GPIO_AF8_UART4;
  HAL_GPIO_Init(GPIOA, &GPIO_InitStruct);

  /*Configure GPIO pin : RS485_DIR_Pin */
  GPIO_InitStruct.Pin = RS485_DIR_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(RS485_DIR_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : PA12 */
  GPIO_InitStruct.Pin = GPIO_PIN_12;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_RISING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(GPIOA, &GPIO_InitStruct);

  /*Configure GPIO pin : MRF_INT_Pin */
  GPIO_InitStruct.Pin = MRF_INT_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_RISING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(MRF_INT_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : SPI1_CS_HILSCHER_Pin */
  GPIO_InitStruct.Pin = SPI1_CS_HILSCHER_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(SPI1_CS_HILSCHER_GPIO_Port, &GPIO_InitStruct);

  /* EXTI interrupt init*/
  HAL_NVIC_SetPriority(EXTI9_5_IRQn, 0, 0);
  HAL_NVIC_EnableIRQ(EXTI9_5_IRQn);
  HAL_NVIC_SetPriority(EXTI15_10_IRQn, 2, 0);
  HAL_NVIC_EnableIRQ(EXTI15_10_IRQn);

/* USER CODE BEGIN MX_GPIO_Init_2 */
/* USER CODE END MX_GPIO_Init_2 */
}

/* USER CODE BEGIN 4 */

static HAL_StatusTypeDef CAN_AddTxMessageWait(CAN_TxHeaderTypeDef *txHeader, uint8_t *txData)
{
	uint32_t start = HAL_GetTick();

	while(HAL_CAN_GetTxMailboxesFreeLevel(&hcan1) == 0U)
	{
		if((HAL_GetTick() - start) > 10U)
		{
			return HAL_TIMEOUT;
		}
	}

	return HAL_CAN_AddTxMessage(&hcan1, txHeader, txData, &pTxMailbox);
}

static int8_t motor_index_from_id(uint8_t motor_ID)
{
	for(int8_t i = 0; i < 4; i++)
	{
		if(motorIDs[i] == motor_ID)
		{
			return i;
		}
	}

	return -1;
}

static void request_motor_status(uint8_t motor_ID)
{
	CAN_TxHeaderTypeDef txHeader;
	uint8_t txData[8] = {0};
	uint32_t extId = ((uint32_t)CMD_GET_STATUS << 24)
				   | ((uint32_t)0xFE << 8)
				   | (uint32_t)motor_ID;

	txHeader.ExtId = extId;
	txHeader.IDE = CAN_ID_EXT;
	txHeader.RTR = CAN_RTR_DATA;
	txHeader.DLC = 8;
	txHeader.TransmitGlobalTime = DISABLE;

	CAN_AddTxMessageWait(&txHeader, txData);
}

static HAL_StatusTypeDef set_mechanical_zero(uint8_t motor_ID)
{
	CAN_TxHeaderTypeDef txHeader;
	uint8_t txData[8] = {0};
	uint32_t extId = ((uint32_t)CMD_SET_MECH_POSITION_TO_ZERO << 24)
				   | ((uint32_t)0xFE << 8)
				   | (uint32_t)motor_ID;

	txHeader.ExtId = extId;
	txHeader.IDE = CAN_ID_EXT;
	txHeader.RTR = CAN_RTR_DATA;
	txHeader.DLC = 8;
	txHeader.TransmitGlobalTime = DISABLE;

	txData[0] = 1;

	return CAN_AddTxMessageWait(&txHeader, txData);
}

void HAL_CAN_RxFifo0MsgPendingCallback(CAN_HandleTypeDef *hcan){
	if(HAL_CAN_GetRxMessage(hcan, CAN_RX_FIFO0, &pRxHeader, received_data) != HAL_OK)
	{
		return;
	}

	uint32_t extId = pRxHeader.ExtId;
	uint8_t type = (uint8_t)((extId >> 24) & 0x1F);
	uint8_t motor_ID = (uint8_t)((extId >> 8) & 0xFF);
	if(type == CMD_REQUEST)
	{
		int8_t motorIndex = motor_index_from_id(motor_ID);

		if(motorIndex >= 0)
		{
			uint16_t rawPosition = ((uint16_t)received_data[0] << 8) | received_data[1];
			float position = ((float)rawPosition * CYBERGEAR_8PI / 65535.0f) - CYBERGEAR_4PI;

			motor_positions[motorIndex] = position;
			motor_position_valid[motorIndex] = 1;
		}
	}
}

void HAL_GPIO_EXTI_Callback(uint16_t GPIO_Pin)
{
	static uint32_t last_button_tick = 0;
	uint32_t now = HAL_GetTick();

	if(GPIO_Pin == B1_Pin && (now - last_button_tick) > 200U)
	{
		last_button_tick = now;
		button_step_requested = 1;
	}
}

void serialWrite(char data[]){
	HAL_UART_Transmit(&huart2, (uint8_t *) data, strlen(data), 10);
	HAL_UART_Transmit(&huart2,(uint8_t *)"\n",1,10);
}

void serialProcessRxData(){
	// Process data
	switch(RxUARTBuffer[0]){
		case '1':
			onFlag = 1;
			CAN_AddTxMessageWait(&pTxHeader, (uint8_t *)&onFlag);
			break;
		case '0':
			onFlag = 0;
			CAN_AddTxMessageWait(&pTxHeader, (uint8_t *)&onFlag);
			break;
		case 'H':
			serialWrite("Ok");
			break;
		default:
			break;
	}
	RxUARTLength = 0;
}

void HAL_UART_RxCpltCallback(UART_HandleTypeDef *huart)
{
	if(huart->Instance == USART2){
		if(RxSingleByte == '\n')
		{
			serialProcessRxData();
		}else
		{
			RxUARTBuffer[RxUARTLength] = RxSingleByte;
			RxUARTLength++;
		}
		HAL_UART_Receive_IT(&huart2, &RxSingleByte, 1);
	}

	//UART5 DMA
	HAL_UART_Receive_DMA(&huart5,RS485_RxBuffer,10);
	memcpy( AllocBuffer,RS485_RxBuffer,RS485_BUFFER_SIZE);
}

void stop_motor (uint8_t motor_ID) {
	// CAN
    // Konfiguracija zaglavja sporocila (posiljanje sporocila)
	 uint8_t cmd_id = CMD_STOP;
	 uint16_t option = 0x00;
	 uint8_t can_id = motor_ID;

	 pTxHeader.ExtId = cmd_id << 24 | option << 8 | can_id; //CAN ID, motor  ima extended ID (29-bitni)
	 pTxHeader.RTR = CAN_RTR_DATA;
	 pTxHeader.IDE = CAN_ID_EXT;
	 pTxHeader.DLC = 8;
	 pTxHeader.TransmitGlobalTime = DISABLE;

	uint8_t data[8] = {0x00};
	//HAL_CAN_Tx(_cybergear_can_id, CMD_STOP, _master_can_id, 8, data);
	CAN_AddTxMessageWait(&pTxHeader, data);
}

void enable_motor (uint8_t motor_ID) {
	// CAN
    // Konfiguracija zaglavja sporocila (posiljanje sporocila)
	uint8_t cmd_id = CMD_ENABLE;
	 uint16_t option = 0x00;
	 uint8_t can_id = motor_ID;

	 pTxHeader.ExtId = cmd_id << 24 | option << 8 | can_id; //CAN ID, motor  ima extended ID (29-bitni)
	 pTxHeader.RTR = CAN_RTR_DATA;
	 pTxHeader.IDE = CAN_ID_EXT;
	 pTxHeader.DLC = 8;//8
	 pTxHeader.TransmitGlobalTime = DISABLE;

	uint8_t data[8] = {0x00};

	CAN_AddTxMessageWait(&pTxHeader, data);
}

void set_limit_speed (float value, uint8_t motor_ID,  float min, float max) {
	// CAN
    // Konfiguracija zaglavja sporocila (posiljanje sporocila)
	uint8_t data[8] = {0x00};
	    data[0] = ADDR_LIMIT_SPEED & 0x00FF;
	    data[1] = ADDR_LIMIT_SPEED >> 8;

	    float val = (max < value) ? max : value;
	    val = (min > value) ? min : value;
	    memcpy(&data[4], &val, 4);

	uint8_t cmd_id = CMD_RAM_WRITE;
	 uint16_t option = 0x00;
	 uint8_t can_id = motor_ID;

	 pTxHeader.ExtId = cmd_id << 24 | option << 8 | can_id; //CAN ID, motor  ima extended ID (29-bitni)
	 pTxHeader.RTR = CAN_RTR_DATA;
	 pTxHeader.IDE = CAN_ID_EXT;
	 pTxHeader.DLC = 8;
	 pTxHeader.TransmitGlobalTime = DISABLE;


	//HAL_CAN_Tx(_cybergear_can_id, CMD_STOP, _master_can_id, 8, data);
	CAN_AddTxMessageWait(&pTxHeader, data);
}

void set_limit_current (float value, uint8_t motor_ID,  float min, float max) {
	// CAN
    // Konfiguracija zaglavja sporocila (posiljanje sporocila)
	uint8_t data[8] = {0x00};
	    data[0] = ADDR_LIMIT_CURRENT & 0x00FF;
	    data[1] = ADDR_LIMIT_CURRENT >> 8;

	    float val = (max < value) ? max : value;
	    val = (min > value) ? min : value;
	    memcpy(&data[4], &val, 4);

	 uint8_t cmd_id =  CMD_RAM_WRITE;
	 uint16_t option = 0x00;
	 uint8_t can_id = motor_ID;

	 pTxHeader.ExtId = cmd_id << 24 | option << 8 | can_id; //CAN ID, motor  ima extended ID (29-bitni)
	 pTxHeader.RTR = CAN_RTR_DATA;
	 pTxHeader.IDE = CAN_ID_EXT;
	 pTxHeader.DLC = 8;
	 pTxHeader.TransmitGlobalTime = DISABLE;


	//HAL_CAN_Tx(_cybergear_can_id, CMD_STOP, _master_can_id, 8, data);
	CAN_AddTxMessageWait(&pTxHeader, data);
}

void set_limit_torque (float value, uint8_t motor_ID,  float min, float max) {
	// CAN
    // Konfiguracija zaglavja sporocila (posiljanje sporocila)
	uint8_t data[8] = {0x00};
	    data[0] = ADDR_LIMIT_TORQUE & 0x00FF;
	    data[1] = ADDR_LIMIT_TORQUE >> 8;

	    float val = (max < value) ? max : value;
	    val = (min > value) ? min : value;
	    memcpy(&data[4], &val, 4);

	uint8_t cmd_id = CMD_RAM_WRITE;
	 uint16_t option = 0x00;
	 uint8_t can_id = motor_ID;

	 pTxHeader.ExtId = cmd_id << 24 | option << 8 | can_id; //CAN ID, motor  ima extended ID (29-bitni)
	 pTxHeader.RTR = CAN_RTR_DATA;
	 pTxHeader.IDE = CAN_ID_EXT;
	 pTxHeader.DLC = 8;
	 pTxHeader.TransmitGlobalTime = DISABLE;


	//HAL_CAN_Tx(_cybergear_can_id, CMD_STOP, _master_can_id, 8, data);
	CAN_AddTxMessageWait(&pTxHeader, data);
}

void send_can_float_package(uint8_t motor_ID, uint16_t addr, float value, float min, float max){
    uint8_t data[8] = {0x00};
    data[0] = addr & 0x00FF;
    data[1] = addr >> 8;

    float val = (max < value) ? max : value;
    val = (min > value) ? min : value;
    memcpy(&data[4], &val, 4);
    //_send_can_package(can_id, CMD_RAM_WRITE, _master_can_id, 8, data);

    // CAN
        // Konfiguracija zaglavja sporocila (posiljanje sporocila)
    	uint8_t cmd_id = CMD_RAM_WRITE;
    	 uint16_t option = 0x00;
    	 uint8_t can_id = motor_ID;

    	 pTxHeader.ExtId = cmd_id << 24 | option << 8 | can_id; //CAN ID, motor  ima extended ID (29-bitni)
    	 pTxHeader.RTR = CAN_RTR_DATA;
    	 pTxHeader.IDE = CAN_ID_EXT;
    	 pTxHeader.DLC = 8;
    	 pTxHeader.TransmitGlobalTime = DISABLE;

    CAN_AddTxMessageWait(&pTxHeader, data);
}

void set_run_mode (uint8_t mode,  uint8_t motor_ID) {

	    uint8_t data[8] = {0x00};
	    data[0] = ADDR_RUN_MODE & 0x00FF;
	    data[1] = ADDR_RUN_MODE >> 8;
	    data[4] = mode;

	    uint8_t cmd_id = CMD_RAM_WRITE;
	    uint16_t option = 0x00;
	    uint8_t can_id = motor_ID;

	    pTxHeader.ExtId = cmd_id << 24 | option << 8 | can_id; //CAN ID, motor  ima extended ID (29-bitni)
	    pTxHeader.RTR = CAN_RTR_DATA;
	    pTxHeader.IDE = CAN_ID_EXT;
	    pTxHeader.DLC = 8;
	    pTxHeader.TransmitGlobalTime = DISABLE;

	    CAN_AddTxMessageWait(&pTxHeader, data);
}

void init_motor (uint8_t mode,uint8_t motor_ID) {
	stop_motor(motor_ID);
	set_run_mode(mode,motor_ID);

}

void set_position_ref(float position,uint8_t motor_ID){
    send_can_float_package(motor_ID, ADDR_POSITION_REF, position, POS_MIN, POS_MAX);
}


HAL_StatusTypeDef writeParameter(uint16_t paramIndex, const volatile void* paramValue,
                                 uint8_t hostID, uint8_t motorID)

{
    CAN_TxHeaderTypeDef txHeader;
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

    return CAN_AddTxMessageWait(&txHeader, txData);
}

HAL_StatusTypeDef getMotorDeviceID(uint8_t hostID, uint8_t motorID)
{
    CAN_TxHeaderTypeDef txHeader;
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
    memset(txData, 0, 8);

    return CAN_AddTxMessageWait(&txHeader, txData);
}


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
