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
#include "DDSM115.h"
#include "mpu6050.h"
#include "CyberGear.h"
#include "MRF24J40.h"
#include <stdio.h>

#include "../../SH2Sensorhub/Inc/demo_app.h"
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */


/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
CAN_HandleTypeDef hcan1;

I2C_HandleTypeDef hi2c3;

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

float MOTangle11=0;      //Current angle    [0-65535] -> [-4pi  4pi]
float MOTvelocity11=0;   //Current velocity [0-65535] -> [-30rad/s  30rad/s]
float MOTtorque11=0;     //Current Torque   [0-65535] -> [-12Nm   12Nm]
float MOTtemp11=0;       //Temperature  10*Celsius
float MOTrpm11 =0;       //Current RPM

/*Drive Wheels*/
extern float DDSangle01;       //Current angle    [0-32767] -> [0 360]   address= 0x01 over RS485
extern float DDSrpm01;         //Current velocity [0-     ] -> [-330RPM  330RPM]
extern float DDScurrent01;     //Current Torque   [-32767   32767] -> [-8A   8A]
extern float DDSvelocityRadial01;    //Velocity

extern float DDSangle30;       //Current angle    [0-32767] -> [0 360]   address= 0x30 over RS485
extern float DDSrpm30;         //Current velocity [0-     ] -> [-330RPM  330RPM]
extern float DDScurrent30;     //Current Torque   [-32767   32767] -> [-8A   8A]
extern float DDSvelocityRadial30;    //Velocity

float desired_angle ;
float desired_angle2 ;
float desired_angle3 ;
float desired_angle4 ;

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


extern double roll_kalman, pitch_kalman;

// We'll do a simple function that sets position mode and moves the motor



/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_DMA_Init(void);
static void MX_USART2_UART_Init(void);
static void MX_CAN1_Init(void);
static void MX_TIM6_Init(void);
static void MX_I2C3_Init(void);
static void MX_UART5_Init(void);
/* USER CODE BEGIN PFP */
//void HAL_CAN_RxFifo0MsgPendingCallback(CAN_HandleTypeDef *hcan);
void serialWrite(char data[]);
void serialProcessRxData();
//void HAL_UART_RxCpltCallback(UART_HandleTypeDef *huart);
static uint16_t float_to_uint(float x, float x_min, float x_max);
void HAL_UART_RxCpltCallback(UART_HandleTypeDef *huart);
//void HAL_UARTEx_RxEventCallback(UART_HandleTypeDef *huart, uint16_t Size);

// Retarget printf to UART
int __io_putchar(int ch)
{
    HAL_UART_Transmit(&huart2, (uint8_t *)&ch, 1, HAL_MAX_DELAY);
    return ch;
}

void Debug_EXTI_Config(void)
{
    // Read EXTI_RTSR and EXTI_FTSR for line 1
    uint32_t rtsr = EXTI->RTSR;  // rising trigger
    uint32_t ftsr = EXTI->FTSR;  // falling trigger

    printf("EXTI RTSR bit1: %lu\r\n", (rtsr >> 1) & 1);  // should be 0
    printf("EXTI FTSR bit1: %lu\r\n", (ftsr >> 1) & 1);  // should be 1
}
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
  MX_DMA_Init();
  MX_USART2_UART_Init();
  MX_CAN1_Init();
  MX_TIM6_Init();
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


  /*  Start CAN */
  HAL_CAN_Start(&hcan1);
  HAL_CAN_ActivateNotification(&hcan1, CAN_IT_RX_FIFO0_MSG_PENDING);
  HAL_NVIC_SetPriority(CAN1_RX0_IRQn, 1, 1);
  HAL_NVIC_EnableIRQ(CAN1_RX0_IRQn);

  //HAL_Delay(10);

  /****************     CyberGear Settings   ****************************/
	  int motor_mode_flag=0;

	  if(motor_mode_flag==0) //MIT mode
	  {
		  // ADD code
		  /* Clear fault */
		    clearMotorFault(/*hostID=*/0x00, /*motorID=*/25);
		    HAL_Delay(10);
		    clearMotorFault(/*hostID=*/0x00, /*motorID=*/26);
		    HAL_Delay(10);
		    clearMotorFault(/*hostID=*/0x00, /*motorID=*/27);
		    HAL_Delay(10);
		    clearMotorFault(/*hostID=*/0x00, /*motorID=*/28);
		    HAL_Delay(10);

		    /*Position Mode - 0x01*/
		    MotorControlMode(/*Mode*/0x01,/*hostID=*/ 0x00, /*motorID=*/ 25);
		    HAL_Delay(50);
		    MotorControlMode(0x01,0x00,26);
		    HAL_Delay(50);
		    MotorControlMode(0x01,0x00,27);
		    HAL_Delay(50);
		    MotorControlMode(0x01,0x00,28);
		    HAL_Delay(50);

		    PositionSpeedLimit(speedLimit, 0x00, 25);
		    HAL_Delay(50);
		    PositionSpeedLimit(speedLimit, 0x00, 26);
		    HAL_Delay(50);
		    PositionSpeedLimit(speedLimit, 0x00, 27);
		    HAL_Delay(50);
		    PositionSpeedLimit(speedLimit, 0x00, 28);
		    HAL_Delay(50);

		    /* Mechanical ZERO */
		    setMechanicalZero(0x00,25);
		    HAL_Delay(50);
		    setMechanicalZero(0x00,26);
		    HAL_Delay(50);
		    setMechanicalZero(0x00,27);
		    HAL_Delay(50);
		    setMechanicalZero(0x00,28);
		    HAL_Delay(100);

		    /* Enable motor */
		    motorEnable(/*hostID=*/0x00, /*motorID=*/25);
		    HAL_Delay(50);
		    motorEnable(/*hostID=*/0x00, /*motorID=*/26);
		    HAL_Delay(50);
		    motorEnable(/*hostID=*/0x00, /*motorID=*/27);
		    HAL_Delay(50);
		    motorEnable(/*hostID=*/0x00, /*motorID=*/28);
		    HAL_Delay(50);

	  }
	  else //Position mode
	  {
		  /* Clear fault */
		    clearMotorFault(/*hostID=*/0x00, /*motorID=*/25);
		    HAL_Delay(10);
		    clearMotorFault(/*hostID=*/0x00, /*motorID=*/26);
		    HAL_Delay(10);
		    clearMotorFault(/*hostID=*/0x00, /*motorID=*/27);
		    HAL_Delay(10);
		    clearMotorFault(/*hostID=*/0x00, /*motorID=*/28);
		    HAL_Delay(10);

		    /* Mechanical ZERO */
		    setMechanicalZero(0x00,25);
		    HAL_Delay(50);
		    setMechanicalZero(0x00,26);
		    HAL_Delay(50);
		    setMechanicalZero(0x00,27);
		    HAL_Delay(50);
		    setMechanicalZero(0x00,28);
		    HAL_Delay(100);


		    /* MIT Mode- 0x00 */
		    MotorControlMode(/*Mode*/0x00,/*hostID=*/ 0x00, /*motorID=*/ 25);
		    HAL_Delay(50);
		    MotorControlMode(0x00,0x00,26);
		    HAL_Delay(50);
		    MotorControlMode(0x00,0x00,27);
		    HAL_Delay(50);
		    MotorControlMode(0x00,0x00,28);
		    HAL_Delay(50);

		    /* Enable motor */
		    motorEnable(/*hostID=*/0x00, /*motorID=*/25);
		    HAL_Delay(50);
		    motorEnable(/*hostID=*/0x00, /*motorID=*/26);
		    HAL_Delay(50);
		    motorEnable(/*hostID=*/0x00, /*motorID=*/27);
		    HAL_Delay(50);
		    motorEnable(/*hostID=*/0x00, /*motorID=*/28);
		    HAL_Delay(50);

		  // ADD code

	  }

   /* LEVEL CyberGear MOtors*/
  HAL_UART_Receive_DMA(&huart5,buffer485,10);

  /****************     END of CyberGear Settings   *******************/

  /****************   DDSM115 Motor settings  ************************/

//int IDsetFlag = 1;
//while(IDsetFlag ==1){
//	ChangeMotorID(0x30);
//	GetMotorID();
//	HAL_Delay(100);
//	if (GetMotorID() == 0x01) break;
//}
//CurrentMode(0x01);
//	  HAL_Delay(4);
//	  CurrentMode(0x30);
//	  HAL_Delay(4);
//
//	  /*Take initial values from DDSM motor encoders*/
//	  sendCurrentCommand(0x30, 0); //Send iq value
//	  HAL_Delay(4);
//	  sendCurrentCommand(0x01, 0);
//	  HAL_Delay(4);
//
//                    // ADD code for DDSM115
// sendCurrentCommand(0x30, /*iq*/ 0.4); //Send iq value
// 	  HAL_Delay(4);
// 	  sendCurrentCommand(0x01, /*iq*/ 0.2);
// 	  HAL_Delay(4);

  /****************  END of  DDSM115 Motor settings  *****************/




   /*Start timer interrupt for MPU6050 and Controller execution*/
	 //HAL_TIM_Base_Start_IT(&htim6);
	 HAL_UART_Transmit(&huart2,"\n\rRun :",7,1000);


	 //BNO settings
  Debug_EXTI_Config();
  BNO_Init(&hi2c3,&huart2,INT_Pin);

  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */

  float vector[4] = {0};

  while (1)
  {
	  sh2_service();
	  BNO_App(vector);

	// izpis 3. in 4. vrednosti iz vektorja

//    sprintf(buffer, "R:%6.2f | Gx:%6.3f\r\n",
//			   vector[2],
//			   vector[3]);
	//HAL_UART_Transmit(&huart2,(uint8_t *)buffer, strlen(buffer), 1000);
	 printf("R:%6.2f | Gx:%6.3f\r\n",vector[2],vector[3]);

	HAL_Delay(100);
  }
  while (0)
  {


    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */

	  /*Motor commands*/



	  HAL_UART_Transmit(&huart2,"\n\rRunD :",7,1000);
	  HAL_Delay(300);
	  //HAL_GPIO_TogglePin(LD2_GPIO_Port,LD2_Pin);
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
  __HAL_RCC_GPIOB_CLK_ENABLE();
  __HAL_RCC_GPIOD_CLK_ENABLE();

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

  /*Configure GPIO pins : MRF_RESET_Pin SPI2_CS_MRF_Pin */
  GPIO_InitStruct.Pin = MRF_RESET_Pin|SPI2_CS_MRF_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOC, &GPIO_InitStruct);

  /*Configure GPIO pin : INT_Pin */
  GPIO_InitStruct.Pin = INT_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_RISING;
  GPIO_InitStruct.Pull = GPIO_PULLUP;
  HAL_GPIO_Init(INT_GPIO_Port, &GPIO_InitStruct);

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
  HAL_NVIC_SetPriority(EXTI1_IRQn, 0, 0);
  HAL_NVIC_EnableIRQ(EXTI1_IRQn);

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
