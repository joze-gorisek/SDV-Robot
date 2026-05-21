/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file    stm32f4xx_it.c
  * @brief   Interrupt Service Routines.
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
#include "stm32f4xx_it.h"
/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include <stdbool.h>
#include "DDSM115.h"
#include "mpu6050.h"
#include "CyberGear.h"


extern void FrontAngle(void);
extern void BackAngle(void);
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN TD */

/* USER CODE END TD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */

/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
/* USER CODE BEGIN PV */
extern uint8_t MOTOR_ID;
extern double iq;
extern double iq01,iq30;
extern double iqmot;
int counter = 0;
int counter_startup = 0;
int counter_dist =0;


/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
/* USER CODE BEGIN PFP */

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */
	extern MPU6050_t MPU6050;
	extern char buffer[100];
	extern double roll_kalman, pitch_kalman;
	double kot = 0;

	int ControllerON       = 0;
	int ControllerONLaterl = 0;

    extern float current;
    extern float RPM;
	extern float angle;
	extern uint8_t RS485_RxBuffer[RS485_BUFFER_SIZE];
	extern uint8_t modeCmd[10];
	extern uint8_t position_mode[10];
	extern uint8_t command[10];
	extern uint8_t  ID_query[10];


	extern float front_angle;
	extern float back_angle;


/* USER CODE END 0 */

/* External variables --------------------------------------------------------*/
extern CAN_HandleTypeDef hcan1;
extern TIM_HandleTypeDef htim6;
extern DMA_HandleTypeDef hdma_uart5_rx;
extern UART_HandleTypeDef huart5;
extern UART_HandleTypeDef huart2;
/* USER CODE BEGIN EV */
#define RS485_BUFFER_SIZE 10
extern uint8_t RS485_RxBuffer[RS485_BUFFER_SIZE];
extern char HEX_Buffer[15];
extern uint8_t motor_status_flag;

extern MPU6050_t MPU6050;

char VCPBuffer[10];
float VPCnum=0.0f;
extern float angle;

uint8_t cnt =0; //Byte counter UART4
uint8_t cnta=0; //Byte counter USART2
int controller_flag=0;
int Cybercounter=0;
/* USER CODE END EV */

/******************************************************************************/
/*           Cortex-M4 Processor Interruption and Exception Handlers          */
/******************************************************************************/
/**
  * @brief This function handles Non maskable interrupt.
  */
void NMI_Handler(void)
{
  /* USER CODE BEGIN NonMaskableInt_IRQn 0 */

  /* USER CODE END NonMaskableInt_IRQn 0 */
  /* USER CODE BEGIN NonMaskableInt_IRQn 1 */
  while (1)
  {
  }
  /* USER CODE END NonMaskableInt_IRQn 1 */
}

/**
  * @brief This function handles Hard fault interrupt.
  */
void HardFault_Handler(void)
{
  /* USER CODE BEGIN HardFault_IRQn 0 */

  /* USER CODE END HardFault_IRQn 0 */
  while (1)
  {
    /* USER CODE BEGIN W1_HardFault_IRQn 0 */
    /* USER CODE END W1_HardFault_IRQn 0 */
  }
}

/**
  * @brief This function handles Memory management fault.
  */
void MemManage_Handler(void)
{
  /* USER CODE BEGIN MemoryManagement_IRQn 0 */

  /* USER CODE END MemoryManagement_IRQn 0 */
  while (1)
  {
    /* USER CODE BEGIN W1_MemoryManagement_IRQn 0 */
    /* USER CODE END W1_MemoryManagement_IRQn 0 */
  }
}

/**
  * @brief This function handles Pre-fetch fault, memory access fault.
  */
void BusFault_Handler(void)
{
  /* USER CODE BEGIN BusFault_IRQn 0 */

  /* USER CODE END BusFault_IRQn 0 */
  while (1)
  {
    /* USER CODE BEGIN W1_BusFault_IRQn 0 */
    /* USER CODE END W1_BusFault_IRQn 0 */
  }
}

/**
  * @brief This function handles Undefined instruction or illegal state.
  */
void UsageFault_Handler(void)
{
  /* USER CODE BEGIN UsageFault_IRQn 0 */

  /* USER CODE END UsageFault_IRQn 0 */
  while (1)
  {
    /* USER CODE BEGIN W1_UsageFault_IRQn 0 */
    /* USER CODE END W1_UsageFault_IRQn 0 */
  }
}

/**
  * @brief This function handles System service call via SWI instruction.
  */
void SVC_Handler(void)
{
  /* USER CODE BEGIN SVCall_IRQn 0 */

  /* USER CODE END SVCall_IRQn 0 */
  /* USER CODE BEGIN SVCall_IRQn 1 */

  /* USER CODE END SVCall_IRQn 1 */
}

/**
  * @brief This function handles Debug monitor.
  */
void DebugMon_Handler(void)
{
  /* USER CODE BEGIN DebugMonitor_IRQn 0 */

  /* USER CODE END DebugMonitor_IRQn 0 */
  /* USER CODE BEGIN DebugMonitor_IRQn 1 */

  /* USER CODE END DebugMonitor_IRQn 1 */
}

/**
  * @brief This function handles Pendable request for system service.
  */
void PendSV_Handler(void)
{
  /* USER CODE BEGIN PendSV_IRQn 0 */

  /* USER CODE END PendSV_IRQn 0 */
  /* USER CODE BEGIN PendSV_IRQn 1 */

  /* USER CODE END PendSV_IRQn 1 */
}

/**
  * @brief This function handles System tick timer.
  */
void SysTick_Handler(void)
{
  /* USER CODE BEGIN SysTick_IRQn 0 */

  /* USER CODE END SysTick_IRQn 0 */
  HAL_IncTick();
  /* USER CODE BEGIN SysTick_IRQn 1 */

  /* USER CODE END SysTick_IRQn 1 */
}

/******************************************************************************/
/* STM32F4xx Peripheral Interrupt Handlers                                    */
/* Add here the Interrupt Handlers for the used peripherals.                  */
/* For the available peripheral interrupt handler names,                      */
/* please refer to the startup file (startup_stm32f4xx.s).                    */
/******************************************************************************/

/**
  * @brief This function handles DMA1 stream0 global interrupt.
  */
void DMA1_Stream0_IRQHandler(void)
{
  /* USER CODE BEGIN DMA1_Stream0_IRQn 0 */

  /* USER CODE END DMA1_Stream0_IRQn 0 */
  HAL_DMA_IRQHandler(&hdma_uart5_rx);
  /* USER CODE BEGIN DMA1_Stream0_IRQn 1 */

  /* USER CODE END DMA1_Stream0_IRQn 1 */
}

/**
  * @brief This function handles CAN1 RX0 interrupt.
  */
void CAN1_RX0_IRQHandler(void)
{
  /* USER CODE BEGIN CAN1_RX0_IRQn 0 */

  /* USER CODE END CAN1_RX0_IRQn 0 */
  HAL_CAN_IRQHandler(&hcan1);
  /* USER CODE BEGIN CAN1_RX0_IRQn 1 */

  /* USER CODE END CAN1_RX0_IRQn 1 */
}

/**
  * @brief This function handles USART2 global interrupt.
  */
void USART2_IRQHandler(void)
{
  /* USER CODE BEGIN USART2_IRQn 0 */

  /* USER CODE END USART2_IRQn 0 */
  HAL_UART_IRQHandler(&huart2);
  /* USER CODE BEGIN USART2_IRQn 1 */


   if((USART2->SR & 0x0525U)!= RESET)
	{

	   VCPBuffer[cnta]=USART2->DR;
	   cnta++;

	   if( VCPBuffer[cnta-1]=='%') //Delimeter
		   {
		   VCPBuffer[cnta-1]='\0';

		      //HAL_UART_Transmit(&huart2,   VCPBuffer, strlen(  VCPBuffer),1000);

              //VPCnum=(float)atof( VCPBuffer); //Convert string to float

              //sprintf(HEX_Buffer,"data %.3f \n\r",VPCnum);
              //HAL_UART_Transmit(&huart2, HEX_Buffer, strlen(HEX_Buffer),1000);
              HAL_UART_Transmit(&huart2, VCPBuffer, strlen( VCPBuffer),10);

              //HAL_UART_Transmit_IT(&huart2, HEX_Buffer, strlen(HEX_Buffer));


		      for(int i=0; i<cnta; i++) //Clear Buffer
		      {
		    	  VCPBuffer[i]='\0';
		      }
		      cnta=0;

		      //__HAL_UART_ENABLE_IT(&huart2, UART_IT_RXNE);
		   }

	}
  /* USER CODE END USART2_IRQn 1 */
}

/**
  * @brief This function handles EXTI line 1 interrupt.
  */
void EXTI1_IRQHandler(void)
{
  /* USER CODE BEGIN EXTI1_IRQn 0 */

  /* USER CODE END EXTI1_IRQn 0 */
  HAL_GPIO_EXTI_IRQHandler(INT_Pin);
  /* USER CODE BEGIN EXTI1_IRQn 1 */

  /* USER CODE END EXTI1_IRQn 1 */
}

/**
  * @brief This function handles EXTI line[15:10] interrupts.
  */
void EXTI15_10_IRQHandler(void)
{
  /* USER CODE BEGIN EXTI15_10_IRQn 0 */

  /* USER CODE END EXTI15_10_IRQn 0 */
  HAL_GPIO_EXTI_IRQHandler(BlueButton_Pin);
  /* USER CODE BEGIN EXTI15_10_IRQn 1 */

  if(controller_flag==0)
  {
	  ControllerON=1;
	  controller_flag=1;
	  HAL_UART_Transmit(&huart2, "GO\n\r", 4,10);

  }else if(controller_flag==1 )
  {
	  ControllerON=0;
      controller_flag=0;
      HAL_UART_Transmit(&huart2, "STOP\n\r", 6,10);

   }



  /* USER CODE END EXTI15_10_IRQn 1 */
}

/**
  * @brief This function handles UART5 global interrupt.
  */
void UART5_IRQHandler(void)
{
  /* USER CODE BEGIN UART5_IRQn 0 */

  /* USER CODE END UART5_IRQn 0 */
  HAL_UART_IRQHandler(&huart5);
  /* USER CODE BEGIN UART5_IRQn 1 */


   /* if((UART5->SR & 0x0525U)!= RESET)
  	{

  	   RS485_RxBuffer[cnt]=UART5->DR;

  	   cnt++;

  	   if(cnt==10)
  		   {
  		      cnt=0;  //Rest counter
  		      for(int i=0; i<RS485_BUFFER_SIZE; i++)
  		      {
  		    	  //sprintf(HEX_Buffer,"0x%x ",RS485_RxBuffer[i]);
  		    	  //HAL_UART_Transmit(&huart2, HEX_Buffer, strlen(HEX_Buffer),1000);

  		      }
  		      motor_status_flag=0;
  		   }

  	}
*/

  /* USER CODE END UART5_IRQn 1 */
}

/**
  * @brief This function handles TIM6 global interrupt and DAC1, DAC2 underrun error interrupts.
  */
void TIM6_DAC_IRQHandler(void)
{
  /* USER CODE BEGIN TIM6_DAC_IRQn 0 */

  /* USER CODE END TIM6_DAC_IRQn 0 */
  HAL_TIM_IRQHandler(&htim6);
  /* USER CODE BEGIN TIM6_DAC_IRQn 1 */
  counter++;



  //HAL_GPIO_TogglePin(LD2_GPIO_Port,LD2_Pin); // fault detection


  if (counter == 2) {


      /*READ BNO085*/

	  /*Controller */


  }

     /*DDSM115 Control MODE*/

	  if(counter == 0  ) {
		  sendCurrentCommand(0x30, -iq);

	  }
	  if(counter == 1  ) {
		  sendCurrentCommand(0x10, iq);

	  }


     /*Distance*/
	  counter_dist++;
	  if(counter_dist==2)
	  {
	    Distnacex01();
	    Distnacex30();
	    counter_dist=0;
	  }

  /* USER CODE END TIM6_DAC_IRQn 1 */
}

/* USER CODE BEGIN 1 */

/* USER CODE END 1 */
