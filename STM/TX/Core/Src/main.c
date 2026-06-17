/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2026 STMicroelectronics.
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
#include "MRF24J40.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */
#define TX_COMMAND_MAX_LEN 32U

/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
SPI_HandleTypeDef hspi1;

UART_HandleTypeDef huart2;

/* USER CODE BEGIN PV */

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_SPI1_Init(void);
static void MX_USART2_UART_Init(void);
/* USER CODE BEGIN PFP */
static void TX_SendPcCommand(uint8_t *cmd, uint8_t len);
static uint8_t TX_CommandIsValid(const uint8_t *cmd, uint8_t len);
static uint8_t TX_WaitMrfTxOk(uint8_t previous_tx_count);
static void TX_ProcessPcCommand(uint8_t *cmd, uint8_t *len);

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

static void TX_SendPcCommand(uint8_t *cmd, uint8_t len)
{
  if (len == 0U)
  {
	  return;
  }

  uint8_t *tx_data = cmd;
  uint8_t tx_len = len;
  uint8_t repeats = 1U;
  uint8_t tx_before = flag_got_tx;
  uint8_t tx_ok = 0U;

  for (uint8_t n = 0; n < repeats; n++)
  {
	  Mrf24j_send16(0xffff, (char*)tx_data, tx_len);
	  tx_ok = TX_WaitMrfTxOk(tx_before);
	  tx_before = flag_got_tx;
  }

  if (tx_ok)
  {
	  HAL_GPIO_TogglePin(LD2_GPIO_Port, LD2_Pin);
	  HAL_UART_Transmit(&huart2, (uint8_t*)"TX OK ", 6, 100);
	  HAL_UART_Transmit(&huart2, cmd, len, 100);
	  HAL_UART_Transmit(&huart2, (uint8_t*)"\r\n", 2, 100);
  }
  else
  {
	  char msg[32];
	  uint8_t txstat = Mrf24j_read_short(MRF_TXSTAT);
	  if ((txstat & ((1 << TXNSTAT) | (1 << CCAFAIL))) == 0U)
	  {
		  HAL_GPIO_TogglePin(LD2_GPIO_Port, LD2_Pin);
		  HAL_UART_Transmit(&huart2, (uint8_t*)"TX OK ", 6, 100);
		  HAL_UART_Transmit(&huart2, cmd, len, 100);
		  HAL_UART_Transmit(&huart2, (uint8_t*)"\r\n", 2, 100);
	  }
	  else
	  {
		  int msg_len = snprintf(msg, sizeof(msg), "TX FAIL stat=0x%02X\r\n", txstat);
		  if (msg_len > 0)
		  {
			  HAL_UART_Transmit(&huart2, (uint8_t*)msg, (uint16_t)msg_len, 100);
		  }
	  }
  }
}

static uint8_t TX_WaitMrfTxOk(uint8_t previous_tx_count)
{
  uint32_t start = HAL_GetTick();

  while ((HAL_GetTick() - start) < 500U)
  {
	  if (HAL_GPIO_ReadPin(MRF_INT_GPIO_Port, MRF_INT_Pin) == GPIO_PIN_RESET)
	  {
		  Mrf24j_interrupt_handler();
	  }

	  uint8_t intstat = Mrf24j_read_short(MRF_INTSTAT);
	  if (intstat & MRF_I_TXNIF)
	  {
		  uint8_t txstat = Mrf24j_read_short(MRF_TXSTAT);
		  flag_got_tx++;
		  tx_info.tx_ok = ((txstat & ((1 << TXNSTAT) | (1 << CCAFAIL))) == 0U);
		  tx_info.retries = txstat >> 6;
		  tx_info.channel_busy = (txstat & (1 << CCAFAIL));
	  }

	  if (flag_got_tx != previous_tx_count)
	  {
		  return tx_info.tx_ok ? 1U : 0U;
	  }
  }

  uint8_t txstat = Mrf24j_read_short(MRF_TXSTAT);
  return ((txstat & ((1 << TXNSTAT) | (1 << CCAFAIL))) == 0U) ? 1U : 0U;
}

static uint8_t TX_CommandIsValid(const uint8_t *cmd, uint8_t len)
{
  char text[TX_COMMAND_MAX_LEN];

  if ((cmd == NULL) || (len == 0U) || (len >= TX_COMMAND_MAX_LEN))
  {
	  return 0U;
  }

  memcpy(text, cmd, len);
  text[len] = '\0';

  char *p = text;
  while (isspace((unsigned char)*p))
  {
	  p++;
  }

  char *end = p;
  (void)strtol(p, &end, 10);
  if (end == p)
  {
	  return 0U;
  }

  while (isspace((unsigned char)*end))
  {
	  end++;
  }

  if (((end[0] != 'm') && (end[0] != 'M')) || ((end[1] != 'm') && (end[1] != 'M')))
  {
	  return 0U;
  }
  end += 2;

  while (isspace((unsigned char)*end))
  {
	  end++;
  }

  if ((end[0] != '0') && (end[0] != '1'))
  {
	  return 0U;
  }
  end++;

  while (isspace((unsigned char)*end))
  {
	  end++;
  }

  return (*end == '\0') ? 1U : 0U;
}

static void TX_ProcessPcCommand(uint8_t *cmd, uint8_t *len)
{
  if (*len == 0U)
  {
	  return;
  }

  HAL_UART_Transmit(&huart2, (uint8_t*)"RX ", 3, 100);
  HAL_UART_Transmit(&huart2, cmd, *len, 100);
  HAL_UART_Transmit(&huart2, (uint8_t*)"\r\n", 2, 100);

  if (TX_CommandIsValid(cmd, *len))
  {
	  TX_SendPcCommand(cmd, *len);
  }
  else
  {
	  HAL_UART_Transmit(&huart2, (uint8_t*)"ERR\r\n", 5, 100);
  }

  *len = 0U;
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

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_SPI1_Init();
  MX_USART2_UART_Init();
  /* USER CODE BEGIN 2 */

  // MRF SETUP
  MRF_InitSPI(&hspi1, &huart2, 0);
  MRF_InitGPIO(MRF_CS_GPIO_Port, MRF_CS_Pin, MRF_CE_GPIO_Port, MRF_CE_Pin, MRF_INT_GPIO_Port, MRF_INT_Pin);

  Mrf24j_init(20u);

  {
	  uint16_t pan = Mrf24j_get_pan();
	  uint8_t rfcon0 = Mrf24j_read_long(MRF_RFCON0);
	  uint8_t rfcon6 = Mrf24j_read_long(MRF_RFCON6);
	  uint8_t rxmcr = Mrf24j_read_short(MRF_RXMCR);
	  char msg[80];
	  int msg_len = snprintf(msg, sizeof(msg),
			  	  	  	  	 "TX MRF PAN=0x%04X RFCON0=0x%02X RFCON6=0x%02X RXMCR=0x%02X\r\n",
							 pan, rfcon0, rfcon6, rxmcr);
	  if (msg_len > 0)
	  {
		  HAL_UART_Transmit(&huart2, (uint8_t*)msg, (uint16_t)msg_len, HAL_MAX_DELAY);
	  }
  }

  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  HAL_UART_Transmit(&huart2, (uint8_t*)"MRF PC-command ready: use 300mm 0!\r\n", 36, HAL_MAX_DELAY);
  uint8_t pc_cmd[TX_COMMAND_MAX_LEN];
  uint8_t pc_cmd_len = 0U;
  while (1)
  {
	  uint8_t ch;

	  while (HAL_UART_Receive(&huart2, &ch, 1, 50) == HAL_OK)
	  {
		  if (ch == '!')
		  {
			  TX_ProcessPcCommand(pc_cmd, &pc_cmd_len);
		  }
		  else if ((ch == '\r') || (ch == '\n'))
		  {
			  /* Commands are ended only with '!'. Ignore terminal line endings. */
		  }
		  else if ((ch == '\b') || (ch == 0x7FU))
		  {
			  if (pc_cmd_len > 0U)
			  {
				  pc_cmd_len--;
			  }
		  }
		  else if (pc_cmd_len < (TX_COMMAND_MAX_LEN - 1U))
		  {
			  pc_cmd[pc_cmd_len++] = ch;
		  }
		  else
		  {
			  pc_cmd_len = 0U;
			  HAL_UART_Transmit(&huart2, (uint8_t*)"ERR len\r\n", 9, 100);
		  }
	  }

	  Mrf24j_poll_uart();

	  HAL_Delay(1);
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
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
  hspi1.Init.BaudRatePrescaler = SPI_BAUDRATEPRESCALER_256;
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

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(LD2_GPIO_Port, LD2_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(MRF_CE_GPIO_Port, MRF_CE_Pin, GPIO_PIN_RESET);
  HAL_GPIO_WritePin(MRF_CS_GPIO_Port, MRF_CS_Pin, GPIO_PIN_SET);

  /*Configure GPIO pin : B1_Pin */
  GPIO_InitStruct.Pin = B1_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_FALLING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(B1_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : LD2_Pin */
  GPIO_InitStruct.Pin = LD2_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(LD2_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : MRF_INT_Pin */
  GPIO_InitStruct.Pin = MRF_INT_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_FALLING;
  GPIO_InitStruct.Pull = GPIO_PULLUP;
  HAL_GPIO_Init(MRF_INT_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pins : MRF_CE_Pin MRF_CS_Pin */
  GPIO_InitStruct.Pin = MRF_CE_Pin|MRF_CS_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOB, &GPIO_InitStruct);

  /* EXTI interrupt init*/
  HAL_NVIC_SetPriority(EXTI15_10_IRQn, 0, 0);
  HAL_NVIC_EnableIRQ(EXTI15_10_IRQn);

  /* USER CODE BEGIN MX_GPIO_Init_2 */

  /* USER CODE END MX_GPIO_Init_2 */
}

/* USER CODE BEGIN 4 */

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
#ifdef USE_FULL_ASSERT
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
