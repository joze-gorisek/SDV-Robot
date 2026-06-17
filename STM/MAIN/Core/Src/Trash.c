/*
 * Trash.c
 *
 *  Created on: Apr 1, 2025
 *      Author: primozs
 */


/**
  * @brief This function handles UART4 global interrupt.
  */
//
//void UART4_IRQHandler(void)
//{
//  /* USER CODE BEGIN UART4_IRQn 0 */
//
//  /* USER CODE END UART4_IRQn 0 */
//  HAL_UART_IRQHandler(&huart4);
//  /* USER CODE BEGIN UART4_IRQn 1 */
//
//   if((UART4->SR & 0x0525U)!= RESET)
//	{
//
//	   RS485_RxBuffer[cnt]=UART4->DR;
//
//	   cnt++;
//
//	   if(cnt==10)
//		   {
//		      cnt=0;  //Rest counter
//		      for(int i=0; i<RS485_BUFFER_SIZE; i++)
//		      {
//		    	 // sprintf(HEX_Buffer,"0x%x ",RS485_RxBuffer[i]);
//		    	 // HAL_UART_Transmit(&huart2, HEX_Buffer, strlen(HEX_Buffer),1000);
//
//		      }
//		      motor_status_flag=0;
//		   }
//
//	}
//
//
//  /* USER CODE END UART4_IRQn 1 */
//}
