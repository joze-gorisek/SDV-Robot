/*
 * demo_app.c
 */

#include "../../SH2Sensorhub/Inc/demo_app.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../../SH2Sensorhub/Inc/sh2.h"
#include "../../SH2Sensorhub/Inc/sh2_err.h"
#include "../../SH2Sensorhub/Inc/sh2_SensorValue.h"

// ----------------------------------------------------------------
// Module-level state
// ----------------------------------------------------------------
static I2C_HandleTypeDef  *_i2c;
static UART_HandleTypeDef *_uart;
static uint16_t            _intPin;
static sh2_Hal_t           _hal;

// Orientation
volatile float bno_qw, bno_qx, bno_qy, bno_qz;
volatile float bno_roll, bno_pitch, bno_yaw;

// Angular velocity (rad/s)
volatile float bno_gx, bno_gy, bno_gz;

volatile bool bnoIntFlag = false;

// Timing
volatile uint32_t t_old = 0;
volatile uint32_t t_new = 0;
volatile uint32_t delta_t    = 0;

// Retarget printf to UART

// ----------------------------------------------------------------
// DWT
// ----------------------------------------------------------------
void DWT_Init(void)
{
    CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
    DWT->CYCCNT = 0;
    DWT->CTRL  |= DWT_CTRL_CYCCNTENA_Msk;
}

uint32_t DWT_GetMicros(void)
{
	return (uint32_t)(DWT->CYCCNT / (SystemCoreClock / 1000000U));
}

// ----------------------------------------------------------------
// HAL stubs
// ----------------------------------------------------------------
static int hal_open(sh2_Hal_t *self)  { (void)self; return SH2_OK; }
static void hal_close(sh2_Hal_t *self) { (void)self; }

static int hal_read(sh2_Hal_t *self, uint8_t *buf, unsigned int len,
                    uint32_t *t_us)
{
	(void)self;

	    // Phase 1: read 4-byte header only
	    HAL_StatusTypeDef st = HAL_I2C_Master_Receive(_i2c, BNO_ADDR, buf, 4, 100);
	    if (st != HAL_OK) return 0;

	    uint16_t pktLen = ((uint16_t)(buf[1] & 0x7F) << 8) | buf[0];
	    uint8_t  chan   = buf[2];

	    //printf("ch=%u pktLen=%u\r\n", chan, pktLen);

	    uint32_t t_start = DWT_GetMicros();
	    while(DWT_GetMicros()-t_start<1000);

	    if (pktLen == 0) return 0;
	    if (pktLen < 4)  return 0;
	    if (pktLen > (uint16_t)len) pktLen = (uint16_t)len;

	    // Phase 2: re-read the ENTIRE packet from byte 0
	    // (BNO086 resets its output pointer on every new I2C START)
	    st = HAL_I2C_Master_Receive(_i2c, BNO_ADDR, buf, pktLen, 100);
	    if (st != HAL_OK) return 0;

	    if (t_us) *t_us = DWT_GetMicros();;
	    return (int)pktLen;
}

static int hal_write(sh2_Hal_t *self, uint8_t *buf, unsigned int len)
{
    (void)self;
    if (HAL_I2C_Master_Transmit(_i2c, BNO_ADDR, buf, len, 100) != HAL_OK)
        return 0;
    return (int)len;
}

static uint32_t hal_getTimeUs(sh2_Hal_t *self)
{
    (void)self;
    return DWT_GetMicros();
}

// ----------------------------------------------------------------
// Sensor callback
// ----------------------------------------------------------------
static void sensorCallback(void *cookie, sh2_SensorEvent_t *event)
{
    (void)cookie;
    sh2_SensorValue_t val;

    if (sh2_decodeSensorEvent(&val, event) != SH2_OK)
        return;

    switch (val.sensorId)
    {
        case SH2_GAME_ROTATION_VECTOR:
            bno_qw = val.un.gameRotationVector.real;  // ✓ correct union
            bno_qx = val.un.gameRotationVector.i;
            bno_qy = val.un.gameRotationVector.j;
            bno_qz = val.un.gameRotationVector.k;

            {
                float sinr_cosp = 2.0f * (bno_qw * bno_qx + bno_qy * bno_qz);
                float cosr_cosp = 1.0f - 2.0f * (bno_qx * bno_qx + bno_qy * bno_qy);
                bno_roll  = atan2f(sinr_cosp, cosr_cosp) * (180.0f / M_PI);

                float sinp = 2.0f * (bno_qw * bno_qy - bno_qz * bno_qx);
                sinp = fmaxf(-1.0f, fminf(1.0f, sinp));
                bno_pitch = asinf(sinp) * (180.0f / M_PI);

                float siny_cosp = 2.0f * (bno_qw * bno_qz + bno_qx * bno_qy);
                float cosy_cosp = 1.0f - 2.0f * (bno_qy * bno_qy + bno_qz * bno_qz);
                bno_yaw   = atan2f(siny_cosp, cosy_cosp) * (180.0f / M_PI);
            }
            break;

        case SH2_GYROSCOPE_CALIBRATED:
            bno_gx = val.un.gyroscope.x;
            bno_gy = val.un.gyroscope.y;
            bno_gz = val.un.gyroscope.z;
            break;
    }
}

// ----------------------------------------------------------------
// Event callback — re-enables reports after BNO reset
// ----------------------------------------------------------------
static void eventCallback(void *cookie, sh2_AsyncEvent_t *event)
{
    (void)cookie;
    if (event->eventId == SH2_RESET)
    {
        sh2_SensorConfig_t cfg = {0};
        cfg.reportInterval_us = 5000;
        sh2_setSensorConfig(SH2_GAME_ROTATION_VECTOR, &cfg);
        sh2_setSensorConfig(SH2_GYROSCOPE_CALIBRATED, &cfg);  // ✓ restore gyro too
    }
}

// ----------------------------------------------------------------
// Init
// ----------------------------------------------------------------
void BNO_Init(I2C_HandleTypeDef *hi2c, UART_HandleTypeDef *huart,
              uint16_t intPin)
{
    DWT_Init();
    _i2c    = hi2c;
    _uart   = huart;
    _intPin = intPin;

    HAL_GPIO_WritePin(RST_GPIO_Port, RST_Pin, GPIO_PIN_RESET);
    HAL_Delay(10);
    HAL_GPIO_WritePin(RST_GPIO_Port, RST_Pin, GPIO_PIN_SET);
    HAL_Delay(300);

    _hal.open      = hal_open;
    _hal.close     = hal_close;
    _hal.read      = hal_read;
    _hal.write     = hal_write;
    _hal.getTimeUs = hal_getTimeUs;

    int rc = sh2_open(&_hal, eventCallback, NULL);
    if (rc != SH2_OK) { printf("sh2_open failed: %d\r\n", rc); return; }

    sh2_setSensorCallback(sensorCallback, NULL);

    for (int i = 0; i < 20; i++) { sh2_service(); HAL_Delay(10); }

    sh2_SensorConfig_t cfg = {0};
    cfg.reportInterval_us = 5000;
    rc  = sh2_setSensorConfig(SH2_GAME_ROTATION_VECTOR, &cfg);
    for (int i = 0; i < 10; i++) { sh2_service(); HAL_Delay(5); }
    rc |= sh2_setSensorConfig(SH2_GYROSCOPE_CALIBRATED, &cfg);
    for (int i = 0; i < 10; i++) { sh2_service(); HAL_Delay(5); }

    if (rc != SH2_OK) printf("setSensorConfig failed: %d\r\n", rc);
    else              printf("BNO086 ready.\r\n");
}

// ----------------------------------------------------------------
// App tick — call from main loop
// ----------------------------------------------------------------


void BNO_App(float vector[4])
{
    if (!bnoIntFlag)
        return;

    uint32_t timeout = DWT_GetMicros();

    do {
        sh2_service();
    }
    while (HAL_GPIO_ReadPin(INT_GPIO_Port, _intPin) == GPIO_PIN_RESET &&
           (DWT_GetMicros() - timeout) < 10000);

    bnoIntFlag = false;

    vector[2] = bno_roll;
    vector[3] = bno_gx;

    printf("R:%6.2f | Gx:%6.3f\r\n",bno_roll,bno_gx);

    HAL_NVIC_ClearPendingIRQ(EXTI1_IRQn);
    HAL_NVIC_EnableIRQ(EXTI1_IRQn);
}



// ----------------------------------------------------------------
// EXTI callback
// ----------------------------------------------------------------
void HAL_GPIO_EXTI_Callback(uint16_t GPIO_Pin)
{
    if (GPIO_Pin == _intPin &&
        HAL_GPIO_ReadPin(INT_GPIO_Port, _intPin) == GPIO_PIN_RESET)
    {
        t_new = DWT_GetMicros();
        delta_t    = t_new - t_old;
        t_old = t_new;

        bnoIntFlag = true;
        HAL_NVIC_DisableIRQ(EXTI1_IRQn);
    }
}
