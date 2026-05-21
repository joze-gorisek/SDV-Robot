# MAIN Change Report

Reference project: `7_Florian`

Target project: `MAIN`

Date: 2026-05-19

## Summary

`MAIN` is based on `7_Florian`, with additional BNO086 support, corrected CyberGear and DDSM115 addressing, shared debug/state arrays, UART debug output, button debounce, and updated communication handling.

The original `7_Florian` motor/control structure is mostly preserved. The main functional additions are sensor integration and cleaner exported data for live debugging or later control work.

## Main Functional Changes

### BNO086 Integration

Added the full SH2/BNO086 sensor hub implementation:

- `MAIN/SH2 Sensorhub/Inc`
- `MAIN/SH2 Sensorhub/Src`

Main BNO086 application files:

- `SH2 Sensorhub/Src/demo_app.c`
- `SH2 Sensorhub/Inc/demo_app.h`

The BNO086 is initialized from `main.c` with:

```c
BNO_Init(&hi2c3, &huart2, INT_Pin);
```

BNO086 I2C address:

```c
#define BNO_ADDR (0x4B << 1)
```

BNO086 data is read in the main loop with:

```c
BNO_App(robot_data);
```

The BNO086 writes its data into `robot_data[10]` through `robot_data[19]`.

### BNO086 GPIO and Interrupts

Added BNO086 interrupt and reset pins in `main.h`:

```c
INT_Pin       = GPIO_PIN_1
INT_GPIO_Port = GPIOA
INT_EXTI_IRQn = EXTI1_IRQn

RST_Pin       = GPIO_PIN_15
RST_GPIO_Port = GPIOB
```

GPIO configuration added:

- `PA1` as falling-edge EXTI input with pull-up
- `PB15` as BNO reset output
- `EXTI1_IRQn` enabled

Added `EXTI1_IRQHandler()` in `stm32f4xx_it.c`:

```c
HAL_GPIO_EXTI_IRQHandler(INT_Pin);
```

The existing blue-button callback and BNO086 callback are merged into one `HAL_GPIO_EXTI_Callback()` because STM32 HAL allows only one callback with this name.

### I2C3 Configuration

`I2C3` is used for the BNO086:

- `PA8` = I2C3 SCL
- `PC9` = I2C3 SDA

The I2C speed was set to fast mode:

```c
hi2c3.Init.ClockSpeed = 400000;
```

This matches the working BNO086 demo behavior.

## CyberGear Changes

### CyberGear IDs

`7_Florian` used old/stale IDs in several places, such as:

```c
100, 9, 10, 11
```

`MAIN` now defines CyberGear IDs in `Core/Inc/CyberGear.h`:

```c
#define CYBER_HOST_ID     0xFEU
#define CYBER_MOTOR_1_ID  17U
#define CYBER_MOTOR_2_ID  18U
#define CYBER_MOTOR_3_ID  19U
#define CYBER_MOTOR_4_ID  20U
```

The helper movement code in `CyberMove.c` now uses these constants instead of hard-coded motor IDs.

### CyberGear Feedback Mapping

CyberGear feedback in `CyberGear.c` maps the new IDs to the existing variables:

```c
CYBER_MOTOR_1_ID -> MOTangle100
CYBER_MOTOR_2_ID -> MOTangle9
CYBER_MOTOR_3_ID -> MOTangle10
CYBER_MOTOR_4_ID -> MOTangle11
```

The old variable names were kept to minimize changes elsewhere in the project.

### CyberGear Sending

`MAIN` uses `CyberGear_Service()` in the main loop.

This service sends one CyberGear command when a CAN transmit mailbox is available:

```c
if (HAL_CAN_GetTxMailboxesFreeLevel(&hcan1) == 0U)
{
    return;
}
```

Then it cycles through the four motors:

```c
SetAngle(desired_angle,  CYBER_HOST_ID, CYBER_MOTOR_1_ID);
SetAngle(desired_angle2, CYBER_HOST_ID, CYBER_MOTOR_2_ID);
SetAngle(desired_angle3, CYBER_HOST_ID, CYBER_MOTOR_3_ID);
SetAngle(desired_angle4, CYBER_HOST_ID, CYBER_MOTOR_4_ID);
```

This avoids fixed blocking delays for CAN communication in `MAIN`.

## DDSM115 Changes

### DDSM115 IDs

The DDSM115 IDs in `MAIN` are:

```c
#define DDSM_LEFT_ID  0x10U
#define DDSM_RIGHT_ID 0x30U
```

The right ID remains `0x30`.

The left ID is updated to `0x10`.

### DDSM115 UART DMA Helper

`MAIN` uses a helper to start UART5 DMA only when UART5 is ready:

```c
static HAL_StatusTypeDef UART5_StartReceiveDmaIfReady(void)
{
    if (huart5.RxState == HAL_UART_STATE_READY)
    {
        return HAL_UART_Receive_DMA(&huart5, (uint8_t *)buffer485, sizeof(buffer485));
    }

    return HAL_BUSY;
}
```

This replaced direct repeated calls like:

```c
HAL_UART_Receive_DMA(&huart5, buffer485, 10);
```

### DDSM115 Feedback-Driven Sending

`MAIN` uses `DDSM_Service()` in the main loop.

The service sends one DDSM command at a time:

1. Send left motor command.
2. Wait for matching left feedback or timeout.
3. Send right motor command.
4. Wait for matching right feedback or timeout.

Timeout:

```c
#define DDSM_FEEDBACK_TIMEOUT_MS 10U
```

Feedback counters were added:

```c
ddsm_left_feedback_count
ddsm_right_feedback_count
ddsm_feedback_timeout_count
```

The UART callback calls:

```c
DDSM_HandleFeedback(DDSM_LEFT_ID);
DDSM_HandleFeedback(DDSM_RIGHT_ID);
```

### DDSM115 Sign Handling

Compared with `7_Florian`, `MAIN/Core/Src/DDSMove.c` has changed signs for mirrored wheel behavior.

Current `MAIN` behavior:

```c
distnacex01 = -1.0 * (...);
distnacex30 = (...);
DDSvelocity30 = -1.0 * DDSvelocityRadial30 * Rwheel;
```

This was done because one DDSM motor is mirror-mounted, so one physical wheel direction must be sign-corrected to represent forward robot motion.

## Button Logic

Button debounce was added:

```c
#define BUTTON_DEBOUNCE_MS 300U
```

The EXTI callback no longer directly changes motor commands. Instead, it sets:

```c
button_step_requested = 1;
```

The main loop processes the button request and updates:

- CyberGear desired angles
- DDSM current command variables

This keeps heavier logic out of the interrupt callback.

## Data Export Arrays

### `robot_data[20]`

Added global array:

```c
float robot_data[ROBOT_DATA_SIZE] = {0};
```

Layout:

```c
robot_data[0]  = DDSM left distance [m]
robot_data[1]  = DDSM right distance [m]
robot_data[2]  = DDSM left speed [m/s]
robot_data[3]  = DDSM right speed [m/s]

robot_data[4]  = CyberGear motor 1 position [rad]
robot_data[5]  = CyberGear motor 2 position [rad]
robot_data[6]  = CyberGear motor 3 position [rad]
robot_data[7]  = CyberGear motor 4 position [rad]
robot_data[8]  = DDSM left position [rad]
robot_data[9]  = DDSM right position [rad]

robot_data[10] = BNO roll [deg]
robot_data[11] = BNO pitch [deg]
robot_data[12] = BNO yaw [deg]
robot_data[13] = BNO gyro X [rad/s]
robot_data[14] = BNO gyro Y [rad/s]
robot_data[15] = BNO gyro Z [rad/s]
robot_data[16] = BNO quaternion W
robot_data[17] = BNO quaternion X
robot_data[18] = BNO quaternion Y
robot_data[19] = BNO quaternion Z
```

`RobotData_Update()` fills indices `0` through `9`.

`BNO_App(robot_data)` fills indices `10` through `19`.

### `state_data[4]`

Added global array:

```c
float state_data[4] = {0};
```

Layout:

```c
state_data[0] = robot position x [m]
state_data[1] = robot speed v [m/s]
state_data[2] = robot angle theta [deg]
state_data[3] = robot angular speed omega [rad/s]
```

Current calculation:

```c
state_data[0] = (distnacex01 + distnacex30) * 0.5f;
state_data[1] = (DDSvelocity01 + DDSvelocity30) * 0.5f;
state_data[2] = robot_data[10]; // BNO roll [deg]
state_data[3] = robot_data[13]; // BNO gyro X [rad/s]
```

This array is intended as the compact state vector:

```text
x = [position, speed, angle, angular speed]
```

It is only for monitoring/export at this stage and does not change controller behavior.

## UART Debug Output

Periodic USART2 debug output was added.

It prints:

- DDSM distances and speeds
- CyberGear positions
- DDSM positions
- BNO roll/pitch/yaw
- BNO gyro XYZ
- BNO quaternion WXYZ
- `state_data[0..3]`

The print interval is:

```c
#define DEBUG_PRINT_INTERVAL_MS 500U
```

The output is timed with `HAL_GetTick()` instead of a blocking debug delay.

## Main Loop Behavior

Current main loop order in `MAIN`:

1. Read/process BNO086:

   ```c
   BNO_App(robot_data);
   ```

2. Update full debug data:

   ```c
   RobotData_Update();
   ```

3. Update compact state vector:

   ```c
   StateData_Update();
   ```

4. Print UART debug data every 500 ms.

5. Start UART5 DMA receive if ready.

6. Process debounced button request.

7. Service DDSM communication:

   ```c
   DDSM_Service();
   ```

8. Service CyberGear communication:

   ```c
   CyberGear_Service();
   ```

## Files Added or Significantly Changed

### Added

- `MAIN/SH2 Sensorhub/Inc/*`
- `MAIN/SH2 Sensorhub/Src/*`
- `MAIN/CHANGE_REPORT.md`

### Significantly changed

- `MAIN/Core/Src/main.c`
- `MAIN/Core/Inc/main.h`
- `MAIN/Core/Src/stm32f4xx_it.c`
- `MAIN/Core/Inc/stm32f4xx_it.h`
- `MAIN/Core/Src/stm32f4xx_hal_msp.c`
- `MAIN/Core/Src/DDSMove.c`
- `MAIN/Core/Inc/DDSM115.h`
- `MAIN/Core/Inc/CyberGear.h`
- `MAIN/Core/Src/CyberGear.c`
- `MAIN/Core/Src/CyberMove.c`
- `MAIN/BipedRobot.ioc`
- `MAIN/.cproject`

## Wiring Assumptions

BNO086 wiring:

```text
PA8   -> SCL
PC9   -> SDA
PA1   -> INT
PB15  -> RST
3.3V  -> VCC
GND   -> GND
```

DDSM115:

```text
Left  ID = 0x10
Right ID = 0x30
UART5 RS485
```

CyberGear:

```text
Host ID = 0xFE
Motor IDs = 17, 18, 19, 20
CAN1
```

## Notes and Open Items

- `MAIN` uses feedback-driven/non-blocking communication for DDSM and CyberGear.
- `MAIN - WithDelay` is a separate project variant that uses explicit `HAL_Delay()` communication timing in the main loop.
- The old variable names such as `MOTangle100`, `MOTangle9`, and `DDSangle01` were kept for minimal code disruption.
- `state_data` and `robot_data` are monitoring/export arrays only. They do not yet drive the controller.
- The DDSM distance still depends on the existing absolute encoder plus revolution-count logic from `7_Florian`, with sign corrections added in `MAIN`.
- If distance should reset exactly to zero at boot or at a button event, a dedicated odometry reset/baseline function should be added.

