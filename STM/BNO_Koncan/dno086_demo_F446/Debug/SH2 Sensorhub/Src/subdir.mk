################################################################################
# Automatically-generated file. Do not edit!
# Toolchain: GNU Tools for STM32 (14.3.rel1)
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../SH2\ Sensorhub/Src/demo_app.c \
../SH2\ Sensorhub/Src/euler.c \
../SH2\ Sensorhub/Src/sh2.c \
../SH2\ Sensorhub/Src/sh2_SensorValue.c \
../SH2\ Sensorhub/Src/sh2_util.c \
../SH2\ Sensorhub/Src/shtp.c 

OBJS += \
./SH2\ Sensorhub/Src/demo_app.o \
./SH2\ Sensorhub/Src/euler.o \
./SH2\ Sensorhub/Src/sh2.o \
./SH2\ Sensorhub/Src/sh2_SensorValue.o \
./SH2\ Sensorhub/Src/sh2_util.o \
./SH2\ Sensorhub/Src/shtp.o 

C_DEPS += \
./SH2\ Sensorhub/Src/demo_app.d \
./SH2\ Sensorhub/Src/euler.d \
./SH2\ Sensorhub/Src/sh2.d \
./SH2\ Sensorhub/Src/sh2_SensorValue.d \
./SH2\ Sensorhub/Src/sh2_util.d \
./SH2\ Sensorhub/Src/shtp.d 


# Each subdirectory must supply rules for building sources it contributes
SH2\ Sensorhub/Src/demo_app.o: ../SH2\ Sensorhub/Src/demo_app.c SH2\ Sensorhub/Src/subdir.mk
	arm-none-eabi-gcc "$<" -mcpu=cortex-m4 -std=gnu11 -g3 -DDEBUG -DUSE_HAL_DRIVER -DSTM32F446xx -c -I../Core/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc/Legacy -I../Drivers/CMSIS/Device/ST/STM32F4xx/Include -I../Drivers/CMSIS/Include -I"../SH2 Sensorhub/Inc" -O0 -ffunction-sections -fdata-sections -Wall -fstack-usage -fcyclomatic-complexity -MMD -MP -MF"SH2 Sensorhub/Src/demo_app.d" -MT"$@" --specs=nano.specs -mfpu=fpv4-sp-d16 -mfloat-abi=hard -mthumb -o "$@"
SH2\ Sensorhub/Src/euler.o: ../SH2\ Sensorhub/Src/euler.c SH2\ Sensorhub/Src/subdir.mk
	arm-none-eabi-gcc "$<" -mcpu=cortex-m4 -std=gnu11 -g3 -DDEBUG -DUSE_HAL_DRIVER -DSTM32F446xx -c -I../Core/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc/Legacy -I../Drivers/CMSIS/Device/ST/STM32F4xx/Include -I../Drivers/CMSIS/Include -I"../SH2 Sensorhub/Inc" -O0 -ffunction-sections -fdata-sections -Wall -fstack-usage -fcyclomatic-complexity -MMD -MP -MF"SH2 Sensorhub/Src/euler.d" -MT"$@" --specs=nano.specs -mfpu=fpv4-sp-d16 -mfloat-abi=hard -mthumb -o "$@"
SH2\ Sensorhub/Src/sh2.o: ../SH2\ Sensorhub/Src/sh2.c SH2\ Sensorhub/Src/subdir.mk
	arm-none-eabi-gcc "$<" -mcpu=cortex-m4 -std=gnu11 -g3 -DDEBUG -DUSE_HAL_DRIVER -DSTM32F446xx -c -I../Core/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc/Legacy -I../Drivers/CMSIS/Device/ST/STM32F4xx/Include -I../Drivers/CMSIS/Include -I"../SH2 Sensorhub/Inc" -O0 -ffunction-sections -fdata-sections -Wall -fstack-usage -fcyclomatic-complexity -MMD -MP -MF"SH2 Sensorhub/Src/sh2.d" -MT"$@" --specs=nano.specs -mfpu=fpv4-sp-d16 -mfloat-abi=hard -mthumb -o "$@"
SH2\ Sensorhub/Src/sh2_SensorValue.o: ../SH2\ Sensorhub/Src/sh2_SensorValue.c SH2\ Sensorhub/Src/subdir.mk
	arm-none-eabi-gcc "$<" -mcpu=cortex-m4 -std=gnu11 -g3 -DDEBUG -DUSE_HAL_DRIVER -DSTM32F446xx -c -I../Core/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc/Legacy -I../Drivers/CMSIS/Device/ST/STM32F4xx/Include -I../Drivers/CMSIS/Include -I"../SH2 Sensorhub/Inc" -O0 -ffunction-sections -fdata-sections -Wall -fstack-usage -fcyclomatic-complexity -MMD -MP -MF"SH2 Sensorhub/Src/sh2_SensorValue.d" -MT"$@" --specs=nano.specs -mfpu=fpv4-sp-d16 -mfloat-abi=hard -mthumb -o "$@"
SH2\ Sensorhub/Src/sh2_util.o: ../SH2\ Sensorhub/Src/sh2_util.c SH2\ Sensorhub/Src/subdir.mk
	arm-none-eabi-gcc "$<" -mcpu=cortex-m4 -std=gnu11 -g3 -DDEBUG -DUSE_HAL_DRIVER -DSTM32F446xx -c -I../Core/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc/Legacy -I../Drivers/CMSIS/Device/ST/STM32F4xx/Include -I../Drivers/CMSIS/Include -I"../SH2 Sensorhub/Inc" -O0 -ffunction-sections -fdata-sections -Wall -fstack-usage -fcyclomatic-complexity -MMD -MP -MF"SH2 Sensorhub/Src/sh2_util.d" -MT"$@" --specs=nano.specs -mfpu=fpv4-sp-d16 -mfloat-abi=hard -mthumb -o "$@"
SH2\ Sensorhub/Src/shtp.o: ../SH2\ Sensorhub/Src/shtp.c SH2\ Sensorhub/Src/subdir.mk
	arm-none-eabi-gcc "$<" -mcpu=cortex-m4 -std=gnu11 -g3 -DDEBUG -DUSE_HAL_DRIVER -DSTM32F446xx -c -I../Core/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc/Legacy -I../Drivers/CMSIS/Device/ST/STM32F4xx/Include -I../Drivers/CMSIS/Include -I"../SH2 Sensorhub/Inc" -O0 -ffunction-sections -fdata-sections -Wall -fstack-usage -fcyclomatic-complexity -MMD -MP -MF"SH2 Sensorhub/Src/shtp.d" -MT"$@" --specs=nano.specs -mfpu=fpv4-sp-d16 -mfloat-abi=hard -mthumb -o "$@"

clean: clean-SH2-20-Sensorhub-2f-Src

clean-SH2-20-Sensorhub-2f-Src:
	-$(RM) ./SH2\ Sensorhub/Src/demo_app.cyclo ./SH2\ Sensorhub/Src/demo_app.d ./SH2\ Sensorhub/Src/demo_app.o ./SH2\ Sensorhub/Src/demo_app.su ./SH2\ Sensorhub/Src/euler.cyclo ./SH2\ Sensorhub/Src/euler.d ./SH2\ Sensorhub/Src/euler.o ./SH2\ Sensorhub/Src/euler.su ./SH2\ Sensorhub/Src/sh2.cyclo ./SH2\ Sensorhub/Src/sh2.d ./SH2\ Sensorhub/Src/sh2.o ./SH2\ Sensorhub/Src/sh2.su ./SH2\ Sensorhub/Src/sh2_SensorValue.cyclo ./SH2\ Sensorhub/Src/sh2_SensorValue.d ./SH2\ Sensorhub/Src/sh2_SensorValue.o ./SH2\ Sensorhub/Src/sh2_SensorValue.su ./SH2\ Sensorhub/Src/sh2_util.cyclo ./SH2\ Sensorhub/Src/sh2_util.d ./SH2\ Sensorhub/Src/sh2_util.o ./SH2\ Sensorhub/Src/sh2_util.su ./SH2\ Sensorhub/Src/shtp.cyclo ./SH2\ Sensorhub/Src/shtp.d ./SH2\ Sensorhub/Src/shtp.o ./SH2\ Sensorhub/Src/shtp.su

.PHONY: clean-SH2-20-Sensorhub-2f-Src

