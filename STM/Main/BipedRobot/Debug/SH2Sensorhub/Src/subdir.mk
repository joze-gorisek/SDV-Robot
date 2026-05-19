################################################################################
# Automatically-generated file. Do not edit!
# Toolchain: GNU Tools for STM32 (14.3.rel1)
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../SH2Sensorhub/Src/demo_app.c \
../SH2Sensorhub/Src/euler.c \
../SH2Sensorhub/Src/sh2.c \
../SH2Sensorhub/Src/sh2_SensorValue.c \
../SH2Sensorhub/Src/sh2_util.c \
../SH2Sensorhub/Src/shtp.c 

OBJS += \
./SH2Sensorhub/Src/demo_app.o \
./SH2Sensorhub/Src/euler.o \
./SH2Sensorhub/Src/sh2.o \
./SH2Sensorhub/Src/sh2_SensorValue.o \
./SH2Sensorhub/Src/sh2_util.o \
./SH2Sensorhub/Src/shtp.o 

C_DEPS += \
./SH2Sensorhub/Src/demo_app.d \
./SH2Sensorhub/Src/euler.d \
./SH2Sensorhub/Src/sh2.d \
./SH2Sensorhub/Src/sh2_SensorValue.d \
./SH2Sensorhub/Src/sh2_util.d \
./SH2Sensorhub/Src/shtp.d 


# Each subdirectory must supply rules for building sources it contributes
SH2Sensorhub/Src/%.o SH2Sensorhub/Src/%.su SH2Sensorhub/Src/%.cyclo: ../SH2Sensorhub/Src/%.c SH2Sensorhub/Src/subdir.mk
	arm-none-eabi-gcc -gdwarf-4 "$<" -mcpu=cortex-m4 -std=gnu11 -g3 -DDEBUG -DUSE_HAL_DRIVER -DSTM32F446xx -c -I../Core/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc/Legacy -I../Drivers/CMSIS/Device/ST/STM32F4xx/Include -I../Drivers/CMSIS/Include -I../SH2Sensorhub/Inc -O0 -ffunction-sections -fdata-sections -Wall -fstack-usage -fcyclomatic-complexity -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" --specs=nano.specs -mfpu=fpv4-sp-d16 -mfloat-abi=hard -mthumb -o "$@"

clean: clean-SH2Sensorhub-2f-Src

clean-SH2Sensorhub-2f-Src:
	-$(RM) ./SH2Sensorhub/Src/demo_app.cyclo ./SH2Sensorhub/Src/demo_app.d ./SH2Sensorhub/Src/demo_app.o ./SH2Sensorhub/Src/demo_app.su ./SH2Sensorhub/Src/euler.cyclo ./SH2Sensorhub/Src/euler.d ./SH2Sensorhub/Src/euler.o ./SH2Sensorhub/Src/euler.su ./SH2Sensorhub/Src/sh2.cyclo ./SH2Sensorhub/Src/sh2.d ./SH2Sensorhub/Src/sh2.o ./SH2Sensorhub/Src/sh2.su ./SH2Sensorhub/Src/sh2_SensorValue.cyclo ./SH2Sensorhub/Src/sh2_SensorValue.d ./SH2Sensorhub/Src/sh2_SensorValue.o ./SH2Sensorhub/Src/sh2_SensorValue.su ./SH2Sensorhub/Src/sh2_util.cyclo ./SH2Sensorhub/Src/sh2_util.d ./SH2Sensorhub/Src/sh2_util.o ./SH2Sensorhub/Src/sh2_util.su ./SH2Sensorhub/Src/shtp.cyclo ./SH2Sensorhub/Src/shtp.d ./SH2Sensorhub/Src/shtp.o ./SH2Sensorhub/Src/shtp.su

.PHONY: clean-SH2Sensorhub-2f-Src

