################################################################################
# Automatically-generated file. Do not edit!
# Toolchain: GNU Tools for STM32 (14.3.rel1)
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../Core/Src/Controllers.c \
../Core/Src/CyberGear.c \
../Core/Src/CyberMove.c \
../Core/Src/DDSM115.c \
../Core/Src/DDSMove.c \
../Core/Src/Trash.c \
../Core/Src/kalman.c \
../Core/Src/main.c \
../Core/Src/mpu6050.c \
../Core/Src/stm32f4xx_hal_msp.c \
../Core/Src/stm32f4xx_it.c \
../Core/Src/syscalls.c \
../Core/Src/sysmem.c \
../Core/Src/system_stm32f4xx.c 

OBJS += \
./Core/Src/Controllers.o \
./Core/Src/CyberGear.o \
./Core/Src/CyberMove.o \
./Core/Src/DDSM115.o \
./Core/Src/DDSMove.o \
./Core/Src/Trash.o \
./Core/Src/kalman.o \
./Core/Src/main.o \
./Core/Src/mpu6050.o \
./Core/Src/stm32f4xx_hal_msp.o \
./Core/Src/stm32f4xx_it.o \
./Core/Src/syscalls.o \
./Core/Src/sysmem.o \
./Core/Src/system_stm32f4xx.o 

C_DEPS += \
./Core/Src/Controllers.d \
./Core/Src/CyberGear.d \
./Core/Src/CyberMove.d \
./Core/Src/DDSM115.d \
./Core/Src/DDSMove.d \
./Core/Src/Trash.d \
./Core/Src/kalman.d \
./Core/Src/main.d \
./Core/Src/mpu6050.d \
./Core/Src/stm32f4xx_hal_msp.d \
./Core/Src/stm32f4xx_it.d \
./Core/Src/syscalls.d \
./Core/Src/sysmem.d \
./Core/Src/system_stm32f4xx.d 


# Each subdirectory must supply rules for building sources it contributes
Core/Src/%.o Core/Src/%.su Core/Src/%.cyclo: ../Core/Src/%.c Core/Src/subdir.mk
	arm-none-eabi-gcc -gdwarf-4 "$<" -mcpu=cortex-m4 -std=gnu11 -g3 -DDEBUG -DUSE_HAL_DRIVER -DSTM32F446xx -c -I../Core/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc -I../Drivers/STM32F4xx_HAL_Driver/Inc/Legacy -I../Drivers/CMSIS/Device/ST/STM32F4xx/Include -I../Drivers/CMSIS/Include -O0 -ffunction-sections -fdata-sections -Wall -fstack-usage -fcyclomatic-complexity -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" --specs=nano.specs -mfpu=fpv4-sp-d16 -mfloat-abi=hard -mthumb -o "$@"

clean: clean-Core-2f-Src

clean-Core-2f-Src:
	-$(RM) ./Core/Src/Controllers.cyclo ./Core/Src/Controllers.d ./Core/Src/Controllers.o ./Core/Src/Controllers.su ./Core/Src/CyberGear.cyclo ./Core/Src/CyberGear.d ./Core/Src/CyberGear.o ./Core/Src/CyberGear.su ./Core/Src/CyberMove.cyclo ./Core/Src/CyberMove.d ./Core/Src/CyberMove.o ./Core/Src/CyberMove.su ./Core/Src/DDSM115.cyclo ./Core/Src/DDSM115.d ./Core/Src/DDSM115.o ./Core/Src/DDSM115.su ./Core/Src/DDSMove.cyclo ./Core/Src/DDSMove.d ./Core/Src/DDSMove.o ./Core/Src/DDSMove.su ./Core/Src/Trash.cyclo ./Core/Src/Trash.d ./Core/Src/Trash.o ./Core/Src/Trash.su ./Core/Src/kalman.cyclo ./Core/Src/kalman.d ./Core/Src/kalman.o ./Core/Src/kalman.su ./Core/Src/main.cyclo ./Core/Src/main.d ./Core/Src/main.o ./Core/Src/main.su ./Core/Src/mpu6050.cyclo ./Core/Src/mpu6050.d ./Core/Src/mpu6050.o ./Core/Src/mpu6050.su ./Core/Src/stm32f4xx_hal_msp.cyclo ./Core/Src/stm32f4xx_hal_msp.d ./Core/Src/stm32f4xx_hal_msp.o ./Core/Src/stm32f4xx_hal_msp.su ./Core/Src/stm32f4xx_it.cyclo ./Core/Src/stm32f4xx_it.d ./Core/Src/stm32f4xx_it.o ./Core/Src/stm32f4xx_it.su ./Core/Src/syscalls.cyclo ./Core/Src/syscalls.d ./Core/Src/syscalls.o ./Core/Src/syscalls.su ./Core/Src/sysmem.cyclo ./Core/Src/sysmem.d ./Core/Src/sysmem.o ./Core/Src/sysmem.su ./Core/Src/system_stm32f4xx.cyclo ./Core/Src/system_stm32f4xx.d ./Core/Src/system_stm32f4xx.o ./Core/Src/system_stm32f4xx.su

.PHONY: clean-Core-2f-Src

