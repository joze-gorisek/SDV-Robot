# Enable compile command to ease indexing with e.g. clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

# Compiler options
target_compile_options(${BUILD_UNIT_0_NAME} PRIVATE
    $<$<COMPILE_LANGUAGE:C>: ${CUBE_CMAKE_C_FLAGS}>
    $<$<COMPILE_LANGUAGE:CXX>: ${CUBE_CMAKE_CXX_FLAGS}>
    $<$<COMPILE_LANGUAGE:ASM>: ${CUBE_CMAKE_ASM_FLAGS}>
)

# Linker options
target_link_options(${BUILD_UNIT_0_NAME} PRIVATE ${CUBE_CMAKE_EXE_LINKER_FLAGS})

# Add sources to executable/library
target_sources(${BUILD_UNIT_0_NAME} PRIVATE
    "Core/Src/Controllers.c"
    "Core/Src/CyberGear.c"
    "Core/Src/CyberMove.c"
    "Core/Src/DDSM115.c"
    "Core/Src/DDSMove.c"
    "Core/Src/MRF24J40.c"
    "Core/Src/Trash.c"
    "Core/Src/kalman.c"
    "Core/Src/main.c"
    "Core/Src/mpu6050.c"
    "Core/Src/stm32f4xx_hal_msp.c"
    "Core/Src/stm32f4xx_it.c"
    "Core/Src/syscalls.c"
    "Core/Src/sysmem.c"
    "Core/Src/system_stm32f4xx.c"
    "Core/Startup/startup_stm32f446retx.s"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_can.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_cortex.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_dma.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_dma_ex.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_exti.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_flash.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_flash_ex.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_flash_ramfunc.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_gpio.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_i2c.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_i2c_ex.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_pwr.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_pwr_ex.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_rcc.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_rcc_ex.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_spi.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_tim.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_tim_ex.c"
    "Drivers/STM32F4xx_HAL_Driver/Src/stm32f4xx_hal_uart.c"
    "SH2 Sensorhub/Src/demo_app.c"
    "SH2 Sensorhub/Src/euler.c"
    "SH2 Sensorhub/Src/sh2.c"
    "SH2 Sensorhub/Src/sh2_SensorValue.c"
    "SH2 Sensorhub/Src/sh2_util.c"
    "SH2 Sensorhub/Src/shtp.c"
)

target_include_directories(${BUILD_UNIT_0_NAME} PRIVATE
    "Core/Inc"
    "Drivers/STM32F4xx_HAL_Driver/Inc"
    "Drivers/STM32F4xx_HAL_Driver/Inc/Legacy"
    "Drivers/CMSIS/Device/ST/STM32F4xx/Include"
    "Drivers/CMSIS/Include"
    "SH2 Sensorhub/Inc"
)

configure_file("${CMAKE_SOURCE_DIR}/STM32F446RETX_FLASH.ld" "${CMAKE_BINARY_DIR}" COPYONLY)

set_target_properties(${BUILD_UNIT_0_NAME} PROPERTIES LINK_DEPENDS "STM32F446RETX_FLASH.ld")

