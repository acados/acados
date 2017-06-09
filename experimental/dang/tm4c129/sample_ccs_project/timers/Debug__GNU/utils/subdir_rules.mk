################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Each subdirectory must supply rules for building sources it contributes
utils/uartstdio.o: /home/dang/ti/TivaWare_C_Series-2.1.3.156/utils/uartstdio.c $(GEN_OPTS) | $(GEN_HDRS)
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Compiler'
	"/home/dang/ti/ccsv7/tools/compiler/gcc-arm-none-eabi-4_9-2015q3/bin/arm-none-eabi-gcc" -c -mcpu=cortex-m4 -march=armv7e-m -mthumb -mfloat-abi=hard -mfpu=fpv4-sp-d16 -DTARGET_IS_TM4C129_RA0 -DPART_TM4C129ENCPDT -I"/home/dang/ti/TivaWare_C_Series-2.1.3.156" -I"/home/dang/ccs/workspace_v7/timers" -I"/home/dang/ti/ccsv7/tools/compiler/gcc-arm-none-eabi-4_9-2015q3/arm-none-eabi/include" -ffunction-sections -fdata-sections -g -gdwarf-3 -gstrict-dwarf -Wall -MD -std=c99 -c -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o"$@" $(GEN_OPTS__FLAG) "$(shell echo $<)"
	@echo 'Finished building: $<'
	@echo ' '


