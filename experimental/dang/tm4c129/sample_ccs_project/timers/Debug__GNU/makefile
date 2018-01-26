################################################################################
# Automatically-generated file. Do not edit!
################################################################################

CG_TOOL_ROOT := /home/dang/ti/ccsv7/tools/compiler/gcc-arm-none-eabi-4_9-2015q3

GEN_OPTS__FLAG := 
GEN_CMDS__FLAG := 

ORDERED_OBJS += \
"./linux_time.o" \
"./startup_ccs.o" \
"./timers.o" \
"./utils/uartstdio.o" \
"../timers_ccs.cmd" \
$(GEN_CMDS__FLAG) \
-Wl,--start-group \
-lc \
-lgcc \
-lnosys \
-ldriver \
-Wl,--end-group \

-include ../makefile.init

RM := rm -rf
RMDIR := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir_vars.mk
-include utils/subdir_vars.mk
-include subdir_rules.mk
-include utils/subdir_rules.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C55_DEPS)),)
-include $(C55_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(S67_DEPS)),)
-include $(S67_DEPS)
endif
ifneq ($(strip $(S62_DEPS)),)
-include $(S62_DEPS)
endif
ifneq ($(strip $(S_DEPS)),)
-include $(S_DEPS)
endif
ifneq ($(strip $(OPT_DEPS)),)
-include $(OPT_DEPS)
endif
ifneq ($(strip $(C??_DEPS)),)
-include $(C??_DEPS)
endif
ifneq ($(strip $(ASM_UPPER_DEPS)),)
-include $(ASM_UPPER_DEPS)
endif
ifneq ($(strip $(S??_DEPS)),)
-include $(S??_DEPS)
endif
ifneq ($(strip $(C64_DEPS)),)
-include $(C64_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(S64_DEPS)),)
-include $(S64_DEPS)
endif
ifneq ($(strip $(INO_DEPS)),)
-include $(INO_DEPS)
endif
ifneq ($(strip $(CLA_DEPS)),)
-include $(CLA_DEPS)
endif
ifneq ($(strip $(S55_DEPS)),)
-include $(S55_DEPS)
endif
ifneq ($(strip $(SV7A_DEPS)),)
-include $(SV7A_DEPS)
endif
ifneq ($(strip $(C62_DEPS)),)
-include $(C62_DEPS)
endif
ifneq ($(strip $(C67_DEPS)),)
-include $(C67_DEPS)
endif
ifneq ($(strip $(PDE_DEPS)),)
-include $(PDE_DEPS)
endif
ifneq ($(strip $(K_DEPS)),)
-include $(K_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C43_DEPS)),)
-include $(C43_DEPS)
endif
ifneq ($(strip $(S43_DEPS)),)
-include $(S43_DEPS)
endif
ifneq ($(strip $(ASM_DEPS)),)
-include $(ASM_DEPS)
endif
ifneq ($(strip $(S_UPPER_DEPS)),)
-include $(S_UPPER_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(SA_DEPS)),)
-include $(SA_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 
EXE_OUTPUTS += \
timers.out \

EXE_OUTPUTS__QUOTED += \
"timers.out" \

BIN_OUTPUTS += \
timers.bin \

BIN_OUTPUTS__QUOTED += \
"timers.bin" \


# All Target
all: $(OBJS) $(CMD_SRCS) $(GEN_CMDS)
	@$(MAKE) --no-print-directory -Onone "timers.out"

# Tool invocations
timers.out: $(OBJS) $(CMD_SRCS) $(GEN_CMDS)
	@echo 'Building target: $@'
	@echo 'Invoking: GNU Linker'
	"/home/dang/ti/ccsv7/tools/compiler/gcc-arm-none-eabi-4_9-2015q3/bin/arm-none-eabi-gcc" -DTARGET_IS_TM4C129_RA0 -DPART_TM4C129ENCPDT -ffunction-sections -fdata-sections -g -gdwarf-3 -gstrict-dwarf -Wall -MD -std=c99 -march=armv7e-m -mthumb -mfloat-abi=hard -mfpu=fpv4-sp-d16 -Wl,-Map,"timers.map" -L"/home/dang/ti/TivaWare_C_Series-2.1.3.156/driverlib/gcc" -o"timers.out" $(ORDERED_OBJS)
	@echo 'Finished building target: $@'
	@echo ' '
	@$(MAKE) --no-print-directory post-build

timers.bin: $(EXE_OUTPUTS)
	@echo 'Invoking: GNU Objcopy Utility'
	"/home/dang/ti/ccsv7/tools/compiler/gcc-arm-none-eabi-4_9-2015q3/bin/arm-none-eabi-objcopy"  $(EXE_OUTPUTS__QUOTED) "timers.bin"
	@echo 'Finished building: $@'
	@echo ' '
	@$(MAKE) --no-print-directory post-build

# Other Targets
clean:
	-$(RM) $(BIN_OUTPUTS__QUOTED)$(EXE_OUTPUTS__QUOTED)
	-$(RM) "linux_time.o" "startup_ccs.o" "timers.o" "utils/uartstdio.o" 
	-$(RM) "linux_time.d" "startup_ccs.d" "timers.d" "utils/uartstdio.d" 
	-@echo 'Finished clean'
	-@echo ' '

post-build:
	-"/home/dang/ti/ccsv7/utils/tiobj2bin/tiobj2bin" "timers.out" "timers.bin" "/home/dang/ti/ccsv7/tools/compiler/gcc-arm-none-eabi-4_9-2015q3/bin/armofd" "/home/dang/ti/ccsv7/tools/compiler/gcc-arm-none-eabi-4_9-2015q3/bin/armhex" "/home/dang/ti/ccsv7/utils/tiobj2bin/mkhex4bin"
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets

