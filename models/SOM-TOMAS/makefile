# ============================================================================================
#                             THIS IS THE MAKEFILE FOR BOX.EXE
# ============================================================================================

# COMPILERS AND FLAGS:
# ============================================================================================
COMP1 = gfortran
COMP2 = g++

FLAGS1 = -O3 -c
FLAGS2 = -std=c++11 -c
FLAGS3 = -fbounds-check -Wall -lstdc++ -g

# PATHS:
# ============================================================================================
PATH0 = objects/
PATH1 = source_ftn/
PATH2 = source_ftn/subroutines/
PATH3 = source_ftn/functions/
PATH4 = source_ftn/modules/
PATH5 = source_cpp/

# OBJECT FILES:
# ============================================================================================
OBJS =  \
$(PATH0)mod.som.o           \
$(PATH0)mod.tomas.o         \
$(PATH0)mod.invars.o        \
$(PATH0)mod.ofiles.o        \
$(PATH0)mod.btab.o          \
$(PATH0)mod.main.o          \
$(PATH0)main.o              \
$(PATH0)r_init_pphase.o     \
$(PATH0)r_step_wallloss.o   \
$(PATH0)r_read_inputs.o     \
$(PATH0)r_init_outfiles.o   \
$(PATH0)r_step_gaschem.o    \
$(PATH0)r_step_hetchem.o    \
$(PATH0)r_step_autoox.o     \
$(PATH0)r_step_aerochem.o   \
$(PATH0)r_solve_cond.o      \
$(PATH0)r_step_coag.o       \
$(PATH0)r_step_nucl.o       \
$(PATH0)r_step_soacond.o    \
$(PATH0)r_write_outputs.o   \
$(PATH0)r_fix_mn.o          \
$(PATH0)r_init_gphase.o     \
$(PATH0)f_get_len.o         \
$(PATH0)f_calc_pwl.o        \
$(PATH0)f_calc_dmdt.o       \
$(PATH0)f_calc_aedens.o     \
$(PATH0)f_get_max.o         \
$(PATH0)f_calc_db.o         \
$(PATH0)f_calc_kpar.o       \
$(PATH0)f_calc_qbulk.o      \
$(PATH0)f_calc_jnuc.o       \
$(PATH0)gas_het_chem.o      \

# MAKE CONTROL:
# ============================================================================================
all: sub1 box.exe

sub1:
	@for i in $(PATH1) $(PATH2) $(PATH3) $(PATH4); do \
	 cp -u $$i/*.f90 $(PATH0) 2>/dev/null || :; cp -u $$i/*.f $(PATH0) 2>/dev/null || :; \
	 done

# EXECUTABLE:
# ============================================================================================
box.exe : $(OBJS) makefile
	@echo "making... $@"
	@$(COMP1) $(FLAGS3) -o $@ $(OBJS)

# OBJECTS:
# ============================================================================================
# FORTRAN CODE:
$(PATH0)%.o : $(PATH0)%.f90 makefile
	@echo "making... $@"
	@$(COMP1) -I$(PATH0) -J$(PATH0) $(FLAGS1) $< -o $@

$(PATH0)%.o : $(PATH0)%.f makefile
	@echo "making... $@"
	@$(COMP1) -I$(PATH0) -J$(PATH0) $(FLAGS1) $< -o $@

# GAS-PHASE CHEM. FROM C++:
$(PATH0)gas_het_chem.o : $(PATH5)main/main.cpp makefile
	@echo "making... $@"
	@$(COMP2) $(FLAGS2) $< -o $@

# ULTILITIES:
# ============================================================================================
clean :
	rm objects/*.o objects/*.mod box.exe

clear :
	rm outputs/*

