CFLAGS = -g -O0 -std=c89 -Wpedantic
CPPFLAGS := -Iinclude -MMD -MP #adds include directory of header files.
LDLIBS   := -lm

SRC_DIR = src
BIN_DIR = bin
OBJ_DIR = obj

EXE :=$(BIN_DIR)/cp_kern
SRC := $(wildcard $(SRC_DIR)/*.c)
OBJ := $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
	

.PHONY: all run clean #so it doesn't create a file or folder with these names 

all: 	$(EXE)

run: 	$(EXE)
		$(EXE)

$(EXE): $(OBJ) | $(BIN_DIR)
	$(CC) $(CFLAGS)  $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

clean:	
	@$(RM) -r $(BIN_DIR) $(OBJ_DIR) 
#The @ disables the echoing of the command

-include $(OBJ:.o=.d) # The dash silences errors when files don't exist (yet)
