# Default target
all: strand_soup
TARGET = main.exe
TARGET_UTILITIES = utilities.exe
TARGET_STRAND_SOUP = strand_soup.exe



# Compilateur
CC = g++
# Options de compilation (can add -g to debug)
CFLAGS  = -Wall -Wextra -std=c++17 -g -I /Users/LPI-USER/miniconda3/envs/viennarna/include
LDFLAGS = -L /Users/LPI-USER/miniconda3/envs/viennarna/lib -lRNA -Wl,-rpath,/Users/LPI-USER/miniconda3/envs/viennarna/lib



# Fichiers source
SRCS = main.cpp global_variables.cpp # Ajoutez ici tous les fichiers source de votre projet
SRCS_UTILITIES = utilities.cpp global_variables.cpp
SRCS_STRAND_SOUP = strand_soup.cpp global_variables.cpp utilities.cpp  RNAEnergyEvaluator.cpp
# SRCS_VALIDATE_COFOLD = validate_cofold.cpp global_variables.cpp utilities.cpp  RNAEnergyEvaluator.cpp

# TARGET_VALIDATE_COFOLD = validate_cofold.exe





# Règle de construction du programme
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o $(TARGET)



utilities: $(SRCS_UTILITIES)
	$(CC) $(CFLAGS) -DUTILITIES_TEST $(SRCS_UTILITIES) -o $(TARGET_UTILITIES)

strand_soup: $(SRCS_STRAND_SOUP)
	$(CC) $(CFLAGS) -DSTRAND_SOUP_TEST $(SRCS_STRAND_SOUP) -o $(TARGET_STRAND_SOUP) $(LDFLAGS)

# validate_cofold: $(SRCS_VALIDATE_COFOLD)
#	$(CC) $(CFLAGS) $(SRCS_VALIDATE_COFOLD) -o $(TARGET_VALIDATE_COFOLD) $(LDFLAGS) 


# Règle pour la construction de tous les fichiers
# all: $(TARGET) $(TARGET_UTILITIES) $(TARGET_STRAND_SOUP) $(TARGET_VALIDATE_COFOLD)

all: $(TARGET_STRAND_SOUP)




# Règle pour nettoyer les fichiers objets et l'exécutable
#clean:
#	rm -f $(TARGET)  $(TARGET_UTILITIES) $(TARGET_STRAND_SOUP) $(TARGET_VALIDATE_COFOLD)
clean:
	rm -f $(TARGET_STRAND_SOUP) $(TARGET)  $(TARGET_UTILITIES)
