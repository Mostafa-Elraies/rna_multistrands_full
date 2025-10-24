# Nom du programme
TARGET = main.exe
TARGET_NUSSINOV = nussinov.exe
TARGET_UTILITIES = utilities.exe
TARGET_STRAND_SOUP = strand_soup.exe



# Compilateur
CC = g++
# Options de compilation (can add -g to debug)
CFLAGS = -Wall -Wextra -std=c++17 -g -I /opt/anaconda3/envs/INF551_env/include #-I ./../../ViennaRNA-2.7.0/include # -I/opt/anaconda3/envs/INF551_env/include

LDFLAGS =  -L /opt/anaconda3/envs/INF551_env/lib -lRNA -Wl,-rpath,/opt/anaconda3/envs/INF551_env/lib # -L ./../../ViennaRNA-2.7.0/lib -lRNA #/opt/anaconda3/envs/INF551_env/lib/libRNA.a -Wl,-rpath,/opt/anaconda3/envs/INF551_env/lib




# Fichiers source
SRCS = main.cpp global_variables.cpp # Ajoutez ici tous les fichiers source de votre projet
SRCS_NUSSINOV = nussinov.cpp global_variables.cpp utilities.cpp
SRCS_UTILITIES = utilities.cpp global_variables.cpp
SRCS_STRAND_SOUP = strand_soup.cpp global_variables.cpp utilities.cpp nussinov.cpp




# Règle de construction du programme
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o $(TARGET)

nussinov: $(SRCS_NUSSINOV)
	$(CC) $(CFLAGS) -DNUSSINOV_TEST $(SRCS_NUSSINOV) -o $(TARGET_NUSSINOV) $(LDFLAGS)

utilities: $(SRCS_UTILITIES)
	$(CC) $(CFLAGS) -DUTILITIES_TEST $(SRCS_UTILITIES) -o $(TARGET_UTILITIES)

strand_soup: $(SRCS_STRAND_SOUP)
	$(CC) $(CFLAGS) -DSTRAND_SOUP_TEST $(SRCS_STRAND_SOUP) -o $(TARGET_STRAND_SOUP) $(LDFLAGS)

# Règle pour la construction de tous les fichiers
all: $(TARGET) $(TARGET_NUSSINOV) $(TARGET_UTILITIES) $(TARGET_STRAND_SOUP)




# Règle pour nettoyer les fichiers objets et l'exécutable
clean:
	rm -f $(TARGET) $(TARGET_NUSSINOV) $(TARGET_UTILITIES) $(TARGET_STRAND_SOUP)
