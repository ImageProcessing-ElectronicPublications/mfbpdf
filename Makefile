PROJECT = mfbpdf
CC = gcc
CFLAGS = -Wall -Isrc
LDFLAGS =  -ltiff -ljpeg -lz -lm
RM = rm -f

all: $(PROJECT)

$(PROJECT): src/$(PROJECT).c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

clean:
	$(RM) $(PROJECT)
