CFLAGS=-O3 -Wuninitialized -Wall -fomit-frame-pointer -ffast-math
LDLIBS=-lm

BINARIES=pitchotron-bin pix-y4m
all: $(BINARIES)

clean:
	rm -f $(BINARIES)

install: pitchotron-bin pix-y4m pitchotron
	install pitchotron-bin /usr/local/bin
	install pix-y4m /usr/local/bin
	install pitchotron /usr/local/bin

pitchotron-bin: pitchotron-bin.c
	gcc -o $@ $(CFLAGS) -o $@ $< $(LDLIBS)  

pix-y4m-notable: pix-y4m.c constants.h
	gcc -o $@ $(CFLAGS) $(CPPFLAGS) -DNOTABLE $< $(LDLIBS) 

pix-y4m.o: pix-y4m.c constants.h
	gcc -c -o $@ $(CFLAGS) $(CPPFLAGS) $<

pix-y4m: pix-y4m.o hyperluma_table.o
	gcc -o $@ $^ $(LDLIBS)
	./$@ -dump_raw_table > raw_table_dumped
	./$@ -calc_raw_table > raw_table_calculated
	if diff raw_table_calculated raw_table_dumped;\
		then rm raw_table_calculated \
		raw_table_dumped; \
		else rm pix-y4m; false; fi

hyperluma_table.o: hyperluma_table.c
	gcc -c -o $@ $<

hyperluma_table.c: pix-y4m-notable
	./$< -calc_compressed_table > $@


	
