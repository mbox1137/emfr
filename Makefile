.PHONY:	all clean

GDB	= -g

CC	= gcc
#clang
CFLAGS	= -m32 -fverbose-asm -Wunused-variable -Werror=return-type -Werror=dangling-else -Werror=pointer-sign -Werror=unused-but-set-variable -pthread
#CFLAGS	= -m32 -fverbose-asm -msse3 -nostdlib
# -std=c99

CXX	= $(CC)

ASM	= nasm
AFLAGS	= -f elf32

AS	= as
ASFLAGS	= -c --32
#AS	= $(CC)
#ASFLAGS	= -m32 -c

LD	= $(CXX)
LDFLAGS	= -no-pie -m32 $(GDB) -lpthread -pthread
#LD	= ld
#LDFLAGS	= -m elf_i386
#LDFLAGS	= -m elf_i386_fbsd	#-m32

EXEC	= main
OBJS	= main.o
SCDRY	= main.s

.SECONDARY: $(SCDRY)
ALL : $(EXEC)

main : main.o Makefile
	$(LD) $(LDFLAGS) -o main main.o

#$(EXEC) : $(OBJS) Makefile
#	$(LD) $(LDFLAGS) $(GDB) -o $@ $(OBJS)

%.o: %.c
#	$(CC) $(CFLAGS) $(GDB) -c -o $@ $<

%.o: %.cpp
#	$(CC) $(CFLAGS) $(GDB) -c -o $@ $<

%.s: %.c Makefile
	$(CC) $(CFLAGS) $(GDB) -S -o $@ $<

%.s: %.cpp Makefile
	$(CXX) $(CFLAGS) $(GDB) -S -o $@ $<

%.o: %.asm Makefile
	$(ASM) $(AFLAGS) $(GDB) -o $@ -l $@.lst $<

%.o: %.s Makefile
	$(AS) $(ASFLAGS) $(GDB) -a=$@.lst -o $@ $<

clean:
	rm -f $(EXEC) $(OBJS) $(SCDRY) *~ *.lst .[^.]*~ *.dat .[^.]*.bak *.bak *.log

