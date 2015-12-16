VERSION=1.0.0
MINOR_VER=20140314
CC=gcc
ifdef DEBUG
CFLAGS=-g3 -W -Wall -O0 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -mpopcnt -mssse3
else
CFLAGS=-W -Wall -O4 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -mpopcnt -mssse3
endif

GLIBS=-lm -lpthread
#GENERIC_SRC=Makefile mem_share.h string.h file_reader.h file_reader.c bitvec.h hashset.h sort.h list.h dna.h thread.h timer.h ksw.h ksw.c kswx.h
GENERIC_SRC=mem_share.h string.h file_reader.h file_reader.c bitvec.h hashset.h sort.h list.h dna.h thread.h timer.h ksw.h ksw.c kswx.h

ifeq (0, ${MAKELEVEL})
UNIQ_ID=$(shell date +"%s")
endif
#all: pairaln wtcyc wtmer wtzmo wtobt wtclp wtext wtlay wtcns wtdif wtcorr wtjnt
all: pairaln wtpre wtcyc wtmer wtzmo wtobt wtclp wtext wtgbo wtlay wtcns wtmsa

pairaln: $(GENERIC_SRC) pairaln.c
	$(CC) $(CFLAGS) -o pairaln file_reader.c ksw.c pairaln.c $(GLIBS)

wtpre: $(GENERIC_SRC) wtpre.c
	$(CC) $(CFLAGS) -o wtpre file_reader.c wtpre.c $(GLIBS)

wtcyc: $(GENERIC_SRC) wtcyc.c
	$(CC) $(CFLAGS) -o wtcyc file_reader.c ksw.c wtcyc.c $(GLIBS)

wtmer: $(GENERIC_SRC) wtmer.c
	$(CC) $(CFLAGS) -o wtmer file_reader.c wtmer.c $(GLIBS)

wtzmo: $(GENERIC_SRC) hzm_aln.h wtzmo.c
	$(CC) $(CFLAGS) -o wtzmo file_reader.c ksw.c wtzmo.c $(GLIBS)

wtobt: $(GENERIC_SRC) wtobt.c
	$(CC) $(CFLAGS) -o wtobt file_reader.c wtobt.c $(GLIBS)

wtclp: $(GENERIC_SRC) wtclp.c
	$(CC) $(CFLAGS) -o wtclp file_reader.c wtclp.c $(GLIBS)

wtext: $(GENERIC_SRC) wtext.c
	$(CC) $(CFLAGS) -o wtext file_reader.c ksw.c wtext.c $(GLIBS)

wtgbo: $(GENERIC_SRC) wtgbo.c wtlay.h
	$(CC) $(CFLAGS) -o wtgbo file_reader.c ksw.c wtgbo.c $(GLIBS)

wtlay: $(GENERIC_SRC) wtlay.h wtlay.c
	$(CC) $(CFLAGS) -o wtlay file_reader.c wtlay.c $(GLIBS)

wtcns: $(GENERIC_SRC) wtcns.c hzm_aln.h bit2vec.h queue.h dagcns.h
	$(CC) $(CFLAGS) -o wtcns file_reader.c ksw.c wtcns.c $(GLIBS)

wtmsa: $(GENERIC_SRC) wtmsa.c hzm_aln.h bit2vec.h bit2vec.h pomsa.h
	$(CC) $(CFLAGS) -o wtmsa file_reader.c ksw.c wtmsa.c $(GLIBS)

wtdif: $(GENERIC_SRC) wtdif.c
	$(CC) $(CFLAGS) -o wtdif file_reader.c ksw.c wtdif.c $(GLIBS)

wtcorr: $(GENERIC_SRC) wtcorr.c bitsvec.h counting_bloom_filter.h
	$(CC) $(CFLAGS) -o wtcorr wtcorr.c file_reader.c $(GLIBS)

wtjnt: $(GENERIC_SRC) wtjnt.c hzm_aln.h dagcns.h
	$(CC) $(CFLAGS) -o wtjnt wtjnt.c file_reader.c ksw.c $(GLIBS)

clean:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out pairaln wtpre wtcyc wtmer wtzmo wtobt wtclp wtext wtlay wtcns wtmsa wtdif wtcorr wtjnt

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out
