newicktree: seqMain.o seqUtil.o Newickform.o
	gcc -o newicktree seqMain.o seqUtil.o Newickform.o
seqMain.o: seqMain.c
	gcc -g -c seqMain.c
seqUtil.o: seqUtil.c
	gcc -g -c seqUtil.c
Newickform.o: Newickform.c
	gcc -g -c Newickform.c

clean:
	-rm *.o newicktree
