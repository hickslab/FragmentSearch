CC=g++
CFLAGS=-pthread
OUT=fragmentsearch
OBJS=Configuration.o Database.o FragmentSearch.o main.o Protein.o Results.o

%.o: %.cpp
	$(CC) $(CFLAGS) $(CLIBS) -c $< -o $@

$(OUT): $(OBJS)
	$(CC) $(CFLAGS) $(CLIBS) $(OBJS) -o $@

default: $(OUT)

clean:
	-rm -f 
	-rm -f $(OUT)
	-rm -f $(OBJS)
