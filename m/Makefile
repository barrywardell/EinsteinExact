all: EinsteinExact.out
%.out: %.m Metrics/*.m Metrics/*/*.m
	../../Kranc/Bin/kranc $*.m 2>&1 | tee $*.out
clean:
	rm -f EinsteinExact.out
.PHONY: all clean
